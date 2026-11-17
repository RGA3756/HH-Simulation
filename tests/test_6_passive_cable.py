"""Test 6: passive cable against the analytic steady state (checks the units conversion)."""
import math

import numpy as np

from hh_energy import cable, membrane
from hh_energy import params as P


def test_6_coupling_constant(record_property):
    """D = a/(2R_i) × 1000 ≈ 0.336 µA/cm² per (mV/cm²)."""
    record_property("measured", f"D = {P.D_COUPLING:.6f} uA/cm2 per (mV/cm2)")
    assert abs(P.D_COUPLING - 0.336) < 1e-3


def test_6_passive_cable_matches_analytic(record_property):
    """Frozen gates + constant current at x = 0: V(x) − V_rest ∝ cosh((L − x)/λ) within 1 %.

    λ and the input resistance are computed here from textbook cable theory in
    ohms and centimetres, independently of D_COUPLING, so a wrong units factor
    in the simulator would make this test fail.
    """
    dx, L, a, R_i = P.DX_CABLE_CM, P.LENGTH_CM, P.RADIUS_CM, P.R_I_OHM_CM
    run = cable.simulate(P.T_REF_C, freeze_gates=True, stim_nodes=[0], stim_amp=P.PASSIVE_TEST_I_UA,
                         stim_start=0.0, stim_duration=P.PASSIVE_TEST_RUN_MS,
                         t_end=P.PASSIVE_TEST_RUN_MS)
    rest = membrane.resting_state()
    G_rest = P.HH.g_Na * rest.m**3 * rest.h + P.HH.g_K * rest.n**4 + P.HH.g_L   # mS/cm²

    # Steady state reached: last two snapshots (0.05 ms apart) agree to 1e-9 mV
    assert np.max(np.abs(run.V_snap[-1] - run.V_snap[-2])) < 1e-9

    R_m = 1000.0 / G_rest                      # Ω·cm²   (1 mS = 1e−3 S)
    lam = math.sqrt(R_m * a / (2.0 * R_i))     # cm
    dV = run.V_snap[-1] - rest.V
    shape_sim = dV / dV[0]
    shape_ana = np.cosh((L - run.x) / lam) / math.cosh(L / lam)
    shape_err = np.max(np.abs(shape_sim / shape_ana - 1.0))
    assert shape_err < P.PASSIVE_TEST_TOL

    # Absolute amplitude: injected current × input resistance of a sealed cable
    I_total_A = P.PASSIVE_TEST_I_UA * 1e-6 * (2.0 * math.pi * a) * (dx / 2.0)   # node 0 is a half cell
    R_in_ohm = R_i / (math.pi * a * a) * lam / math.tanh(L / lam)
    V0_mV = I_total_A * R_in_ohm * 1000.0
    record_property("measured", f"lambda = {lam:.4f} cm; max shape error = {shape_err:.1e}; "
                    f"V(0) - V_rest = {dV[0]:.5f} mV vs analytic {V0_mV:.5f} mV (limit 1%)")
    assert abs(dV[0] / V0_mV - 1.0) < P.PASSIVE_TEST_TOL
