"""Tests 4 and 5: membrane charge balance and RK4 accuracy against DOP853."""
import numpy as np
import pytest

from hh_energy import membrane


@pytest.mark.parametrize("T_C", [6.3, 25.0])
def test_4_membrane_charge_balance(amplitudes, T_C, record_property):
    """∫ I_stim dt = ∫ I_ion dt + C_m (V_end − V_start) within 0.1 % once the spike is over."""
    tr = membrane.simulate_one(T_C, amplitudes["amp2x"])
    assert tr.back_to_rest
    q_stim = tr.stim_charge
    q_ion = np.trapezoid(tr.I_Na + tr.I_K + tr.I_L, tr.t)
    q_cap = 1.0 * (tr.V[-1] - tr.V[0])
    record_property("measured", f"relative imbalance = {abs(q_stim - (q_ion + q_cap)) / q_stim:.2e} (limit 1e-3)")
    assert abs(q_stim - (q_ion + q_cap)) / q_stim < 1e-3


@pytest.mark.parametrize("T_C", [6.3, 25.0])
def test_5_rk4_vs_dop853(amplitudes, T_C, record_property):
    """Fixed-step RK4 (Δt = 5 µs) differs from DOP853 (rtol 1e−10) by < 0.01 mV everywhere."""
    tr = membrane.simulate_one(T_C, amplitudes["amp2x"])
    ref = membrane.simulate_reference(T_C, amplitudes["amp2x"], tr.t[-1], tr.dt)
    record_property("measured", f"max |dV| = {np.max(np.abs(tr.V - ref.V)):.2e} mV (limit 0.01 mV)")
    assert np.max(np.abs(tr.V - ref.V)) < 0.01
