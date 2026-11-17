"""
Hodgkin–Huxley cable (unmyelinated axon), sealed ends.

    C_m ∂V/∂t = D ∂²V/∂x² − (I_Na + I_K + I_L) + I_stim(x, t)
    D = (a / 2R_i) × 1000          (µA/cm² per mV/cm², see params.D_COUPLING)

Units: V mV, t ms, x and a cm, R_i Ω·cm, I µA/cm², g mS/cm², C_m µF/cm².

Grid: vertex-centred nodes x_i = iΔx, i = 0..N (N = L/Δx). Sealed ends use a
mirrored ghost node (V₋₁ = V₁, V_{N+1} = V_{N−1}), so the end Laplacian is
2(V₁ − V₀)/Δx².

Each time step:
 1. gates by exponential Euler at the old voltage Vⁿ:
        z^{n+1} = z∞(Vⁿ) + (zⁿ − z∞(Vⁿ)) exp(−Δt/τ_z(Vⁿ, T))
 2. voltage by backward Euler with the new gates (the ionic current is
    linear in V once the gates are fixed):
        (C_m/Δt + G) V^{n+1} − D Lap(V^{n+1}) = (C_m/Δt) Vⁿ + Σ g_k E_k + I_stim
    solved with scipy.linalg.solve_banded (tridiagonal).

Why backward Euler: an explicit (forward Euler) voltage step is only stable
if D Δt / (C_m Δx²) ≤ 1/2. Here that number is ≈ 16.8 (see
explicit_stability_number), so an explicit step would blow up.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded

from . import params as P
from .membrane import ionic_currents, resting_state
from .rates import all_rates, phi


def explicit_stability_number(dx=P.DX_CABLE_CM, dt=P.DT_CABLE_MS, C_m=P.HH.C_m, D=P.D_COUPLING):
    """D Δt / (C_m Δx²); forward Euler needs this ≤ 0.5."""
    return D * dt / (C_m * dx * dx)


def grid(dx=P.DX_CABLE_CM, L=P.LENGTH_CM):
    """Node positions x_i = iΔx (cm), i = 0..N; raises if L/Δx is not an integer."""
    N = int(round(L / dx))
    if abs(N * dx - L) > 1e-9:
        raise ValueError(f"L = {L} cm is not a whole number of Δx = {dx} cm")
    return np.arange(N + 1) * dx


def node_index(x_cm, dx):
    """Index of the node at x (must be an exact node)."""
    i = int(round(x_cm / dx))
    if abs(i * dx - x_cm) > 1e-9:
        raise ValueError(f"x = {x_cm} cm is not a node for Δx = {dx} cm")
    return i


def banded_laplacian_matrix(n_nodes, dx, D):
    """Banded storage of −D·Lap (sealed ends), for solve_banded((1, 1), ...).

    Row i:  −k V_{i−1} + 2k V_i − k V_{i+1},  k = D/Δx²;
    row 0:  2k V_0 − 2k V_1;   row N: −2k V_{N−1} + 2k V_N.
    """
    k = D / (dx * dx)
    ab = np.zeros((3, n_nodes))
    ab[0, 1:] = -k            # super-diagonal A[i, i+1]
    ab[0, 1] = -2.0 * k       # A[0, 1] (mirror at x = 0)
    ab[2, :-1] = -k           # sub-diagonal A[i+1, i]
    ab[2, -2] = -2.0 * k      # A[N, N−1] (mirror at x = L)
    ab[1, :] = 2.0 * k        # diagonal part of −D·Lap
    return ab


@dataclass
class CableRun:
    """Result of one cable simulation."""
    T_C: float
    dx: float
    dt: float
    x: np.ndarray                    # node positions (cm)
    t: np.ndarray                    # every time step (ms)
    record_x: tuple                  # positions with full-resolution traces
    V_rec: dict                      # x -> V(t)
    I_Na_rec: dict                   # x -> I_Na(t)  (negative = Na⁺ entering)
    I_K_rec: dict
    I_L_rec: dict
    t_snap: np.ndarray               # snapshot times (every SNAPSHOT_EVERY_MS)
    V_snap: np.ndarray               # (n_snap, n_nodes)
    rest: object
    back_to_rest: bool               # x = 2.5 cm returned within 0.5 mV of rest after its AHP
    stop_reason: str
    info: dict = field(default_factory=dict)


def simulate(T_C, dx=P.DX_CABLE_CM, dt=P.DT_CABLE_MS, p=P.HH, D=P.D_COUPLING,
             stim_amp=P.CABLE_STIM_AMP_UA, stim_x_max=P.CABLE_STIM_X_MAX_CM,
             stim_start=P.STIM_START_MS, stim_duration=P.STIM_DURATION_MS,
             stim_nodes=None, t_end=None, freeze_gates=False,
             record_x=P.RECORD_X_CM, snapshot_every=P.SNAPSHOT_EVERY_MS,
             L=P.LENGTH_CM) -> CableRun:
    """Simulate the cable at temperature T_C (°C).

    Stimulus: stim_amp (µA/cm²) on every node with x < stim_x_max (or on the
    explicit index array stim_nodes) during [stim_start, stim_start + stim_duration).
    The stimulus is on during step n → n+1 when stim_start ≤ t_n < stim_end, so it
    delivers exactly stim_amp × stim_duration of charge per node.

    t_end: fixed run length (ms). If None: run until the node at x = 2.5 cm is
    back within 0.5 mV of V_rest after its AHP (and the node at 4.0 cm has had
    time to fire), capped at T_MAX_MS. A run in which the spike never reaches
    2.5 cm stops at NO_SPIKE_STOP_MS.

    freeze_gates: keep m, h, n at rest (passive cable, used by the units test).
    """
    rest = resting_state(p)
    x = grid(dx, L)
    n_nodes = x.size
    phi_T = float(phi(T_C))

    # State
    V = np.full(n_nodes, rest.V)
    m = np.full(n_nodes, rest.m)
    h = np.full(n_nodes, rest.h)
    n = np.full(n_nodes, rest.n)

    # Stimulus
    I_stim_on = np.zeros(n_nodes)
    if stim_nodes is None:
        I_stim_on[x < stim_x_max - 1e-12] = stim_amp
    else:
        I_stim_on[np.asarray(stim_nodes)] = stim_amp
    k_on = int(round(stim_start / dt))
    k_off = int(round((stim_start + stim_duration) / dt))
    assert abs(k_on * dt - stim_start) < 1e-9 and abs(k_off * dt - stim_start - stim_duration) < 1e-9

    # Matrix pieces
    ab_lap = banded_laplacian_matrix(n_nodes, dx, D)
    ab = ab_lap.copy()
    c_dt = p.C_m / dt

    # Recording
    rec_idx = {xr: node_index(xr, dx) for xr in record_x}
    cap = int(round(P.T_MAX_MS / dt))
    n_steps = cap if t_end is None else int(round(t_end / dt))
    rec_V = {xr: np.empty(n_steps + 1) for xr in record_x}
    rec_Na = {xr: np.empty(n_steps + 1) for xr in record_x}
    rec_K = {xr: np.empty(n_steps + 1) for xr in record_x}
    rec_L = {xr: np.empty(n_steps + 1) for xr in record_x}
    snap_every = int(round(snapshot_every / dt))
    snaps = []

    def record(k):
        for xr, i in rec_idx.items():
            rec_V[xr][k] = V[i]
            I_Na, I_K, I_L = ionic_currents(V[i], m[i], h[i], n[i], p)
            rec_Na[xr][k], rec_K[xr][k], rec_L[xr][k] = I_Na, I_K, I_L
        if k % snap_every == 0:
            snaps.append(V.copy())

    # Run-length bookkeeping at x = 2.5 cm (measurement node) and x = 4.0 cm
    adaptive = t_end is None
    i_meas = node_index(P.MEASURE_X_CM, dx) if adaptive else None
    i_far = node_index(P.PROPAGATION_X_CM, dx) if adaptive else None
    arrived = peaked = in_ahp = False
    v_max_meas = -np.inf
    k_recovered = None
    far_fired = False
    stop_reason = "fixed length" if not adaptive else "cap reached"
    tol = P.REST_RETURN_TOL_MV

    record(0)
    k = 0
    while k < n_steps:
        # 1. gates by exponential Euler at the old voltage
        if not freeze_gates:
            am, bm, ah, bh, an, bn = all_rates(V)
            for z, a_, b_ in ((m, am, bm), (h, ah, bh), (n, an, bn)):
                s = a_ + b_
                z_inf = a_ / s
                z[:] = z_inf + (z - z_inf) * np.exp(-dt * phi_T * s)
        # 2. voltage by backward Euler with the new gates
        g_Na = p.g_Na * m**3 * h
        g_K = p.g_K * n**4
        G = g_Na + g_K + p.g_L
        rhs = c_dt * V + g_Na * p.E_Na + g_K * p.E_K + p.g_L * p.E_L
        if k_on <= k < k_off:
            rhs = rhs + I_stim_on
        ab[1] = ab_lap[1] + c_dt + G
        V = solve_banded((1, 1), ab, rhs, overwrite_ab=False, overwrite_b=True, check_finite=False)
        k += 1
        record(k)

        if adaptive:
            v = V[i_meas]
            if not arrived:
                arrived = v > rest.V + P.CABLE_ARRIVAL_MARGIN_MV
                if not arrived and k * dt >= P.NO_SPIKE_STOP_MS:
                    stop_reason = "no spike at x = 2.5 cm"
                    break
            if V[i_far] > P.PROPAGATION_V_MV:
                far_fired = True
            if arrived and k_recovered is None:
                if v > v_max_meas:
                    v_max_meas = v
                elif not peaked:
                    peaked = True
                if peaked and v < rest.V - tol:
                    in_ahp = True
                if in_ahp and v >= rest.V - tol:
                    k_recovered = k
            if k_recovered is not None and (
                    far_fired or (k - k_recovered) * dt >= P.CABLE_FAR_NODE_WAIT_MS):
                stop_reason = "x = 2.5 cm back within 0.5 mV of rest"
                break

    t = np.arange(k + 1) * dt
    trim = lambda d: {xr: arr[:k + 1].copy() for xr, arr in d.items()}
    return CableRun(
        T_C=float(T_C), dx=dx, dt=dt, x=x, t=t, record_x=tuple(record_x),
        V_rec=trim(rec_V), I_Na_rec=trim(rec_Na), I_K_rec=trim(rec_K), I_L_rec=trim(rec_L),
        t_snap=np.arange(len(snaps)) * snap_every * dt, V_snap=np.array(snaps), rest=rest,
        back_to_rest=k_recovered is not None, stop_reason=stop_reason,
        info={"k_recovered": k_recovered, "explicit_number": explicit_stability_number(dx, dt, p.C_m, D)})
