"""
Space-clamped Hodgkin–Huxley membrane patch.

    C_m dV/dt = −I_Na − I_K − I_L + I_stim
    dz/dt     = φ(T) [α_z(V)(1 − z) − β_z(V) z],   z ∈ {m, h, n}

Units: V mV, t ms, I µA/cm², C_m µF/cm², g mS/cm², T °C.
Sign convention: I_Na = ḡ_Na m³h (V − E_Na) is NEGATIVE when Na⁺ flows in.

Integration: classical RK4 at a fixed step. The stimulus is a list of
constant-current pieces [(t_start, t_end, amplitude), ...]. Every piece is
integrated separately and its boundaries must fall exactly on step
boundaries, so no RK4 step ever straddles a jump in the stimulus (a jump
inside a step would drop RK4 to first order).

Several simulations can run side by side ("batch"): T, the stimulus
amplitude and ḡ_Na may be arrays of the same length B. The runs do not
interact; batching only shares the Python loop.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import params as P
from .rates import all_rates, phi, steady_state


# ---------------------------------------------------------------------------
# Currents and resting state
# ---------------------------------------------------------------------------
def ionic_currents(V, m, h, n, p=P.HH, g_Na=None):
    """(I_Na, I_K, I_L) in µA/cm², outward positive. g_Na overrides p.g_Na if given."""
    gNa = p.g_Na if g_Na is None else g_Na
    I_Na = gNa * m**3 * h * (V - p.E_Na)
    I_K = p.g_K * n**4 * (V - p.E_K)
    I_L = p.g_L * (V - p.E_L)
    return I_Na, I_K, I_L


def total_steady_state_current(V, p=P.HH):
    """Total ionic current (µA/cm²) with every gate at its steady state z∞(V)."""
    return sum(ionic_currents(V, *steady_state(V), p=p))


@dataclass(frozen=True)
class RestState:
    V: float
    m: float
    h: float
    n: float
    I_Na: float   # µA/cm², negative (small resting Na⁺ influx)
    I_K: float
    I_L: float


def resting_state(p=P.HH) -> RestState:
    """V_rest from the zero of the total steady-state ionic current; gates at z∞(V_rest).

    Does not depend on temperature (z∞ does not). Asserts |V_rest − (−65)| < 0.5 mV.
    """
    lo, hi = P.V_REST_BRACKET_MV
    V = brentq(total_steady_state_current, lo, hi, args=(p,), xtol=1e-13, rtol=1e-15, maxiter=200)
    assert abs(V - P.V_REST_EXPECTED_MV) < P.V_REST_TOLERANCE_MV, f"V_rest = {V} mV"
    m, h, n = (float(z) for z in steady_state(V))
    I_Na, I_K, I_L = (float(i) for i in ionic_currents(V, m, h, n, p))
    return RestState(float(V), m, h, n, I_Na, I_K, I_L)


# ---------------------------------------------------------------------------
# Right-hand side and RK4
# ---------------------------------------------------------------------------
def rhs(y, I_stim, phi_T, g_Na, p=P.HH):
    """dy/dt for y = (V, m, h, n); each entry is a float or an array of shape (B,)."""
    V, m, h, n = y
    am, bm, ah, bh, an, bn = all_rates(V)
    I_Na, I_K, I_L = ionic_currents(V, m, h, n, p, g_Na)
    return np.array([
        (-(I_Na + I_K + I_L) + I_stim) / p.C_m,
        phi_T * (am * (1.0 - m) - bm * m),
        phi_T * (ah * (1.0 - h) - bh * h),
        phi_T * (an * (1.0 - n) - bn * n),
    ])


def rk4_step(y, dt, I_stim, phi_T, g_Na, p=P.HH):
    k1 = rhs(y, I_stim, phi_T, g_Na, p)
    k2 = rhs(y + 0.5 * dt * k1, I_stim, phi_T, g_Na, p)
    k3 = rhs(y + 0.5 * dt * k2, I_stim, phi_T, g_Na, p)
    k4 = rhs(y + dt * k3, I_stim, phi_T, g_Na, p)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def steps_in(t0, t1, dt):
    """Number of steps of size dt in [t0, t1]; raises if the boundary is off the step grid."""
    n = int(round((t1 - t0) / dt))
    if abs(n * dt - (t1 - t0)) > 1e-9 * max(1.0, abs(t1 - t0)):
        raise ValueError(f"piece [{t0}, {t1}] is not a whole number of steps of {dt} ms")
    return n


def pulse_pieces(amp, t_end, start=P.STIM_START_MS, duration=P.STIM_DURATION_MS):
    """Single rectangular pulse as constant pieces: off, on (amp), off until t_end."""
    return [(0.0, start, 0.0), (start, start + duration, amp), (start + duration, t_end, 0.0)]


# ---------------------------------------------------------------------------
# Simulation result
# ---------------------------------------------------------------------------
@dataclass
class MembraneTrace:
    """One simulation (1-D arrays over time). Units as in the module docstring."""
    t: np.ndarray
    V: np.ndarray
    m: np.ndarray
    h: np.ndarray
    n: np.ndarray
    I_Na: np.ndarray
    I_K: np.ndarray
    I_L: np.ndarray
    I_stim: np.ndarray
    T_C: float
    amp: float
    g_Na: float
    dt: float
    rest: RestState
    stim_start: float
    stim_end: float
    back_to_rest: bool   # run-length rule satisfied (V back within tol of rest after the AHP)

    @property
    def stim_charge(self):
        """∫ I_stim dt (nC/cm²), exact from the piece definition."""
        return self.amp * (self.stim_end - self.stim_start)


def _end_index(V, t, V_rest, tol=P.REST_RETURN_TOL_MV):
    """Index where a run ends, and whether the back-to-rest rule was met.

    Spiking run (max V > 0 mV): first sample after the after-hyperpolarisation
    (AHP) minimum at which V is back within `tol` of V_rest. "After the AHP"
    means V has gone more than `tol` below rest after the peak and come back.
    Run with no spike: the sample at NO_SPIKE_STOP_MS.
    Returns (index or None if the rule is not met yet, rule_met).
    """
    i_peak = int(np.argmax(V))
    if V[i_peak] <= P.SPIKE_V_MIN_MV:
        k = int(np.searchsorted(t, P.NO_SPIKE_STOP_MS - 1e-9))
        return (k if k < len(t) else None), False
    after = V[i_peak:]
    below = np.nonzero(after < V_rest - tol)[0]
    if below.size == 0:
        return None, False
    j0 = i_peak + below[0]
    back = np.nonzero(V[j0:] >= V_rest - tol)[0]
    if back.size == 0:
        return None, False
    return j0 + int(back[0]), True


def simulate(T_C, amp, p=P.HH, g_Na=None, dt=P.DT_MEMBRANE_MS, t_end=None,
             start=P.STIM_START_MS, duration=P.STIM_DURATION_MS):
    """Run one or several membrane-patch simulations with RK4.

    T_C, amp, g_Na: float or arrays (broadcast to one batch).
    t_end: fixed run length in ms. If None, each run lasts until V is back
    within 0.5 mV of V_rest after its AHP (extended in CHUNK_MS chunks,
    capped at T_MAX_MS); spikeless runs stop at NO_SPIKE_STOP_MS.
    Returns a list of MembraneTrace (one per batch member).
    """
    rest = resting_state(p)
    T_C, amp, gNa = np.broadcast_arrays(
        np.atleast_1d(np.asarray(T_C, float)),
        np.atleast_1d(np.asarray(amp, float)),
        np.atleast_1d(np.asarray(p.g_Na if g_Na is None else g_Na, float)))
    B = T_C.size
    phi_T = phi(T_C)
    adaptive = t_end is None
    t_stop = (start + duration + P.CHUNK_MS) if adaptive else t_end
    pieces = pulse_pieces(amp, t_stop, start, duration)

    n_total = sum(steps_in(a, b, dt) for a, b, _ in pieces)
    cap = steps_in(0.0, P.T_MAX_MS, dt)
    Y = np.empty((max(n_total, cap) + 1, 4, B))
    I_st = np.zeros((Y.shape[0], B))
    Y[0] = np.array([[rest.V], [rest.m], [rest.h], [rest.n]])
    k = 0
    y = Y[0].copy()

    def advance(n_steps, I):
        nonlocal k, y
        I = np.broadcast_to(I, (B,))
        for _ in range(n_steps):
            I_st[k] = I               # stimulus active during step k -> k+1
            y = rk4_step(y, dt, I, phi_T, gNa, p)
            k += 1
            Y[k] = y

    for (a, b, I) in pieces:
        advance(steps_in(a, b, dt), I)

    ends = [None] * B
    met = [False] * B
    if adaptive:
        chunk = steps_in(0.0, P.CHUNK_MS, dt)
        while True:
            t_now = np.arange(k + 1) * dt
            for i in range(B):
                if ends[i] is None:
                    ends[i], met[i] = _end_index(Y[:k + 1, 0, i], t_now, rest.V)
            if all(e is not None for e in ends) or k >= cap:
                break
            advance(min(chunk, cap - k), 0.0)
        ends = [k if e is None else e for e in ends]
    else:
        ends = [k] * B
        t_now = np.arange(k + 1) * dt
        met = [_end_index(Y[:k + 1, 0, i], t_now, rest.V)[1] for i in range(B)]

    traces = []
    for i in range(B):
        e = ends[i] + 1
        V, m, h, n = (Y[:e, j, i].copy() for j in range(4))
        I_Na, I_K, I_L = ionic_currents(V, m, h, n, p, gNa[i])
        traces.append(MembraneTrace(
            t=np.arange(e) * dt, V=V, m=m, h=h, n=n, I_Na=I_Na, I_K=I_K, I_L=I_L,
            I_stim=I_st[:e, i].copy(), T_C=float(T_C[i]), amp=float(amp[i]),
            g_Na=float(gNa[i]), dt=dt, rest=rest, stim_start=start,
            stim_end=start + duration, back_to_rest=bool(met[i])))
    return traces


def simulate_one(T_C, amp, **kw) -> MembraneTrace:
    return simulate(T_C, amp, **kw)[0]


# ---------------------------------------------------------------------------
# Reference solution (scipy DOP853), integrated piece by piece
# ---------------------------------------------------------------------------
def simulate_reference(T_C, amp, t_end, dt, p=P.HH, g_Na=None,
                       start=P.STIM_START_MS, duration=P.STIM_DURATION_MS) -> MembraneTrace:
    """Same model with solve_ivp(DOP853, rtol=1e-10, atol=1e-12), sampled on the RK4 time grid."""
    rest = resting_state(p)
    gNa = p.g_Na if g_Na is None else g_Na
    phi_T = float(phi(T_C))
    n_total = steps_in(0.0, t_end, dt)
    t_grid = np.arange(n_total + 1) * dt
    Y = np.empty((n_total + 1, 4))
    Y[0] = [rest.V, rest.m, rest.h, rest.n]
    I_st = np.zeros(n_total + 1)
    y0 = Y[0].copy()
    k0 = 0
    for (a, b, I) in pulse_pieces(amp, t_end, start, duration):
        n = steps_in(a, b, dt)
        t_eval = t_grid[k0:k0 + n + 1]
        sol = solve_ivp(lambda t, y: rhs(y, I, phi_T, gNa, p), (t_eval[0], t_eval[-1]), y0,
                        method=P.REF_METHOD, rtol=P.REF_RTOL, atol=P.REF_ATOL, t_eval=t_eval)
        assert sol.success, sol.message
        Y[k0:k0 + n + 1] = sol.y.T
        I_st[k0:k0 + n] = I
        y0 = sol.y[:, -1]
        k0 += n
    V, m, h, n = Y.T
    I_Na, I_K, I_L = ionic_currents(V, m, h, n, p, gNa)
    met = _end_index(V, t_grid, rest.V)[1]
    return MembraneTrace(t=t_grid, V=V, m=m, h=h, n=n, I_Na=I_Na, I_K=I_K, I_L=I_L,
                         I_stim=I_st, T_C=float(T_C), amp=float(amp), g_Na=float(gNa), dt=dt,
                         rest=rest, stim_start=start, stim_end=start + duration,
                         back_to_rest=bool(met))
