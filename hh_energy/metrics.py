"""
Measurements on one recorded spike (membrane patch or one node of the cable).

Inputs are 1-D arrays on a time grid t (ms): V (mV) and the ionic currents
I_Na, I_K (µA/cm², outward positive, so I_Na < 0 when sodium enters).
Only here is the sign flipped: Ĩ_Na = −I_Na (positive = sodium entering).
All charges are in nC/cm² (µA/cm² × ms).

Integrals use the trapezoid rule on the samples. A partial integral that
ends between two samples (e.g. at t_peak) uses linear interpolation of the
cumulative integral.
"""
from dataclasses import dataclass, asdict

import numpy as np
from scipy.integrate import cumulative_trapezoid

from . import params as P


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def cumulative(t, y):
    """Cumulative trapezoid integral of y over t, starting at 0."""
    return cumulative_trapezoid(y, t, initial=0.0)


def peak(t, V):
    """(t_peak, V_peak, i_peak): sample maximum refined by a parabola through 3 samples."""
    i = int(np.argmax(V))
    if 0 < i < len(V) - 1:
        y0, y1, y2 = V[i - 1], V[i], V[i + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom < 0:
            d = 0.5 * (y0 - y2) / denom          # offset in samples, |d| ≤ 0.5
            dt = t[i + 1] - t[i]
            return t[i] + d * dt, y1 - 0.25 * (y0 - y2) * d, i
    return float(t[i]), float(V[i]), i


def crossing(t, V, level, start=0, direction="up"):
    """First time at or after sample `start` where V crosses `level` (linear interpolation).

    direction "up": V goes from below to ≥ level; "down": from above to ≤ level.
    Returns nan if there is no crossing.
    """
    s = V[start:] - level
    if direction == "up":
        idx = np.nonzero((s[:-1] < 0) & (s[1:] >= 0))[0]
    else:
        idx = np.nonzero((s[:-1] > 0) & (s[1:] <= 0))[0]
    if idx.size == 0:
        return np.nan
    j = start + int(idx[0])
    f = (level - V[j]) / (V[j + 1] - V[j])
    return float(t[j] + f * (t[j + 1] - t[j]))


def last_crossing_before(t, V, level, stop, direction="up"):
    """Last crossing of `level` before sample `stop` (linear interpolation); nan if none."""
    s = V[:stop + 1] - level
    if direction == "up":
        idx = np.nonzero((s[:-1] < 0) & (s[1:] >= 0))[0]
    else:
        idx = np.nonzero((s[:-1] > 0) & (s[1:] <= 0))[0]
    if idx.size == 0:
        return np.nan
    j = int(idx[-1])
    f = (level - V[j]) / (V[j + 1] - V[j])
    return float(t[j] + f * (t[j + 1] - t[j]))


def atp_per_cm2(Q_Na_nC):
    """Number of ATP molecules per cm² to pump out Q_Na (3 Na⁺ per ATP)."""
    return Q_Na_nC * 1e-9 / (P.NA_PER_ATP * P.ELEMENTARY_CHARGE_C)


def membrane_spiked(V_peak, V_peak_passive):
    """Membrane spike criterion: V_peak > 0 mV and > passive-control peak + 20 mV."""
    return bool(V_peak > P.SPIKE_V_MIN_MV and V_peak - V_peak_passive > P.SPIKE_MARGIN_OVER_PASSIVE_MV)


def cable_propagated(V_far):
    """Travelling spike criterion: V at x = 4.0 cm exceeds 0 mV at some time."""
    return bool(np.max(V_far) > P.PROPAGATION_V_MV)


def velocity_m_s(t, V_a, V_b, x_a, x_b):
    """θ = (x_b − x_a)/(t_b − t_a) from upward 0 mV crossings; cm/ms × 10 = m/s."""
    t_a = crossing(t, V_a, 0.0, direction="up")
    t_b = crossing(t, V_b, 0.0, direction="up")
    return 10.0 * (x_b - x_a) / (t_b - t_a), t_a, t_b


# ---------------------------------------------------------------------------
# Main measurement
# ---------------------------------------------------------------------------
@dataclass
class ChargeMetrics:
    """Charges for one definition of Q_Na (nC/cm²)."""
    Q_Na: float
    r: float              # Q_Na / Q_min
    eta_pct: float        # efficiency 100 % / r
    Q_dep: float          # overlap A: Na⁺ charge from window start to t_peak
    Q_ov_A: float         # overlap A: Q_Na − Q_dep
    ov_frac_A: float      # Q_ov_A / Q_Na
    Q_ov_B: float         # overlap B: ∫ min(max(ΔĨ_Na,0), max(ΔI_K,0)) dt
    ov_frac_B: float      # Q_ov_B / Q_Na
    N_ATP_per_cm2: float
    N_ATP_per_um2: float
    window_start: float   # ms
    window_end: float     # ms


@dataclass
class SpikeMetrics:
    V_rest: float
    V_peak: float
    t_peak: float
    Q_min: float          # C_m (V_peak − V_rest)
    width: float          # ms, full width at half amplitude
    V_ahp: float          # after-hyperpolarisation minimum (mV)
    t_ahp: float
    t_end: float          # end of the recorded run (ms)
    primary: ChargeMetrics
    raw: ChargeMetrics

    def flat(self):
        """Dict with primary/raw fields prefixed, for tables."""
        d = {k: v for k, v in asdict(self).items() if k not in ("primary", "raw")}
        return d, asdict(self.primary), asdict(self.raw)


def _charges(t, Q_cum, ovB_integrand, t0, t1, t_peak, Q_min):
    """Charge metrics for the window [t0, t1] given the cumulative Na⁺ integral."""
    Q = lambda ta, tb: float(np.interp(tb, t, Q_cum) - np.interp(ta, t, Q_cum))
    Q_Na = Q(t0, t1)
    Q_dep = Q(t0, t_peak)
    ovB_cum = cumulative(t, ovB_integrand)
    Q_ov_B = float(np.interp(t1, t, ovB_cum) - np.interp(t0, t, ovB_cum))
    r = Q_Na / Q_min
    N = atp_per_cm2(Q_Na)
    return ChargeMetrics(Q_Na=Q_Na, r=r, eta_pct=100.0 / r, Q_dep=Q_dep, Q_ov_A=Q_Na - Q_dep,
                         ov_frac_A=(Q_Na - Q_dep) / Q_Na, Q_ov_B=Q_ov_B, ov_frac_B=Q_ov_B / Q_Na,
                         N_ATP_per_cm2=N, N_ATP_per_um2=N / P.UM2_PER_CM2,
                         window_start=float(t0), window_end=float(t1))


def spike_metrics(t, V, I_Na, I_K, rest, C_m=P.HH.C_m, raw_start=None) -> SpikeMetrics:
    """Every per-spike measurement, for both Q_Na definitions.

    rest: RestState (V, I_Na, I_K at the resting steady state).
    raw_start: start of the raw window in ms (membrane: stimulus onset). If None
    (cable), it is the first time V > V_rest + 1 mV.

    Primary Q_Na = ∫ (Ĩ_Na − Ĩ_Na,rest) dt over the whole run.
    Raw Q_Na     = ∫ Ĩ_Na dt from raw_start to the first time after the AHP
                   minimum that V is within 1 mV of V_rest.
    """
    t = np.asarray(t, float)
    V = np.asarray(V, float)
    Na_in = -np.asarray(I_Na, float)            # Ĩ_Na, positive = entering
    Na_in_rest = -rest.I_Na
    t_peak, V_peak, i_peak = peak(t, V)
    Q_min = C_m * (V_peak - rest.V)

    # Width at half amplitude
    V_half = rest.V + 0.5 * (V_peak - rest.V)
    width = (crossing(t, V, V_half, start=i_peak, direction="down")
             - last_crossing_before(t, V, V_half, i_peak, direction="up"))

    # After-hyperpolarisation minimum
    i_ahp = i_peak + int(np.argmin(V[i_peak:]))
    V_ahp, t_ahp = float(V[i_ahp]), float(t[i_ahp])

    # Overlap-B integrand (rest-subtracted, positive parts only)
    dNa = Na_in - Na_in_rest
    dK = np.asarray(I_K, float) - rest.I_K
    ovB = np.minimum(np.maximum(dNa, 0.0), np.maximum(dK, 0.0))

    # Primary definition: whole run, rest-subtracted
    primary = _charges(t, cumulative(t, dNa), ovB, t[0], t[-1], t_peak, Q_min)

    # Raw definition: adaptive window, not rest-subtracted
    if raw_start is None:
        raw_start = crossing(t, V, rest.V + P.RAW_ONSET_MARGIN_MV, direction="up")
    level = rest.V - P.RAW_END_TOL_MV
    raw_end = t_ahp if V_ahp >= level else crossing(t, V, level, start=i_ahp, direction="up")
    raw = _charges(t, cumulative(t, Na_in), ovB, raw_start, raw_end, t_peak, Q_min)

    return SpikeMetrics(V_rest=rest.V, V_peak=float(V_peak), t_peak=float(t_peak), Q_min=float(Q_min),
                        width=float(width), V_ahp=V_ahp, t_ahp=t_ahp, t_end=float(t[-1]),
                        primary=primary, raw=raw)
