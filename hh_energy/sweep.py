"""
Experiments: stimulus threshold, temperature sweeps, failure-temperature
bisection, grid refinement, validation, membrane-vs-travelling comparison.

Units: T °C, amplitudes µA/cm², charges nC/cm², velocity m/s, times ms.
"""
import json
import math
from dataclasses import asdict

import numpy as np
import pandas as pd

from . import params as P
from . import cable, membrane, metrics
from .rates import phi


# ---------------------------------------------------------------------------
# Membrane spike helpers
# ---------------------------------------------------------------------------
def passive_peaks(T_C, amp, dt=P.DT_MEMBRANE_MS):
    """Peak V (mV) of the passive control (ḡ_Na = 0, same start state and stimulus)."""
    runs = membrane.simulate(T_C, amp, g_Na=0.0, dt=dt, t_end=P.PASSIVE_RUN_MS)
    return np.array([metrics.peak(r.t, r.V)[1] for r in runs])


def membrane_runs(T_C, amp, dt=P.DT_MEMBRANE_MS):
    """Spike runs plus their passive controls (batched). Returns (traces, spiked flags, passive peaks)."""
    traces = membrane.simulate(T_C, amp, dt=dt)
    T_b, a_b = np.broadcast_arrays(np.atleast_1d(T_C), np.atleast_1d(amp))
    pp = passive_peaks(T_b, a_b, dt)
    spiked = [metrics.membrane_spiked(metrics.peak(tr.t, tr.V)[1], p_) for tr, p_ in zip(traces, pp)]
    return traces, spiked, pp


def membrane_spikes(T_C, amp, dt=P.DT_MEMBRANE_MS):
    return membrane_runs(T_C, amp, dt)[1][0]


def find_threshold(T_C=P.THRESHOLD_T_C, tol=P.THRESHOLD_TOL_UA, dt=P.DT_MEMBRANE_MS):
    """Smallest pulse amplitude (µA/cm²) that gives a membrane spike, by bisection.

    Returns (threshold, lo, hi): lo does not spike, hi = threshold spikes, hi − lo ≤ tol.
    """
    lo, hi = P.THRESHOLD_BRACKET_UA
    assert not membrane_spikes(T_C, lo, dt) and membrane_spikes(T_C, hi, dt)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if membrane_spikes(T_C, mid, dt):
            hi = mid
        else:
            lo = mid
    return hi, lo, hi


def stimulus_amplitudes(threshold):
    return {"amp2x": P.AMP_FACTOR_PRIMARY * threshold, "amp4x": P.AMP_FACTOR_VARIANT * threshold}


# ---------------------------------------------------------------------------
# One temperature: all spike types and definitional variants
# ---------------------------------------------------------------------------
def _rows(base, sm):
    """Two rows (primary and raw Q_Na definitions) for one spike."""
    common, prim, raw = sm.flat()
    return [{**base, **common, "qna_def": name, **charges}
            for name, charges in (("primary", prim), ("raw", raw))]


def _failed_rows(base):
    return [{**base, "qna_def": name} for name in ("primary", "raw")]


def membrane_at(T_C, amps, dt=P.DT_MEMBRANE_MS):
    """Membrane spikes at both amplitudes. Returns (rows, {variant: spiked}, {variant: trace})."""
    names = list(amps)
    traces, spiked, pp = membrane_runs(T_C, np.array([amps[k] for k in names]), dt)
    rows, ok, keep = [], {}, {}
    for name, tr, s, p_ in zip(names, traces, spiked, pp):
        base = {"T_C": T_C, "phi": float(phi(T_C)), "spike_type": "membrane", "stim_variant": name,
                "stim_amp_uA_cm2": amps[name], "spiked": s, "dx_um": np.nan, "dt_us": dt * 1000,
                "passive_peak_mV": p_, "V_peak_sample_mV": float(np.max(tr.V)),
                "back_to_rest": tr.back_to_rest}
        if s:
            assert tr.back_to_rest, f"membrane run at {T_C} °C did not return to rest"
            sm = metrics.spike_metrics(tr.t, tr.V, tr.I_Na, tr.I_K, tr.rest, raw_start=tr.stim_start)
            rows += _rows(base, sm)
        else:
            rows += _failed_rows(base)
        ok[name], keep[name] = s, tr
    return rows, ok, keep


def cable_at(T_C, dx=P.DX_CABLE_CM, dt=P.DT_CABLE_MS):
    """Travelling spike at T_C. Returns (rows, propagated, run)."""
    run = cable.simulate(T_C, dx=dx, dt=dt)
    xm = P.MEASURE_X_CM
    ok = metrics.cable_propagated(run.V_rec[P.PROPAGATION_X_CM])
    base = {"T_C": T_C, "phi": float(phi(T_C)), "spike_type": "travelling", "stim_variant": "cable",
            "stim_amp_uA_cm2": P.CABLE_STIM_AMP_UA, "spiked": ok, "dx_um": dx * 1e4, "dt_us": dt * 1000,
            "V_max_4cm_mV": float(np.max(run.V_rec[P.PROPAGATION_X_CM])),
            "V_peak_sample_mV": float(np.max(run.V_rec[xm])), "back_to_rest": run.back_to_rest}
    if not ok:
        return _failed_rows(base), ok, run
    assert run.back_to_rest, f"cable run at {T_C} °C: x = 2.5 cm did not return to rest ({run.stop_reason})"
    xa, xb = P.VELOCITY_X_CM
    base["velocity_m_s"] = metrics.velocity_m_s(run.t, run.V_rec[xa], run.V_rec[xb], xa, xb)[0]
    sm = metrics.spike_metrics(run.t, run.V_rec[xm], run.I_Na_rec[xm], run.I_K_rec[xm], run.rest)
    return _rows(base, sm), ok, run


# ---------------------------------------------------------------------------
# 1. Coarse sweep
# ---------------------------------------------------------------------------
def coarse_temperatures():
    """6.3, then 8, 10, 12, ... (2 °C steps) with 18.5 inserted in order."""
    extra = sorted(P.SWEEP_T_EXTRA_C)
    yield P.SWEEP_T_FIRST_C
    T = P.SWEEP_T_GRID_START_C
    while T <= P.SWEEP_T_CAP_C:
        while extra and extra[0] < T:
            yield extra.pop(0)
        if extra and abs(extra[0] - T) < 1e-9:
            extra.pop(0)
        yield T
        T += P.SWEEP_T_STEP_C


def coarse_sweep(amps, log=print):
    """Sweep upward until the membrane spike (both amplitudes) and the travelling spike have all failed.

    Returns (rows, status DataFrame, membrane traces {T: {variant: trace}}, cable runs {T: run}).
    """
    rows, status, mem_traces, cab_runs = [], [], {}, {}
    for T in coarse_temperatures():
        r_m, ok_m, tr_m = membrane_at(T, amps)
        r_c, ok_c, run_c = cable_at(T)
        rows += r_m + r_c
        mem_traces[T], cab_runs[T] = tr_m, run_c
        status.append({"T_C": T, **{f"membrane_{k}": v for k, v in ok_m.items()}, "travelling": ok_c})
        log(f"  T = {T:5.1f} °C  membrane 2x: {ok_m['amp2x']!s:5}  4x: {ok_m['amp4x']!s:5}  "
            f"travelling: {ok_c!s:5}")
        if not any(ok_m.values()) and not ok_c:
            break
    else:
        raise RuntimeError("spikes did not fail below SWEEP_T_CAP_C")
    return rows, pd.DataFrame(status), mem_traces, cab_runs


# ---------------------------------------------------------------------------
# 2. Failure temperatures
# ---------------------------------------------------------------------------
def _bracket(status, column):
    """(last success, first failure) on the coarse grid; checks failure is monotonic."""
    ok = status[column].to_numpy(bool)
    T = status["T_C"].to_numpy(float)
    first_fail = int(np.argmin(ok))
    assert not ok[first_fail], f"{column}: never failed"
    monotonic = bool(ok[:first_fail].all() and not ok[first_fail:].any())
    return T[first_fail - 1], T[first_fail], monotonic


def bisect_failure(succeeds, lo, hi, tol=P.FAILURE_BISECT_TOL_C):
    """Bisection on T between a success (lo) and a failure (hi) until hi − lo ≤ tol."""
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if succeeds(mid):
            lo = mid
        else:
            hi = mid
    return lo, hi


def failure_temperatures(status, amps, log=print):
    tests = {
        "membrane_amp2x": lambda T: membrane_spikes(T, amps["amp2x"]),
        "membrane_amp4x": lambda T: membrane_spikes(T, amps["amp4x"]),
        "travelling": lambda T: cable_at(T)[1],
    }
    out = []
    for name, succeeds in tests.items():
        lo0, hi0, mono = _bracket(status, name)
        lo, hi = bisect_failure(succeeds, lo0, hi0)
        out.append({"spike": name, "coarse_last_success_C": lo0, "coarse_first_failure_C": hi0,
                    "last_success_C": lo, "first_failure_C": hi, "T_fail_C": 0.5 * (lo + hi),
                    "half_width_C": 0.5 * (hi - lo), "monotonic_on_coarse_grid": mono})
        log(f"  {name}: fails between {lo:.4f} and {hi:.4f} °C")
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# 2b. Diagnostic failure temperatures (extra rules, reported separately and labelled as such)
# ---------------------------------------------------------------------------
def membrane_active(T_C, amp, dt=P.DT_MEMBRANE_MS):
    """Diagnostic: active response > passive peak + 20 mV, without the 0 mV requirement."""
    traces, _, pp = membrane_runs(T_C, amp, dt)
    return bool(np.max(traces[0].V) - pp[0] > P.SPIKE_MARGIN_OVER_PASSIVE_MV)


def wave_nondecremental(run):
    """Diagnostic: the wave keeps its size from x = 2.0 to 4.0 cm (does not die out)."""
    xa, xb = P.DIAG_WAVE_X_CM
    d_a = np.max(run.V_rec[xa]) - run.rest.V
    d_b = np.max(run.V_rec[xb]) - run.rest.V
    return bool(d_a >= P.DIAG_WAVE_MIN_DEPOL_MV and d_b >= P.DIAG_WAVE_RATIO * d_a)


SPEC_CRITERIA = {
    "membrane_amp2x": "V_peak > 0 mV and > passive peak + 20 mV (main rule)",
    "membrane_amp4x": "V_peak > 0 mV and > passive peak + 20 mV (main rule)",
    "travelling": "V(x = 4.0 cm) > 0 mV (main rule)",
}
DIAG_CRITERIA = {
    "membrane_amp2x": "DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement)",
    "membrane_amp4x": "DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement)",
    "travelling": "DIAGNOSTIC: depolarisation at 4.0 cm >= 90% of that at 2.0 cm (>= 20 mV)",
}


def diagnostic_failure_temperatures(amps, log=print):
    """Coarse scan (same T grid) then bisection, for the diagnostic rules."""
    tests = {
        "membrane_amp2x": lambda T: membrane_active(T, amps["amp2x"]),
        "membrane_amp4x": lambda T: membrane_active(T, amps["amp4x"]),
        "travelling": lambda T: wave_nondecremental(cable.simulate(T)),
    }
    out = []
    for name, succeeds in tests.items():
        last_ok = None
        for T in coarse_temperatures():
            if succeeds(T):
                last_ok = T
            else:
                break
        assert last_ok is not None, f"{name}: diagnostic fails already at {T} °C"
        lo, hi = bisect_failure(succeeds, last_ok, T)
        out.append({"spike": name, "criterion": DIAG_CRITERIA[name], "is_main_rule": False,
                    "coarse_last_success_C": last_ok, "coarse_first_failure_C": T,
                    "last_success_C": lo, "first_failure_C": hi, "T_fail_C": 0.5 * (lo + hi),
                    "half_width_C": 0.5 * (hi - lo), "monotonic_on_coarse_grid": True})
        log(f"  {name} (diagnostic): response disappears between {lo:.4f} and {hi:.4f} °C")
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# 3. Validation: HH 1952 velocity and RK4 vs DOP853
# ---------------------------------------------------------------------------
def validation(amps, log=print):
    """Rows for validation.csv, plus the 6.3 °C RK4/DOP853 traces for the figure."""
    rows, keep = [], {}
    for i, (dx, dt) in enumerate(P.CABLE_REFINEMENT):
        _, ok, run = cable_at(P.VALIDATION_T_C, dx, dt)
        xa, xb = P.VELOCITY_X_CM
        theta = metrics.velocity_m_s(run.t, run.V_rec[xa], run.V_rec[xb], xa, xb)[0] if ok else np.nan
        rows.append({
            "check": "conduction velocity vs Hodgkin & Huxley 1952", "T_C": P.VALIDATION_T_C,
            "grid": f"dx = {dx * 1e4:g} um, dt = {dt * 1000:g} us", "is_base_grid": i == 0,
            "model_value": theta, "unit": "m/s",
            "reference_value": P.HH1952_VELOCITY_COMPUTED_M_S, "reference": "HH 1952 computed",
            "pct_difference": 100.0 * (theta - P.HH1952_VELOCITY_COMPUTED_M_S) / P.HH1952_VELOCITY_COMPUTED_M_S,
            "measured_value": P.HH1952_VELOCITY_MEASURED_M_S,
            "pct_difference_vs_measured": 100.0 * (theta - P.HH1952_VELOCITY_MEASURED_M_S) / P.HH1952_VELOCITY_MEASURED_M_S,
        })
        log(f"  velocity at {P.VALIDATION_T_C} °C, dx = {dx * 1e4:g} um: {theta:.3f} m/s")
    for T in (6.3, 25.0):
        tr = membrane.simulate_one(T, amps["amp2x"])
        ref = membrane.simulate_reference(T, amps["amp2x"], tr.t[-1], tr.dt)
        sm = metrics.spike_metrics(tr.t, tr.V, tr.I_Na, tr.I_K, tr.rest, raw_start=tr.stim_start)
        sr = metrics.spike_metrics(ref.t, ref.V, ref.I_Na, ref.I_K, ref.rest, raw_start=ref.stim_start)
        dV = np.abs(tr.V - ref.V)
        in_spike = (tr.t >= sm.raw.window_start) & (tr.t <= sm.raw.window_end)
        common = {"T_C": T, "grid": f"RK4 dt = {tr.dt * 1000:g} us vs DOP853 rtol=1e-10 atol=1e-12",
                  "is_base_grid": True, "reference": "DOP853"}
        rows += [
            {**common, "check": "RK4 vs DOP853: max |dV| during the spike (raw window)",
             "model_value": float(dV[in_spike].max()), "unit": "mV"},
            {**common, "check": "RK4 vs DOP853: max |dV| over the whole run",
             "model_value": float(dV.max()), "unit": "mV"},
            {**common, "check": "RK4 vs DOP853: relative difference in Q_Na (primary)",
             "model_value": (sm.primary.Q_Na - sr.primary.Q_Na) / sr.primary.Q_Na, "unit": "fraction",
             "reference_value": sr.primary.Q_Na},
            {**common, "check": "RK4 vs DOP853: relative difference in r (primary)",
             "model_value": (sm.primary.r - sr.primary.r) / sr.primary.r, "unit": "fraction",
             "reference_value": sr.primary.r},
        ]
        log(f"  RK4 vs DOP853 at {T} °C: max |dV| = {dV.max():.2e} mV, "
            f"rel dQ_Na = {(sm.primary.Q_Na - sr.primary.Q_Na) / sr.primary.Q_Na:.2e}")
        keep[T] = (tr, ref)
    return pd.DataFrame(rows), keep


# ---------------------------------------------------------------------------
# 4. Grid refinement
# ---------------------------------------------------------------------------
REFINED_QUANTITIES = {"travelling": ("r", "Q_Na", "Q_min", "velocity_m_s"),
                      "membrane": ("r", "Q_Na", "Q_min")}


def refinement(amps, temperatures, log=print):
    """Long table: model, T, level, dx, dt, quantity, value, |value − value(base)|."""
    rows = []
    for T in temperatures:
        for i, (dx, dt) in enumerate(P.CABLE_REFINEMENT):
            r_c, ok, _ = cable_at(T, dx, dt)
            assert ok, f"travelling spike failed at {T} °C, dx = {dx}"
            prim = r_c[0]
            for q in REFINED_QUANTITIES["travelling"]:
                rows.append({"model": "travelling", "T_C": T, "level": i, "dx_um": dx * 1e4,
                             "dt_us": dt * 1000, "quantity": q, "value": prim[q]})
            log(f"  travelling {T} °C dx = {dx * 1e4:g} um dt = {dt * 1000:g} us: "
                f"r = {prim['r']:.5f}, theta = {prim['velocity_m_s']:.4f} m/s")
        for i, dt in enumerate(P.MEMBRANE_REFINEMENT):
            r_m, ok, _ = membrane_at(T, {"amp2x": amps["amp2x"]}, dt)
            assert ok["amp2x"], f"membrane spike failed at {T} °C, dt = {dt}"
            prim = r_m[0]
            for q in REFINED_QUANTITIES["membrane"]:
                rows.append({"model": "membrane", "T_C": T, "level": i, "dx_um": np.nan,
                             "dt_us": dt * 1000, "quantity": q, "value": prim[q]})
            log(f"  membrane   {T} °C dt = {dt * 1000:g} us: r = {prim['r']:.6f}")
    df = pd.DataFrame(rows)
    base = df[df.level == 0].set_index(["model", "T_C", "quantity"])["value"]
    df["abs_diff_from_base"] = [abs(v - base[(m, T, q)]) for m, T, q, v in
                                zip(df.model, df.T_C, df.quantity, df.value)]
    return df


def epsilons(ref_df):
    """ε per model, quantity and T, and the value used everywhere (max over the two T)."""
    eps = ref_df.groupby(["model", "quantity", "T_C"])["abs_diff_from_base"].max().unstack("T_C")
    eps.columns = [f"eps_at_{c:g}C" for c in eps.columns]
    eps["eps_used"] = eps.max(axis=1)
    return eps.reset_index()


# ---------------------------------------------------------------------------
# 4b. Extra robustness checks near the failure temperature
# ---------------------------------------------------------------------------
def edge_checks(amps, status, fail_spec, log=print):
    """r of both impulse types close to the failure temperature, and grid refinement there.

    T_edge: warmest coarse temperature at which the membrane spike (2×) and the travelling spike
    both exist. T_last_both: warmest temperature at which the failure bisections found both.
    Returns (near-failure table, refinement table at T_edge in the format of refinement()).
    """
    both = status[status.membrane_amp2x.astype(bool) & status.travelling.astype(bool)]
    T_edge = float(both.T_C.max())
    last = fail_spec.set_index("spike").last_success_C
    T_last = float(min(last["membrane_amp2x"], last["travelling"]))
    rows = []
    for T, which in ((T_edge, "warmest sweep temperature with both impulses"),
                     (T_last, "warmest temperature with both impulses (from the failure bisections)")):
        r_m, ok_m, _ = membrane_at(T, {"amp2x": amps["amp2x"]})
        r_c, ok_c, _ = cable_at(T)
        assert ok_m["amp2x"] and ok_c, f"both impulses should exist at {T} °C"
        m, c = r_m[0], r_c[0]
        pct = 100.0 * (c["r"] - m["r"]) / m["r"]
        rows.append({"T_C": T, "which": which, "r_membrane": m["r"], "r_travel": c["r"], "pct_difference": pct,
                     "Q_Na_membrane": m["Q_Na"], "Q_Na_travel": c["Q_Na"],
                     "Q_min_membrane": m["Q_min"], "Q_min_travel": c["Q_min"],
                     "V_peak_membrane_mV": m["V_peak"], "V_peak_travel_mV": c["V_peak"]})
        log(f"  {T:.5g} °C ({which}): r_membrane = {m['r']:.4f}, r_travel = {c['r']:.4f}, difference = {pct:.2f} %")
    log(f"  grid refinement at {T_edge:g} °C:")
    ref_edge = refinement(amps, (T_edge,), log)
    return pd.DataFrame(rows), ref_edge, T_edge, T_last


def failure_vs_grid(amps, fail_spec, log=print):
    """Re-find the main-rule failure temperatures on every refinement grid (bisection to FAILURE_GRID_TOL_C).

    Uses the coarse bracket from the base-grid sweep; asserts that it still brackets the failure.
    """
    rows = []
    br = fail_spec.set_index("spike")
    cases = [("travelling", lvl, dx, dt, (lambda T, dx=dx, dt=dt: cable_at(T, dx, dt)[1]))
             for lvl, (dx, dt) in enumerate(P.CABLE_REFINEMENT)]
    cases += [("membrane_amp2x", lvl, np.nan, dt, (lambda T, dt=dt: membrane_spikes(T, amps["amp2x"], dt)))
              for lvl, dt in enumerate(P.MEMBRANE_REFINEMENT)]
    for spike, lvl, dx, dt, ok in cases:
        lo, hi = float(br.loc[spike].coarse_last_success_C), float(br.loc[spike].coarse_first_failure_C)
        assert ok(lo) and not ok(hi), f"{spike}: [{lo}, {hi}] °C no longer brackets the failure on grid {lvl}"
        lo, hi = bisect_failure(ok, lo, hi, P.FAILURE_GRID_TOL_C)
        rows.append({"spike": spike, "level": lvl, "dx_um": dx * 1e4, "dt_us": dt * 1000,
                     "last_success_C": lo, "first_failure_C": hi, "T_fail_C": 0.5 * (lo + hi),
                     "half_width_C": 0.5 * (hi - lo)})
        log(f"  {spike} grid {lvl} (dx = {dx * 1e4:g} um, dt = {dt * 1000:g} us): fails at {0.5 * (lo + hi):.3f} °C")
    df = pd.DataFrame(rows)
    base = df[df.level == 0].set_index("spike").T_fail_C
    df["shift_from_base_C"] = [T - base[s] for s, T in zip(df.spike, df.T_fail_C)]
    return df


# ---------------------------------------------------------------------------
# 5. Membrane vs travelling comparison
# ---------------------------------------------------------------------------
def comparison(sweep_df, eps_df):
    """Δr(T) = r_travel − r_membrane with the flag |Δr| > 2ε, ε = sqrt(ε_mem² + ε_cab²)."""
    e = eps_df.set_index(["model", "quantity"])["eps_used"]
    eps = math.sqrt(e[("membrane", "r")] ** 2 + e[("travelling", "r")] ** 2)
    ok = sweep_df[sweep_df.spiked.astype(bool)]
    pick = lambda typ, var, d: ok[(ok.spike_type == typ) & (ok.stim_variant == var) & (ok.qna_def == d)] \
        .set_index("T_C")["r"]
    variants = {"": ("amp2x", "primary"), "_raw": ("amp2x", "raw"), "_amp4x": ("amp4x", "primary")}
    out = None
    for suffix, (var, d) in variants.items():
        df = pd.DataFrame({f"r_membrane{suffix}": pick("membrane", var, d),
                           f"r_travel{suffix}": pick("travelling", "cable", d)}).dropna()
        df[f"delta_r{suffix}"] = df[f"r_travel{suffix}"] - df[f"r_membrane{suffix}"]
        df[f"real_difference{suffix}"] = df[f"delta_r{suffix}"].abs() > P.DELTA_R_FLAG_FACTOR * eps
        out = df if out is None else out.join(df, how="outer")
    out.insert(0, "phi", phi(out.index.to_numpy()))
    out["eps_combined"] = eps
    out["two_eps"] = P.DELTA_R_FLAG_FACTOR * eps
    out["flag"] = np.where(out["real_difference"], "real difference", "within numerical uncertainty")
    return out.reset_index()


# ---------------------------------------------------------------------------
# Traces for figures / animation / web export
# ---------------------------------------------------------------------------
def save_traces(path, mem_traces, cab_runs, val_keep):
    """Store the simulated traces the figures need (npz; keys documented in load_traces)."""
    d = {}
    for T, by_var in mem_traces.items():
        for var, tr in by_var.items():
            for f in ("t", "V", "m", "h", "n", "I_Na", "I_K", "I_L"):
                d[f"mem|{T:g}|{var}|{f}"] = getattr(tr, f)
    for T, run in cab_runs.items():
        d[f"cab|{T:g}|t"] = run.t
        for xr in run.record_x:
            d[f"cab|{T:g}|V|{xr:g}"] = run.V_rec[xr]
            d[f"cab|{T:g}|I_Na|{xr:g}"] = run.I_Na_rec[xr]
            d[f"cab|{T:g}|I_K|{xr:g}"] = run.I_K_rec[xr]
        d[f"cab|{T:g}|x"] = run.x
        d[f"cab|{T:g}|t_snap"] = run.t_snap
        d[f"cab|{T:g}|V_snap"] = run.V_snap.astype(np.float32)
    for T, (tr, ref) in val_keep.items():
        d[f"val|{T:g}|t"] = tr.t
        d[f"val|{T:g}|V_rk4"] = tr.V
        d[f"val|{T:g}|V_ref"] = ref.V
    np.savez_compressed(path, **d)


def load_traces(path):
    """Nested dict: mem[T][variant][field], cab[T][field] (V/I_Na/I_K keyed by x), val[T][field]."""
    out = {"mem": {}, "cab": {}, "val": {}}
    with np.load(path) as z:
        for key in z.files:
            parts = key.split("|")
            kind, T = parts[0], float(parts[1])
            node = out[kind].setdefault(T, {})
            for p_ in parts[2:-1]:
                node = node.setdefault(float(p_) if p_[0].isdigit() else p_, {})
            last = parts[-1]
            node[float(last) if last[0].isdigit() else last] = z[key]
    return out


# ---------------------------------------------------------------------------
# Everything
# ---------------------------------------------------------------------------
def run_experiments(results_dir, log=print):
    """Run every experiment of section 7 and write the result files. Returns a summary dict."""
    from pathlib import Path
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    log(f"Explicit-Euler stability number D·dt/(C_m·dx²) = {cable.explicit_stability_number():.2f} "
        "(explicit stepping needs ≤ 0.5, so backward Euler is used)")
    log("Threshold at 6.3 °C (bisection):")
    threshold, th_lo, th_hi = find_threshold()
    amps = stimulus_amplitudes(threshold)
    log(f"  threshold = {threshold:.4f} uA/cm2 (bracket {th_lo:.4f}–{th_hi:.4f}); "
        f"primary 2x = {amps['amp2x']:.4f}, variant 4x = {amps['amp4x']:.4f}")

    log("Coarse temperature sweep:")
    rows, status, mem_traces, cab_runs = coarse_sweep(amps, log)
    sweep_df = pd.DataFrame(rows)

    log("Failure temperatures (main 0 mV rules):")
    fail_df = failure_temperatures(status, amps, log)
    fail_df.insert(1, "criterion", [SPEC_CRITERIA[s] for s in fail_df.spike])
    fail_df.insert(2, "is_main_rule", True)
    log("Failure temperatures (diagnostic rules):")
    diag_df = diagnostic_failure_temperatures(amps, log)
    fail_all = pd.concat([fail_df, diag_df], ignore_index=True)

    # Diagnostic flags for every coarse row (labelled, not used by the main analysis)
    diag_flag = []
    for _, row in sweep_df.iterrows():
        if row.spike_type == "membrane":
            diag_flag.append(bool(row.V_peak_sample_mV - row.passive_peak_mV > P.SPIKE_MARGIN_OVER_PASSIVE_MV))
        else:
            diag_flag.append(wave_nondecremental(cab_runs[row.T_C]))
    sweep_df["diagnostic_response_present"] = diag_flag

    T_prop = status.loc[status.travelling, "T_C"].max()
    T_warm = float(T_prop - P.T_WARM_OFFSET_C)
    log(f"Highest coarse T with a travelling spike: {T_prop:g} °C -> T_warm = {T_warm:g} °C")

    log("Validation:")
    val_df, val_keep = validation(amps, log)

    log("Grid refinement:")
    ref_df = refinement(amps, (P.SWEEP_T_FIRST_C, T_warm), log)
    eps_df = epsilons(ref_df)
    cmp_df = comparison(sweep_df, eps_df)

    log("Extra robustness checks near the failure temperature:")
    near_df, ref_edge_df, T_edge, T_last = edge_checks(amps, status, fail_df, log)
    log("Failure temperature on each refinement grid (extra robustness check):")
    fail_grid_df = failure_vs_grid(amps, fail_df, log)

    # Write files
    sweep_df.to_csv(results_dir / "sweep.csv", index=False)
    ref_df.to_csv(results_dir / "refinement.csv", index=False)
    eps_df.to_csv(results_dir / "epsilon.csv", index=False)
    fail_all.to_csv(results_dir / "failure_temps.csv", index=False)
    val_df.to_csv(results_dir / "validation.csv", index=False)
    cmp_df.to_csv(results_dir / "comparison.csv", index=False)
    status.to_csv(results_dir / "coarse_status.csv", index=False)
    near_df.to_csv(results_dir / "near_failure.csv", index=False)
    fail_grid_df.to_csv(results_dir / "failure_vs_grid.csv", index=False)
    ref_edge_df.to_csv(results_dir / "refinement_edge.csv", index=False)
    save_traces(results_dir / "traces.npz", mem_traces, cab_runs, val_keep)

    summary = {
        "threshold_uA_cm2": threshold, "threshold_bracket_uA_cm2": [th_lo, th_hi],
        "amp_primary_uA_cm2": amps["amp2x"], "amp_variant_uA_cm2": amps["amp4x"],
        "temperatures_C": [float(t) for t in status.T_C], "T_warm_C": T_warm,
        "T_highest_propagating_coarse_C": float(T_prop),
        "T_edge_C": T_edge, "T_last_both_C": T_last,
        "explicit_stability_number": cable.explicit_stability_number(),
        "rest": asdict(membrane.resting_state()),
        "n_nodes_base": int(cable.grid().size),
    }
    config = {**P.config_dict(), "derived": summary}
    (results_dir / "config.json").write_text(json.dumps(config, indent=2, default=float))
    return summary
