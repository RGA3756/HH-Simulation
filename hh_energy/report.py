"""
Writes results/env_info.txt, RESULTS_SUMMARY.md and METHODS_FILL_IN.md.

Every number in the generated documents is read from the files in results/
(written by the simulations and tests); nothing is typed in by hand.
Units: T °C, V mV, Q nC/cm², velocity m/s.
"""
import datetime
import json
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from . import params as P

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

TESTS = {
    "1": ("Temperature factor", "φ(6.3) = 1 and φ(16.3) = 3 to 1e−12"),
    "2": ("Removable singularities", "α_m(−40) = 1, α_n(−55) = 0.1; rates finite and continuous within 1e−6 mV"),
    "3": ("Resting state", "|V_rest − (−65)| < 0.5 mV and total ionic current < 1e−9 µA/cm²"),
    "4": ("Membrane charge balance", "∫I_stim dt = ∫I_ion dt + C_m ΔV within 0.1 % (6.3 and 25 °C)"),
    "5": ("RK4 vs DOP853", "max |ΔV| < 0.01 mV at 6.3 and 25 °C"),
    "6": ("Passive cable / units", "steady state matches cosh((L − x)/λ) and the input resistance within 1 %"),
    "7": ("Propagation", "travelling spike reaches x = 4 cm above 0 mV at 6.3 and 18.5 °C"),
}
VELOCITY_TOL_PCT = 5.0   # my criterion for the PASS/FAIL line of the HH 1952 benchmark (not in the spec)


def f(x, d=3):
    """Number formatting with a proper minus sign."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:.{d}f}".replace("-", "−")


def cell(text):
    """Escape '|' so text can sit inside a Markdown table cell."""
    return str(text).replace("|", "\\|")


def g(x, d=2):
    return f"{x:.{d}e}".replace("-", "−")


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
def _sysctl(key):
    try:
        return subprocess.run(["sysctl", "-n", key], capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return ""


def write_env_info(path=RESULTS / "env_info.txt"):
    import matplotlib
    import scipy
    info = {
        "OS": platform.platform(),
        "Machine (platform.machine)": platform.machine(),
        "Processor (platform.processor)": platform.processor() or "n/a",
        "Machine model (sysctl hw.model)": _sysctl("hw.model") or "n/a",
        "CPU (sysctl machdep.cpu.brand_string)": _sysctl("machdep.cpu.brand_string") or "n/a",
        "Python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "pandas": pd.__version__,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{k}: {v}\n" for k, v in info.items()))
    return info


def read_env_info(path=RESULTS / "env_info.txt"):
    return dict(line.split(": ", 1) for line in path.read_text().splitlines() if ": " in line)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def read_tests(path=RESULTS / "pytest_junit.xml"):
    cases = []
    for tc in ET.parse(path).getroot().iter("testcase"):
        name = tc.get("name")
        status = "FAIL" if (tc.find("failure") is not None or tc.find("error") is not None) else \
                 "SKIP" if tc.find("skipped") is not None else "PASS"
        measured = [p.get("value") for p in tc.iter("property") if p.get("name") == "measured"]
        cases.append({"name": name, "number": name.split("_")[1], "status": status, "measured": measured})
    return cases


def load():
    d = {n: pd.read_csv(RESULTS / f"{n}.csv") for n in
         ("sweep", "comparison", "epsilon", "failure_temps", "validation", "refinement", "js_consistency",
          "coarse_status", "near_failure", "refinement_edge", "failure_vs_grid")}
    d["config"] = json.loads((RESULTS / "config.json").read_text())
    d["env"] = read_env_info()
    d["tests"] = read_tests()
    return d


def pick(sweep, spike_type, variant, qna_def="primary"):
    s = sweep[(sweep.spike_type == spike_type) & (sweep.stim_variant == variant)
              & (sweep.qna_def == qna_def) & sweep.spiked.astype(bool)]
    return s.sort_values("T_C").set_index("T_C")


def eps(d, model, quantity, col="eps_used"):
    e = d["epsilon"]
    return float(e[(e.model == model) & (e.quantity == quantity)][col].iloc[0])


def fail_row(d, spike, spec=True):
    ft = d["failure_temps"]
    return ft[(ft.spike == spike) & (ft.is_spec_definition.astype(bool) == spec)].iloc[0]


# ---------------------------------------------------------------------------
# RESULTS_SUMMARY.md
# ---------------------------------------------------------------------------
def results_summary(d):
    sw, cmp_, val, cfg = d["sweep"], d["comparison"].set_index("T_C"), d["validation"], d["config"]["derived"]
    mem, mem4, trav = pick(sw, "membrane", "amp2x"), pick(sw, "membrane", "amp4x"), pick(sw, "travelling", "cable")
    T_warm = cfg["T_warm_C"]
    L = []
    L.append("# Results summary\n")
    L.append(f"*Generated automatically by `run_all.py` on {datetime.datetime.now():%Y-%m-%d %H:%M}. "
             "Do not edit by hand: every number below is read from the files in `results/`.*\n")

    # --- checklist
    L.append("## 1. Checklist\n")
    L.append("| Check | What it verifies | Result | Measured |")
    L.append("|---|---|---|---|")
    for num, (title, what) in TESTS.items():
        cases = [c for c in d["tests"] if c["number"] == num]
        status = "PASS" if cases and all(c["status"] == "PASS" for c in cases) else "FAIL"
        meas = "; ".join(m for c in cases for m in c["measured"])
        L.append(f"| Test {num}: {title} | {cell(what)} ({len(cases)} case{'s' if len(cases) != 1 else ''}) | "
                 f"**{status}** | {cell(meas)} |")
    vb = val[val.check.str.startswith("conduction velocity") & val.is_base_grid.astype(bool)].iloc[0]
    v_ok = abs(vb.pct_difference) < VELOCITY_TOL_PCT
    L.append(f"| Validation: HH 1952 velocity | θ at {f(vb.T_C, 1)} °C vs HH's computed 18.8 m/s "
             f"(PASS if within ±{VELOCITY_TOL_PCT:g} %, a criterion I chose; the spec only asks for the % difference) | "
             f"**{'PASS' if v_ok else 'FAIL'}** | θ = {f(vb.model_value, 3)} m/s, difference {f(vb.pct_difference, 2)} % "
             f"(vs measured 21.2 m/s: {f(vb.pct_difference_vs_measured, 1)} %) |")
    js = d["js_consistency"]
    js_ok = bool(js["pass"].all())
    L.append(f"| Explorer: JavaScript vs Python | r from the in-browser model within 0.5 % of Python at "
             f"{', '.join(f(t, 1) for t in js.T_C)} °C | **{'PASS' if js_ok else 'FAIL'}** | "
             + "; ".join(f"{f(r.T_C, 1)} °C: r_JS = {f(r.r_js, 6)}, r_Py = {f(r.r_python, 6)} "
                         f"(rel. diff. {g(r.rel_diff_r)})" for r in js.itertuples()) + " |")
    L.append("")

    # --- main table
    L.append("## 2. Results by temperature\n")
    L.append("Primary definitions: Q_Na = ∫(Ĩ_Na − Ĩ_Na,rest) dt over the whole run; membrane stimulus = 2 × threshold. "
             "Travelling spike measured at x = 2.5 cm. Charges in nC/cm². "
             "Overlap A = share of Q_Na that enters after the voltage peak.\n")
    L.append("| T (°C) | φ | r membrane | r travel | Δr | flag | θ (m/s) | Q_Na mem | Q_Na trav | "
             "Q_min mem | Q_min trav | overlap A mem | overlap A trav |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for T, c in cmp_.iterrows():
        m, t = mem.loc[T], trav.loc[T]
        L.append(f"| {f(T, 1)} | {f(c.phi, 3)} | {f(c.r_membrane, 3)} | {f(c.r_travel, 3)} | {f(c.delta_r, 3)} | "
                 f"{c.flag} | {f(t.velocity_m_s, 2)} | {f(m.Q_Na, 1)} | {f(t.Q_Na, 1)} | {f(m.Q_min, 1)} | "
                 f"{f(t.Q_min, 1)} | {f(100 * m.ov_frac_A, 1)} % | {f(100 * t.ov_frac_A, 1)} % |")
    extra = [T for T in mem4.index if T not in cmp_.index]
    if extra:
        L.append("\nAbove that, only the membrane spike at 4 × threshold still meets the spec's spike rule: "
                 + "; ".join(f"{f(T, 1)} °C: r = {f(mem4.loc[T].r, 3)}" for T in extra) + ".")
    L.append(f"\nΔr flag: \"real difference\" if |Δr| > 2ε with ε = sqrt(ε_mem² + ε_cab²) = "
             f"{f(float(cmp_.eps_combined.iloc[0]), 4)}, so 2ε = {f(float(cmp_.two_eps.iloc[0]), 4)}. "
             "Every variant (raw Q_Na, overlap B, 4 × stimulus) is in `results/sweep.csv` and `results/comparison.csv`.\n")

    # --- failure temperatures and epsilon
    L.append("## 3. Failure temperatures\n")
    L.append("| Spike | Rule | Is the spec's rule? | T_fail (°C) |")
    L.append("|---|---|---|---|")
    for r in d["failure_temps"].itertuples():
        L.append(f"| {r.spike} | {r.criterion} | {'yes' if r.is_spec_definition else 'no (diagnostic)'} | "
                 f"{f(r.T_fail_C, 3)} ± {f(r.half_width_C, 3)} (between {f(r.last_success_C, 4)} and {f(r.first_failure_C, 4)}) |")
    L.append("\n**How to read this.** The spec's rules mark the temperature at which the spike *peak* falls below 0 mV. "
             "The heat block in this model is gradual: just above those temperatures the response is still active and "
             "the travelling wave still runs the full length of the axon with a peak below 0 mV. The diagnostic rows "
             "(not part of the spec) mark where the response really disappears. Write \"the spike peak falls below 0 mV "
             "at …\", not \"conduction fails at …\", when quoting the spec values.\n")
    L.append("## 4. Numerical uncertainty ε (grid refinement)\n")
    L.append(f"ε = max |value(refined) − value(base)| over the two refined grids. Computed at 6.3 °C and at "
             f"T_warm = {f(T_warm, 1)} °C (highest coarse T with a travelling spike, {f(cfg['T_highest_propagating_coarse_C'], 1)} °C, minus 2 °C). "
             "The larger value is used at every temperature.\n")
    e = d["epsilon"]
    cols = [c for c in e.columns if c.startswith("eps_at_")]
    L.append("| Model | Quantity | " + " | ".join(c.replace("eps_at_", "ε at ").replace("C", " °C") for c in cols) + " | ε used |")
    L.append("|---|---|" + "---|" * len(cols) + "---|")
    for _, r in e.iterrows():
        L.append(f"| {r.model} | {r.quantity} | " + " | ".join(g(r[c], 2) for c in cols) + f" | {g(r.eps_used, 2)} |")
    L.append("")

    # --- soft expectations
    L.append("## 5. Soft sanity expectations (matches / differs, not pass/fail)\n")
    L.append("| Expectation | Model value | Verdict |")
    L.append("|---|---|---|")
    m63, t63 = mem.loc[P.SWEEP_T_FIRST_C], trav.loc[P.SWEEP_T_FIRST_C]
    vp = 30 <= m63.V_peak <= 50
    L.append(f"| V_peak at 6.3 °C roughly +30 to +50 mV | membrane {f(m63.V_peak, 2)} mV, travelling {f(t63.V_peak, 2)} mV | "
             f"{'matches' if vp else 'differs'} |")
    qm = 80 <= m63.Q_min <= 115
    L.append(f"| Q_min at 6.3 °C about 80–115 nC/cm² | membrane {f(m63.Q_min, 1)}, travelling {f(t63.Q_min, 1)} nC/cm² | "
             f"{'matches' if qm else 'differs'} |")
    order10 = 10 ** 0.5 <= m63.r <= 10 ** 1.5
    L.append(f"| r at 6.3 °C of order 10 (Sengupta et al. 2010: about 11, i.e. about 9 % efficiency, under their own "
             f"definitions) | membrane r = {f(m63.r, 2)} (η = {f(m63.eta_pct, 1)} %), travelling r = {f(t63.r, 2)} "
             f"(η = {f(t63.eta_pct, 1)} %) | {'matches' if order10 else 'differs'} (order of magnitude; numerically "
             f"{f(100 * (m63.r / 11 - 1), 0)} % / {f(100 * (t63.r / 11 - 1), 0)} % above 11) |")
    T18 = 18.0
    seg = mem.loc[[T for T in mem.index if T <= T18]].r.to_numpy()
    seg_t = trav.loc[[T for T in trav.index if T <= T18]].r.to_numpy()
    dec = bool(np.all(np.diff(seg) < 0) and np.all(np.diff(seg_t) < 0))
    L.append(f"| r decreases from 6.3 °C to about 18 °C (Sengupta: efficiency rises to about 24 % at 18 °C) | "
             f"r falls monotonically: membrane {f(seg[0], 2)} → {f(seg[-1], 2)} (η at 18 °C = {f(mem.loc[T18].eta_pct, 1)} %), "
             f"travelling {f(seg_t[0], 2)} → {f(seg_t[-1], 2)} (η = {f(trav.loc[T18].eta_pct, 1)} %) | "
             f"{'matches' if dec else 'differs'} |")
    vok = 16 <= vb.model_value <= 21
    L.append(f"| Velocity at 18.5 °C about 16–21 m/s | {f(vb.model_value, 2)} m/s | {'matches' if vok else 'differs'} |")
    L.append("\nRules used for the verdicts: V_peak and Q_min inside the stated ranges; \"order 10\" = between 10^0.5 and "
             "10^1.5; \"decreases\" = every step on the coarse grid from 6.3 to 18 °C goes down; velocity inside 16–21 m/s.\n")

    # --- answers
    L.append("## 6. Answers to the research question\n")
    L.extend(answers(d, mem, mem4, trav, cmp_))
    L.append("")
    L.extend(paper_values(d, mem, trav))
    return "\n".join(L) + "\n"


CLEARLY_BIGGER = 2.0   # "[A] and [B] clearly bigger than [E]" = more than 2 × [E] (same factor as the 2ε flag)


def _rel_changes(ref):
    """Relative change of r against the base grid, per model / T / level (percent, signed)."""
    r = ref[ref.quantity == "r"].copy()
    base = r[r.level == 0].set_index(["model", "T_C"]).value
    r["rel_pct"] = [100.0 * (v - base[(m, t)]) / base[(m, t)] for m, t, v in zip(r.model, r.T_C, r.value)]
    return r


def _finest_pct(ref, T):
    """(r_travel − r_membrane)/r_membrane × 100 using the finest grid of a refinement table."""
    r = ref[(ref.quantity == "r") & (ref.T_C == T)]
    top = r[r.level == r.level.max()].set_index("model").value
    return 100.0 * (top["travelling"] - top["membrane"]) / top["membrane"]


def _decimals(uncertainty):
    """Decimals whose last digit is no finer than the uncertainty (10^-d >= uncertainty), at least 0."""
    if not np.isfinite(uncertainty) or uncertainty <= 0:
        return 2
    return max(0, int(np.floor(-np.log10(uncertainty))))


def _rounded(x, uncertainty):
    return f"{x:.{_decimals(uncertainty)}f}".replace("-", "−")


def paper_values(d, mem, trav):
    """Section 7: the bracketed values of the paper sentence, each derived from results/."""
    der = d["config"]["derived"]
    near = d["near_failure"].set_index("T_C")
    fg = d["failure_vs_grid"]
    T0, T_edge, T_last = P.SWEEP_T_FIRST_C, der["T_edge_C"], der["T_last_both_C"]
    T_both = [T for T in mem.index if T in trav.index]
    pcts = {T: 100.0 * (trav.loc[T].r - mem.loc[T].r) / mem.loc[T].r for T in T_both}
    A, B = pcts[T0], pcts[T_edge]
    higher = all(trav.loc[T].r > mem.loc[T].r for T in T_both)
    lower = all(trav.loc[T].r < mem.loc[T].r for T in T_both)
    hl = "higher" if higher else "lower" if lower else "mixed (higher at some temperatures, lower at others)"

    ft, fm, fd = fail_row(d, "travelling"), fail_row(d, "membrane_amp2x"), fail_row(d, "travelling", spec=False)
    fg_t = fg[fg.spike == "travelling"].sort_values("level")
    fg_m = fg[fg.spike == "membrane_amp2x"].sort_values("level")
    shift_t = float(fg_t.shift_from_base_C.abs().max())
    shift_m = float(fg_m.shift_from_base_C.abs().max())
    T_fine = float(fg_t.T_fail_C.iloc[-1])
    if shift_t > P.T_FAIL_ABOUT_SHIFT_C:
        T_quote = f"about {T_fine:.0f} °C"
    else:
        T_quote = f"{ft.T_fail_C:.1f} °C"

    rc = _rel_changes(d["refinement"])
    big = rc.loc[rc.rel_pct.abs().idxmax()]
    E_plan = abs(big.rel_pct)
    E_first = rc[rc.level == 1].rel_pct.abs().max()
    E_mem = rc[rc.model == "membrane"].rel_pct.abs().max()
    rce = _rel_changes(d["refinement_edge"])
    big_e = rce.loc[rce.rel_pct.abs().idxmax()]
    E_edge = abs(big_e.rel_pct)
    E = max(E_plan, E_edge)                          # largest change in r over every grid check
    E_quote = f"{E:.{max(0, -int(np.floor(np.log10(E))))}f}%"   # one significant figure
    A_fine, B_fine = _finest_pct(d["refinement"], T0), _finest_pct(d["refinement_edge"], T_edge)
    A_unc, B_unc = abs(A - A_fine), abs(B - B_fine)  # grid sensitivity at that temperature (percentage points)
    A_quote, B_quote = _rounded(A, A_unc) + "%", _rounded(B, B_unc) + "%"
    real = A > CLEARLY_BIGGER * E and B > CLEARLY_BIGGER * E

    # Overlap: r = Q_ov/Q_min + (Q_Na − Q_ov)/Q_min, so Δr splits exactly into an overlap part and the rest.
    ov = {}
    for key, name in (("Q_ov_A", "A: Na⁺ entering after the voltage peak"),
                      ("Q_ov_B", "B: Na⁺ entering while K⁺ flows out at the same time")):
        part = {T: trav.loc[T][key] / trav.loc[T].Q_min - mem.loc[T][key] / mem.loc[T].Q_min for T in T_both}
        expl = {T: 100.0 * part[T] / (trav.loc[T].r - mem.loc[T].r) for T in T_both}
        share = {T: (trav.loc[T][key] / trav.loc[T].Q_Na, mem.loc[T][key] / mem.loc[T].Q_Na) for T in T_both}
        word = "more" if all(v > 0 for v in part.values()) else "less" if all(v < 0 for v in part.values()) else "mixed"
        ov[key] = dict(name=name, part=part, expl=expl, share=share, word=word)

    # Check 1: split the r gap into "more sodium" and "smaller C·ΔV": r_t/r_m = (Q_Na,t/Q_Na,m) / (Q_min,t/Q_min,m)
    split = []
    for T in T_both:
        split.append((T, "sweep", trav.loc[T].Q_Na, mem.loc[T].Q_Na, trav.loc[T].Q_min, mem.loc[T].Q_min))
    for T, row in near.iterrows():
        if T not in T_both:
            split.append((T, "extra run", row.Q_Na_travel, row.Q_Na_membrane, row.Q_min_travel, row.Q_min_membrane))
    split.sort()
    split_rows = []
    for T, src, qt, qm_, mt, mm in split:
        qn, qmin = qt / qm_, mt / mm
        rr = qn / qmin
        split_rows.append(dict(T=T, src=src, gap=100 * (rr - 1), dNa=100 * (qn - 1), dQmin=100 * (qmin - 1),
                               sNa=100 * np.log(qn) / np.log(rr), sQmin=-100 * np.log(qmin) / np.log(rr)))
    sp = {r["T"]: r for r in split_rows}
    edge_from_Na = min(sp[T_edge]["sNa"], sp[T_last]["sNa"])
    steps = [(T_both[i], T_both[i + 1], pcts[T_both[i + 1]] - pcts[T_both[i]]) for i in range(len(T_both) - 1)]
    big_step = max(steps, key=lambda s_: s_[2])
    last_step = steps[-1]
    T_peak = max(pcts, key=pcts.get)
    jump_at_end = last_step[2] >= big_step[2] and last_step[2] > 0

    L = ["## 7. Values for the paper sentence\n",
         "Every value below is computed from `results/` (sweep, refinement and the extra near-failure runs). "
         "Percentage difference = (r_propagating − r_membrane) ÷ r_membrane × 100, primary Q_Na definition, "
         "membrane stimulus 2 × threshold, base grid unless stated.\n",
         "**Rounding rule.** [A] and [B] are rounded so that their last digit is no finer than how much they change "
         "between the base and the finest grid at that temperature. [E] is given to one significant figure. "
         f"[T_fail] is quoted as \"about N °C\" if refining the grid moves it by more than "
         f"{P.T_FAIL_ABOUT_SHIFT_C:g} °C.\n",
         "| Slot | Value to quote | Unrounded | Note |", "|---|---|---|---|",
         f"| [T_fail] | **{T_quote}** | {f(ft.T_fail_C, 3)} °C (base grid) | Planned rule: V at x = 4.0 cm no longer "
         f"exceeds 0 mV. Refining the grid moves it by up to {f(shift_t, 2)} °C (finest grid: {f(T_fine, 2)} °C). The wave "
         f"keeps travelling, with a peak below 0 mV, until {f(fd.T_fail_C, 2)} °C (diagnostic). |",
         f"| [A] | **{A_quote}** | {f(A, 3)} % | at {f(T0, 1)} °C; finest grid {f(A_fine, 2)} % |",
         f"| [B] | **{B_quote}** | {f(B, 3)} % | at {f(T_edge, 1)} °C, the warmest sweep temperature where both impulses exist "
         f"(r cannot be measured at T_fail itself); finest grid {f(B_fine, 2)} % |",
         f"| [higher/lower] | **{hl}** | – | r_propagating vs r_membrane at all {len(T_both)} sweep temperatures with both impulses |",
         f"| [E] | **{E_quote}** | {f(E, 3)} % | largest change in r over every grid check ({big_e.model if E_edge >= E_plan else big.model}, "
         f"{f(T_edge if E_edge >= E_plan else big.T_C, 1)} °C, finest vs base grid) |",
         f"| Real difference? | **{'yes' if real else 'NO'}** | – | [A] = {f(A / E, 1)} × [E], [B] = {f(B / E, 1)} × [E]; "
         f"\"clearly bigger\" taken as more than {CLEARLY_BIGGER:g} × [E] |",
         f"| [more/less] | **{ov['Q_ov_B']['word']}** (definition B), **{ov['Q_ov_A']['word']}** (definition A) | – | "
         f"overlap sodium per unit of capacitor charge (Q_ov/Q_min) in the propagating impulse; how much of the "
         f"difference it explains is below |",
         "",
         "Copy-paste list:\n",
         "```",
         f"T_fail: {T_quote}", f"A: {A_quote}", f"B: {B_quote}", f"higher/lower: {hl}", f"E: {E_quote}",
         f"more/less: {ov['Q_ov_B']['word']}",
         "```",
         ""]
    if not real:
        L.append(f"**There is no clear difference:** [A] and/or [B] are not more than {CLEARLY_BIGGER:g} × [E], so the "
                 "difference between the two impulse types cannot be separated from numerical error. Do not use the "
                 "sentence as written.\n")

    # --- Check 1
    L.append("### Check 1: is [B] an edge effect?\n")
    L.append("Near failure both impulses get smaller, which lowers Q_min = C_m ΔV and raises r even if sodium entry does "
             "not change. Since r_prop/r_mem = (Q_Na,prop/Q_Na,mem) ÷ (Q_min,prop/Q_min,mem), the gap splits exactly "
             "into a \"more sodium\" part and a \"smaller C_m ΔV\" part (shares computed on a log scale).\n")
    L.append("| T (°C) | source | r gap | Q_Na prop vs mem | C_m ΔV prop vs mem | gap from more sodium | gap from smaller C_m ΔV |")
    L.append("|---|---|---|---|---|---|---|")
    for r in split_rows:
        L.append(f"| {f(r['T'], 3 if r['src'] != 'sweep' else 1)} | {r['src']} | {f(r['gap'], 2)} % | {f(r['dNa'], 1)} % | "
                 f"{f(r['dQmin'], 1)} % | {f(r['sNa'], 0)} % | {f(r['sQmin'], 0)} % |")
    L.append("")
    if edge_from_Na >= 50:
        L.append(f"- **Not an edge effect.** Near T_fail the gap comes from extra sodium ({f(sp[T_edge]['sNa'], 0)} % at "
                 f"{f(T_edge, 1)} °C, {f(sp[T_last]['sNa'], 0)} % at {f(T_last, 3)} °C), not from a smaller C_m ΔV: close "
                 "to failure both impulses shrink by about the same amount. The smaller-C_m ΔV part is largest at the cold "
                 f"end ({f(sp[T0]['sQmin'], 0)} % of the gap at {f(T0, 1)} °C).")
    else:
        L.append(f"- **Partly an edge effect.** Near T_fail only {f(edge_from_Na, 0)} % of the gap comes from extra sodium; "
                 "the rest comes from a smaller C_m ΔV, so do not attribute the whole gap to overlap.")
    if jump_at_end:
        L.append(f"- **The gap jumps at the end** ({f(last_step[2], 2)} percentage points between {f(last_step[0], 1)} and "
                 f"{f(last_step[1], 1)} °C, the largest step). Consider quoting {f(pcts[18.5], 2)} % at 18.5 °C instead.")
    else:
        L.append(f"- **The trend is smooth.** The gap rises from {f(A, 2)} % at {f(T0, 1)} °C to a maximum of "
                 f"{f(pcts[T_peak], 2)} % at {f(T_peak, 1)} °C (largest single step {f(big_step[2], 2)} percentage points, "
                 f"{f(big_step[0], 1)} → {f(big_step[1], 1)} °C), then changes by {f(last_step[2], 2)} points to "
                 f"{f(B, 2)} % at {f(T_edge, 1)} °C and is {f(near.loc[T_last].pct_difference, 2)} % at {f(T_last, 3)} °C. "
                 f"No jump at the end. For reference, the value at 18.5 °C is {f(pcts[18.5], 2)} % (not grid-checked at that "
                 "temperature).")
    L.append("")

    # --- Check 2
    L.append("### Check 2: where was [E] measured, and does T_fail move with the grid?\n")
    L.append(f"- The planned check ran at 6.3 °C and T_warm = {f(der['T_warm_C'], 1)} °C: largest change in r "
             f"{f(E_plan, 2)} % ({big.model}, {f(big.T_C, 1)} °C); after the first halving only {f(E_first, 2)} %; membrane "
             f"patch at most {g(E_mem, 1)} %.")
    L.append(f"- Repeated at {f(T_edge, 1)} °C, next to T_fail (extra check): largest change {f(E_edge, 2)} % ({big_e.model}). "
             f"[E] above uses the larger of the two, {f(E, 2)} %.")
    L.append(f"- Failure temperature re-found on each grid (bisection to {P.FAILURE_GRID_TOL_C:g} °C):\n")
    L.append("| Impulse | Grid | T_fail (°C) | shift from base grid (°C) |")
    L.append("|---|---|---|---|")
    for r in fg.sort_values(["spike", "level"]).itertuples():
        grid = (f"Δx = {r.dx_um:g} µm, Δt = {r.dt_us:g} µs" if r.spike == "travelling" else f"Δt = {r.dt_us:g} µs")
        L.append(f"| {r.spike} | {grid} | {f(r.T_fail_C, 3)} ± {f(r.half_width_C, 3)} | {f(r.shift_from_base_C, 3)} |")
    steps_t = np.diff(fg_t.T_fail_C.to_numpy())
    L.append("")
    if shift_t > P.T_FAIL_ABOUT_SHIFT_C:
        L.append(f"- The propagating impulse's failure temperature moves by up to {f(shift_t, 2)} °C (more than "
                 f"{P.T_FAIL_ABOUT_SHIFT_C:g} °C), so quote **{T_quote}**, not {f(ft.T_fail_C, 1)} °C. The successive shifts "
                 f"({', '.join(f(s_, 3) for s_ in steps_t)} °C) get smaller, so the value is settling as the grid is refined.")
    else:
        L.append(f"- The propagating impulse's failure temperature moves by at most {f(shift_t, 3)} °C, so {T_quote} is fine.")
    L.append(f"- The membrane impulse's failure temperature does not move at this resolution "
             f"(largest shift {f(shift_m, 3)} °C, bisection tolerance {P.FAILURE_GRID_TOL_C:g} °C)."
             if shift_m < P.FAILURE_GRID_TOL_C else
             f"- The membrane impulse's failure temperature moves by at most {f(shift_m, 3)} °C.")
    L.append("")

    # --- [more/less]
    L.append("### [more/less]: does overlap explain the difference?\n")
    L.append("Since r = Q_ov/Q_min + (Q_Na − Q_ov)/Q_min, the difference in r splits exactly into an overlap part and "
             "the rest.\n")
    L.append("| Overlap definition | T (°C) | Q_ov membrane | Q_ov propagating | share of Q_Na (membrane / propagating) | "
             "Q_ov/Q_min (membrane / propagating) | part of the r difference explained |")
    L.append("|---|---|---|---|---|---|---|")
    for key in ("Q_ov_B", "Q_ov_A"):
        o = ov[key]
        for T in (T0, T_edge):
            m, t = mem.loc[T], trav.loc[T]
            L.append(f"| {cell(o['name'])} | {f(T, 1)} | {f(m[key], 1)} nC/cm² | {f(t[key], 1)} nC/cm² | "
                     f"{f(100 * o['share'][T][1], 1)} % / {f(100 * o['share'][T][0], 1)} % | "
                     f"{f(m[key] / m.Q_min, 2)} / {f(t[key] / t.Q_min, 2)} | {f(o['expl'][T], 0)} % |")
    L.append("")
    for key in ("Q_ov_B", "Q_ov_A"):
        o = ov[key]
        e0, e1 = o["expl"][T0], o["expl"][T_edge]
        shares_bigger = [o["share"][T][0] > o["share"][T][1] for T in T_both]
        share_txt = ("bigger at every temperature" if all(shares_bigger) else
                     "smaller at every temperature" if not any(shares_bigger) else "bigger only at some temperatures")
        if e0 >= 50 and e1 >= 50:
            verdict = "overlap explains most of the difference at both ends, so the sentence is supported"
        elif e0 < 50 and e1 < 50:
            verdict = "overlap explains less than half of the difference at both ends: delete the sentence rather than claim it"
        else:
            verdict = (f"overlap explains most of the difference at one end but only {f(min(e0, e1), 0)} % at the other, "
                       "so say it explains the difference partly, or only at the temperature where it does")
        L.append(f"- Definition {o['name']}: the propagating impulse has **{o['word']}** overlap sodium per unit of "
                 f"capacitor charge; it explains {f(e0, 0)} % of the r difference at {f(T0, 1)} °C and {f(e1, 0)} % at "
                 f"{f(T_edge, 1)} °C; {verdict}. As a share of Q_Na, the propagating impulse's overlap is {share_txt}, so "
                 "write \"more overlap sodium\", not \"a larger overlap fraction\".")
    return L


def answers(d, mem, mem4, trav, cmp_):
    T0, T1 = mem.index.min(), mem.index.max()
    Tt1 = trav.index.max()
    real = cmp_.flag == "real difference"
    rel0 = 100 * cmp_.delta_r.iloc[0] / cmp_.r_membrane.iloc[0]
    rel1 = 100 * cmp_.delta_r.iloc[-1] / cmp_.r_membrane.iloc[-1]
    ft_m, ft_t = fail_row(d, "membrane_amp2x"), fail_row(d, "travelling")
    dg_t = fail_row(d, "travelling", spec=False)
    s = []
    s.append(f"1. **r falls steeply as the axon warms.** For the membrane spike r drops from {f(mem.loc[T0].r, 2)} at "
             f"{f(T0, 1)} °C to {f(mem.loc[T1].r, 2)} at {f(T1, 1)} °C (efficiency {f(mem.loc[T0].eta_pct, 1)} % → "
             f"{f(mem.loc[T1].eta_pct, 1)} %); for the travelling spike from {f(trav.loc[T0].r, 2)} to "
             f"{f(trav.loc[Tt1].r, 2)} at {f(Tt1, 1)} °C.")
    s.append(f"2. **The fall in r comes from Q_Na, not from Q_min.** For the membrane spike Q_Na falls "
             f"{f(mem.loc[T0].Q_Na / mem.loc[T1].Q_Na, 1)}-fold ({f(mem.loc[T0].Q_Na, 0)} → {f(mem.loc[T1].Q_Na, 0)} nC/cm²) "
             f"while Q_min falls only {f(mem.loc[T0].Q_min / mem.loc[T1].Q_min, 2)}-fold ({f(mem.loc[T0].Q_min, 0)} → "
             f"{f(mem.loc[T1].Q_min, 0)} nC/cm², which on its own would raise r); the share of Na⁺ entering after the "
             f"peak (overlap A) drops from {f(100 * mem.loc[T0].ov_frac_A, 0)} % to {f(100 * mem.loc[T1].ov_frac_A, 0)} %.")
    if real.all():
        cmp_txt = (f"at every temperature from {f(cmp_.index.min(), 1)} to {f(cmp_.index.max(), 1)} °C")
    else:
        cmp_txt = (f"at {int(real.sum())} of {len(real)} temperatures ("
                   + ", ".join(f(T, 1) for T in cmp_.index[real]) + " °C)")
    s.append(f"3. **The membrane patch does not reproduce the travelling spike beyond numerical uncertainty.** "
             f"r_travel − r_membrane is positive {cmp_txt}: Δr = {f(cmp_.delta_r.iloc[0], 3)} at "
             f"{f(cmp_.index[0], 1)} °C and {f(cmp_.delta_r.iloc[-1], 3)} at {f(cmp_.index[-1], 1)} °C, against "
             f"2ε = {f(float(cmp_.two_eps.iloc[0]), 3)}; the patch underestimates r by {f(rel0, 1)} % at "
             f"{f(cmp_.index[0], 1)} °C and {f(rel1, 1)} % at {f(cmp_.index[-1], 1)} °C.")
    s.append(f"4. **The size of that gap depends on the patch stimulus.** With 4 × threshold instead of 2 ×, "
             f"Δr is {f(cmp_.delta_r_amp4x.min(), 3)}–{f(cmp_.delta_r_amp4x.max(), 3)} over the same range, so the "
             f"membrane estimate is not unique.")
    s.append(f"5. **Warm limit.** The spike peak falls below 0 mV at {f(ft_m.T_fail_C, 2)} °C (membrane, 2 ×) and "
             f"{f(ft_t.T_fail_C, 2)} °C (travelling, spec rule), but the travelling wave keeps running without decaying "
             f"until {f(dg_t.T_fail_C, 2)} °C (diagnostic rule).")
    return s


# ---------------------------------------------------------------------------
# METHODS_FILL_IN.md
# ---------------------------------------------------------------------------
def methods_fill_in(d):
    cfg, env, val = d["config"], d["env"], d["validation"]
    der = cfg["derived"]
    T_warm = der["T_warm_C"]
    js = d["js_consistency"]
    L = ["# Values for the Methods section\n",
         f"*Generated automatically by `run_all.py` on {datetime.datetime.now():%Y-%m-%d %H:%M}. "
         "Every value is read from `results/`.*\n",
         "## Computer and software\n",
         f"- Computer / machine model: {env.get('Machine model (sysctl hw.model)', 'n/a')}; "
         f"CPU: {env.get('CPU (sysctl machdep.cpu.brand_string)', 'n/a')} ({env.get('Machine (platform.machine)')})",
         f"- Operating system: {env.get('OS')}",
         f"- Python {env.get('Python')}; numpy {env.get('numpy')}; scipy {env.get('scipy')}; "
         f"matplotlib {env.get('matplotlib')}; pandas {env.get('pandas')}",
         f"- Node.js (explorer consistency check only): {js.node_version.iloc[0]}",
         "",
         "## Model and numerics\n",
         f"- Membrane: C_m = {cfg['membrane']['C_m']:g} µF/cm², ḡ_Na = {cfg['membrane']['g_Na']:g}, ḡ_K = "
         f"{cfg['membrane']['g_K']:g}, g_L = {cfg['membrane']['g_L']:g} mS/cm², E_Na = {cfg['membrane']['E_Na']:g}, "
         f"E_K = {cfg['membrane']['E_K']:g}, E_L = {cfg['membrane']['E_L']:g} mV (NEURON hh.mod rate functions).",
         f"- Temperature factor φ(T) = {cfg['Q10']:g}^((T − {cfg['T_REF_C']:g})/10), applied to the gating rates only.",
         f"- Resting potential V_rest = {f(der['rest']['V'], 4)} mV (root of the total steady-state ionic current).",
         f"- Membrane patch: classical RK4, Δt = {cfg['DT_MEMBRANE_MS'] * 1000:g} µs, stimulus integrated piece by piece; "
         f"reference solver scipy DOP853 (rtol {cfg['REF_RTOL']:g}, atol {cfg['REF_ATOL']:g}).",
         f"- Cable: a = {cfg['RADIUS_CM'] * 1e4:g} µm, R_i = {cfg['R_I_OHM_CM']:g} Ω·cm, L = {cfg['LENGTH_CM']:g} cm, sealed ends, "
         f"vertex-centred grid with Δx = {cfg['DX_CABLE_CM'] * 1e4:g} µm ({der['n_nodes_base']} nodes), Δt = "
         f"{cfg['DT_CABLE_MS'] * 1000:g} µs; gates by exponential Euler, voltage by backward Euler "
         f"(tridiagonal solve). Coupling coefficient a/(2R_i) = {f(cfg['D_coupling_uA_per_mV'], 4)} µA/cm² per (mV/cm²). "
         f"Explicit stability number DΔt/(C_mΔx²) = {f(der['explicit_stability_number'], 1)} (explicit stepping would need ≤ 0.5).",
         f"- Travelling-spike stimulus: {cfg['CABLE_STIM_AMP_UA']:g} µA/cm² for {cfg['STIM_DURATION_MS']:g} ms from "
         f"t = {cfg['STIM_START_MS']:g} ms on all nodes with x < {cfg['CABLE_STIM_X_MAX_CM']:g} cm; measured at "
         f"x = {cfg['MEASURE_X_CM']:g} cm; velocity from the 0 mV crossings at x = 2.0 and 3.0 cm.",
         f"- Run length: until V (membrane, or the cable node at 2.5 cm) is back within "
         f"{cfg['REST_RETURN_TOL_MV']:g} mV of V_rest after the after-hyperpolarisation (cap {cfg['T_MAX_MS']:g} ms). "
         "Charges integrated with the trapezoid rule.",
         "",
         "## Stimulus amplitude (membrane spike)\n",
         f"- Threshold at 6.3 °C (0.3 ms pulse, bisection to {cfg['THRESHOLD_TOL_UA']:g} µA/cm²): "
         f"**{f(der['threshold_uA_cm2'], 3)} µA/cm²** (bracket {f(der['threshold_bracket_uA_cm2'][0], 3)}–"
         f"{f(der['threshold_bracket_uA_cm2'][1], 3)})",
         f"- Primary amplitude (2 × threshold): **{f(der['amp_primary_uA_cm2'], 3)} µA/cm²**; "
         f"variant (4 ×): **{f(der['amp_variant_uA_cm2'], 3)} µA/cm²**; same amplitude at every temperature.",
         "",
         "## RK4 vs DOP853\n"]
    for T in sorted(val.T_C.unique()):
        v = val[(val.T_C == T) & val.check.str.startswith("RK4")].set_index("check")["model_value"]
        if v.empty:
            continue
        L.append(f"- {f(T, 1)} °C: max |V_RK4 − V_DOP853| = **{g(v['RK4 vs DOP853: max |dV| over the whole run'])} mV** "
                 f"(during the spike: {g(v['RK4 vs DOP853: max |dV| during the spike (raw window)'])} mV); "
                 f"relative difference in Q_Na = **{g(v['RK4 vs DOP853: relative difference in Q_Na (primary)'])}**")
    L += ["", "## Numerical uncertainty ε (grid refinement)\n",
          f"- Cable grids (Δx, Δt): " + ", ".join(f"({dx * 1e4:g} µm, {dt * 1000:g} µs)" for dx, dt in P.CABLE_REFINEMENT)
          + "; membrane Δt: " + ", ".join(f"{dt * 1000:g} µs" for dt in P.MEMBRANE_REFINEMENT) + ".",
          f"- Refinement temperatures: 6.3 °C and T_warm = {f(T_warm, 1)} °C."]
    for model, q, unit in (("membrane", "r", ""), ("travelling", "r", ""), ("travelling", "velocity_m_s", " m/s"),
                           ("membrane", "Q_Na", " nC/cm²"), ("travelling", "Q_Na", " nC/cm²"),
                           ("membrane", "Q_min", " nC/cm²"), ("travelling", "Q_min", " nC/cm²")):
        L.append(f"- ε for {q.replace('velocity_m_s', 'θ')} ({model}): {g(eps(d, model, q, 'eps_at_6.3C'))}{unit} at 6.3 °C, "
                 f"{g(eps(d, model, q, f'eps_at_{T_warm:g}C'))}{unit} at {f(T_warm, 1)} °C → used: "
                 f"**{g(eps(d, model, q))}{unit}**")
    c = d["comparison"]
    rce = _rel_changes(d["refinement_edge"])
    big_e = rce.loc[rce.rel_pct.abs().idxmax()]
    L.append(f"- Extra check, not part of the planned method: the same refinement at {f(der['T_edge_C'], 1)} °C (warmest "
             f"sweep temperature with both impulses) changes r by at most {f(abs(big_e.rel_pct), 2)} % ({big_e.model}).")
    L.append(f"- Combined ε for Δr = sqrt(ε_mem² + ε_cab²) = **{f(float(c.eps_combined.iloc[0]), 4)}**; "
             f"flag threshold 2ε = {f(float(c.two_eps.iloc[0]), 4)}.")
    L += ["", "## Temperature grid\n",
          f"- Coarse sweep: {', '.join(f(t, 1) for t in der['temperatures_C'])} °C "
          f"(6.3 °C, then 2 °C steps from 8 °C, plus 18.5 °C, until every spike type had failed)."]
    L += ["", "## Failure temperatures (bisection to 0.05 °C)\n"]
    for r in d["failure_temps"].itertuples():
        tag = "" if r.is_spec_definition else " — DIAGNOSTIC, not the spec's rule"
        L.append(f"- {r.spike}: **{f(r.T_fail_C, 3)} ± {f(r.half_width_C, 3)} °C** ({r.criterion}){tag}")
    fg = d["failure_vs_grid"]
    for spike in ("travelling", "membrane_amp2x"):
        rows = fg[fg.spike == spike].sort_values("level")
        L.append(f"- Extra check, not part of the planned method: {spike} failure temperature re-found on each refinement "
                 f"grid (bisection to {P.FAILURE_GRID_TOL_C:g} °C): " + ", ".join(f"{f(T, 3)} °C" for T in rows.T_fail_C)
                 + f" (largest shift from the base grid {f(rows.shift_from_base_C.abs().max(), 3)} °C)")
    L += ["", "## Validation against Hodgkin & Huxley (1952)\n"]
    for r in val[val.check.str.startswith("conduction velocity")].itertuples():
        L.append(f"- {r.grid}: θ(18.5 °C) = **{f(r.model_value, 3)} m/s**; vs HH computed 18.8 m/s: "
                 f"{f(r.pct_difference, 2)} %; vs HH measured 21.2 m/s: {f(r.pct_difference_vs_measured, 2)} %"
                 + (" (base grid)" if r.is_base_grid else ""))
    L += ["", "## Explorer consistency check\n"]
    for r in js.itertuples():
        L.append(f"- {f(r.T_C, 1)} °C: r (JavaScript) = {f(r.r_js, 6)}, r (Python) = {f(r.r_python, 6)}, "
                 f"relative difference {g(r.rel_diff_r)} (tolerance {r.tolerance * 100:g} %)")
    return "\n".join(L) + "\n"


def write_all(log=print):
    d = load()
    (ROOT / "RESULTS_SUMMARY.md").write_text(results_summary(d))
    (ROOT / "METHODS_FILL_IN.md").write_text(methods_fill_in(d))
    log("  wrote RESULTS_SUMMARY.md and METHODS_FILL_IN.md")
