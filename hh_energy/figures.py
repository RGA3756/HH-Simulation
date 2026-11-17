"""
All paper figures (PNG at 300 dpi + PDF) from the files in results/.

Every plotted number comes from results/*.csv or results/traces.npz, which
sweep.run_experiments() wrote. Units: V mV, t ms, I µA/cm², Q nC/cm², T °C.
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from . import params as P
from .sweep import load_traces

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
CM = 1 / 2.54

C_NA, C_K = P.COLOR_NA, P.COLOR_K
C_INK = "#1a1a1a"          # main neutral ink (membrane / travelling lines)
C_GREY = "#8a8a8a"         # secondary series (variants)
C_ACCENT = "#a8326e"       # reference values (DOP853, HH 1952, overlap variant A)
C_BAND = "#d9d9d9"         # uncertainty bands and overlap shading
# Voltage colour map: one hue (violet), light = rest/AHP, dark = spike peak
V_CMAP = LinearSegmentedColormap.from_list("v_violet", ["#f3f1fa", "#b9b1de", "#6a58b8", "#2a1866"])
V_RANGE = (-80.0, 50.0)

MEM = dict(ls="-", marker="o", ms=3.5, mfc=None)          # membrane: solid, filled markers
TRAV = dict(ls="--", marker="o", ms=3.5, mfc="white")     # travelling: dashed, open markers


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": P.FONT_SIZE_PT, "axes.labelsize": P.FONT_SIZE_PT,
        "xtick.labelsize": P.FONT_SIZE_PT, "ytick.labelsize": P.FONT_SIZE_PT,
        "legend.fontsize": P.FONT_SIZE_PT, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "lines.linewidth": 1.3, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": P.FIG_DPI, "figure.dpi": 100,
    })


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{name}.png", dpi=P.FIG_DPI)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)


def letter(ax, s, dx=-30, dy=4):
    ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                fontsize=10, fontweight="bold", ha="left", va="bottom")


def load():
    d = {name: pd.read_csv(RESULTS / f"{name}.csv")
         for name in ("sweep", "comparison", "epsilon", "failure_temps", "validation", "refinement")}
    d["config"] = json.loads((RESULTS / "config.json").read_text())
    d["traces"] = load_traces(RESULTS / "traces.npz")
    return d


def pick(sweep, spike_type, variant, qna_def="primary"):
    s = sweep[(sweep.spike_type == spike_type) & (sweep.stim_variant == variant)
              & (sweep.qna_def == qna_def) & sweep.spiked.astype(bool)]
    return s.sort_values("T_C")


def eps_used(eps, model, quantity):
    return float(eps[(eps.model == model) & (eps.quantity == quantity)].eps_used.iloc[0])


def spec_fail(fail, spike):
    row = fail[(fail.spike == spike) & fail.is_spec_definition.astype(bool)]
    return float(row.T_fail_C.iloc[0])


# ---------------------------------------------------------------------------
# 1. Validation: RK4 vs DOP853
# ---------------------------------------------------------------------------
def fig_validation(d):
    v = d["traces"]["val"][6.3]
    t, V1, V2 = v["t"], v["V_rk4"], v["V_ref"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(P.FIG_WIDTH_FULL_CM * CM, 6.0 * CM), layout="constrained")
    a.plot(t, V1, color=C_INK, lw=2.0, label="RK4, Δt = 5 µs")
    a.plot(t, V2, color=C_ACCENT, lw=1.0, ls="--", label="DOP853 reference")
    a.set_xlabel("Time (ms)")
    a.set_ylabel("Membrane potential V (mV)")
    a.legend(loc="upper right")
    diff = np.abs(V1 - V2)
    m = diff > 0
    b.semilogy(t[m], diff[m], color=C_INK, lw=0.9)
    b.axhline(0.01, color=C_GREY, lw=0.8)
    b.text(t[-1], 0.01 * 1.3, "test limit 0.01 mV", color=C_GREY, ha="right", va="bottom")
    b.set_ylim(1e-12, 0.1)
    b.set_xlabel("Time (ms)")
    b.set_ylabel("|V$_{RK4}$ − V$_{DOP853}$| (mV)")
    letter(a, "A")
    letter(b, "B")
    save(fig, "fig_validation")


# ---------------------------------------------------------------------------
# 2. Spike anatomy (teaching figure), membrane spike at 6.3 °C
# ---------------------------------------------------------------------------
def fig_spike_anatomy(d):
    tr = d["traces"]["mem"][6.3]["amp2x"]
    row = pick(d["sweep"], "membrane", "amp2x").set_index("T_C").loc[6.3]
    rest = d["config"]["derived"]["rest"]
    t, V = tr["t"], tr["V"]
    na_in = -tr["I_Na"]
    I_K = tr["I_K"]
    win = (0.0, 7.0)
    k = (t >= win[0]) & (t <= win[1])

    fig, axs = plt.subplots(3, 1, figsize=(P.FIG_WIDTH_FULL_CM * CM, 14.0 * CM), sharex=True,
                            layout="constrained", gridspec_kw={"height_ratios": [1, 0.8, 1.2]})
    a, b, c = axs
    # A: voltage
    a.plot(t[k], V[k], color=C_INK)
    a.axhline(rest["V"], color=C_GREY, lw=0.6)
    a.text(win[1], rest["V"] + 2, f"V$_{{rest}}$ = {rest['V']:.1f} mV".replace("-", "−"), color=C_GREY, ha="right", va="bottom")
    a.annotate(f"V$_{{peak}}$ = {row.V_peak:.1f} mV", xy=(row.t_peak, row.V_peak),
               xytext=(row.t_peak + 0.5, row.V_peak - 8), arrowprops=dict(arrowstyle="-", lw=0.6, color=C_INK),
               va="center")
    x_br = 5.0
    a.plot([row.t_peak, x_br + 0.15], [row.V_peak, row.V_peak], color=C_GREY, lw=0.6, ls=":")
    a.annotate("", xy=(x_br, row.V_peak), xytext=(x_br, rest["V"]),
               arrowprops=dict(arrowstyle="<->", lw=0.7, color=C_INK, shrinkA=0, shrinkB=0))
    a.text(x_br + 0.12, 0.5 * (row.V_peak + rest["V"]),
           "Q$_{min}$ = C$_m$(V$_{peak}$ − V$_{rest}$)\n" + f"= {row.Q_min:.1f} nC/cm²",
           ha="left", va="center")
    a.set_ylabel("V (mV)")
    a.set_ylim(-85, 55)

    # B: gates
    b.plot(t[k], tr["m"][k], color=C_NA, label="m (Na$^+$ activation)")
    b.plot(t[k], tr["h"][k], color=C_NA, ls="--", label="h (Na$^+$ inactivation)")
    b.plot(t[k], tr["n"][k], color=C_K, label="n (K$^+$ activation)")
    b.set_ylabel("Gate value")
    b.set_ylim(0, 1.32)
    b.set_yticks([0, 0.5, 1.0])
    b.legend(loc="upper center", ncols=3, handlelength=1.6, columnspacing=1.2)

    # C: currents
    dNa = na_in - (-rest["I_Na"])
    dK = I_K - rest["I_K"]
    Qcum = np.concatenate([[0.0], np.cumsum(0.5 * (dNa[1:] + dNa[:-1]) * np.diff(t))])
    t_qmin = np.interp(row.Q_min, Qcum, t)          # time at which Q_min of sodium has entered
    kq = t <= t_qmin
    overlap = np.minimum(np.maximum(dNa, 0), np.maximum(dK, 0))
    c.fill_between(t[k], 0, overlap[k], color=C_BAND, lw=0, label="overlap (Na$^+$ in and K$^+$ out together)")
    c.fill_between(t[kq], 0, na_in[kq], color=C_NA, alpha=0.55, lw=0, label="first Q$_{min}$ of Na$^+$ (fills the capacitor)")
    c.plot(t[k], na_in[k], color=C_NA, label="Na$^+$ current entering, Ĩ$_{Na}$")
    c.plot(t[k], I_K[k], color=C_K, label="K$^+$ current leaving, I$_K$")
    for ax in axs:
        ax.axvline(row.t_peak, color=C_GREY, lw=0.6, ls=":")
    c.set_ylabel("Current (µA/cm²)")
    c.set_xlabel("Time (ms)")
    c.set_xlim(*win)
    c.legend(loc="upper right")
    txt = (f"Q$_{{Na}}$ = {row.Q_Na:.0f} nC/cm²\nQ$_{{min}}$ = {row.Q_min:.0f} nC/cm²\n"
           f"r = Q$_{{Na}}$/Q$_{{min}}$ = {row.r:.1f}\n"
           f"Na$^+$ entering after the peak: {100 * row.ov_frac_A:.0f}%")
    c.text(0.985, 0.42, txt, transform=c.transAxes, ha="right", va="top")
    a.text(row.t_peak + 0.05, -82, "peak", color=C_GREY, ha="left", va="bottom")
    for ax, s in zip(axs, "ABC"):
        letter(ax, s, dx=-40)
    save(fig, "fig_spike_anatomy")


# ---------------------------------------------------------------------------
# 3. Currents at three temperatures, time normalised to spike width
# ---------------------------------------------------------------------------
def three_temperatures(d):
    T_warm = d["config"]["derived"]["T_warm_C"]
    coarse = np.array(d["config"]["derived"]["temperatures_C"])
    mid_target = 0.5 * (P.SWEEP_T_FIRST_C + T_warm)
    T_mid = float(coarse[np.argmin(np.abs(coarse - mid_target))])
    return (P.SWEEP_T_FIRST_C, T_mid, T_warm)


def fig_temperature_currents(d):
    temps = three_temperatures(d)
    rest = d["config"]["derived"]["rest"]
    sw = d["sweep"]
    fig, axs = plt.subplots(2, 3, figsize=(P.FIG_WIDTH_FULL_CM * CM, 10.0 * CM), sharex=True, sharey=True,
                            layout="constrained")
    for j, T in enumerate(temps):
        for i, kind in enumerate(("membrane", "travelling")):
            ax = axs[i, j]
            if kind == "membrane":
                tr = d["traces"]["mem"][T]["amp2x"]
                t, V, I_Na, I_K = tr["t"], tr["V"], tr["I_Na"], tr["I_K"]
                row = pick(sw, "membrane", "amp2x").set_index("T_C").loc[T]
                ls = "-"
            else:
                c = d["traces"]["cab"][T]
                t, V, I_Na, I_K = c["t"], c["V"][2.5], c["I_Na"][2.5], c["I_K"][2.5]
                row = pick(sw, "travelling", "cable").set_index("T_C").loc[T]
                ls = "--"
            s = (t - row.t_peak) / row.width
            dNa = -I_Na + rest["I_Na"]
            dK = I_K - rest["I_K"]
            ov = np.minimum(np.maximum(dNa, 0), np.maximum(dK, 0))
            ax.fill_between(s, 0, ov, color=C_BAND, lw=0)
            ax.plot(s, -I_Na, color=C_NA, ls=ls)
            ax.plot(s, I_K, color=C_K, ls=ls)
            ax.axvline(0, color=C_GREY, lw=0.6, ls=":")
            ax.text(0.98, 0.72, f"{T:g} °C\nwidth {row.width:.2f} ms\noverlap A {100 * row.ov_frac_A:.0f}%\n"
                    f"overlap B {100 * row.ov_frac_B:.0f}%",
                    transform=ax.transAxes, ha="right", va="top", linespacing=1.25)
            ax.set_xlim(-1.5, 5.0)
            if j == 0:
                ax.set_ylabel(("Membrane spike" if kind == "membrane" else "Travelling spike")
                              + "\ncurrent (µA/cm²)")
            if i == 1:
                ax.set_xlabel("(t − t$_{peak}$) / spike width")
    handles = [plt.Line2D([], [], color=C_NA, label="Na$^+$ entering, Ĩ$_{Na}$"),
               plt.Line2D([], [], color=C_K, label="K$^+$ leaving, I$_K$"),
               matplotlib.patches.Patch(color=C_BAND, label="overlap (both flowing)"),
               plt.Line2D([], [], color=C_INK, label="membrane spike (solid)"),
               plt.Line2D([], [], color=C_INK, ls="--", label="travelling spike (dashed)")]
    fig.legend(handles=handles, loc="outside lower center", ncols=3)
    for ax, s in zip(axs.flat, "ABCDEF"):
        letter(ax, s, dx=-14 if ax not in axs[:, 0] else -44)
    save(fig, "fig_temperature_currents")


# ---------------------------------------------------------------------------
# 4. r vs T, Q_Na and Q_min vs T
# ---------------------------------------------------------------------------
def fail_lines(ax, fail, label=True):
    specs = [("membrane_amp2x", "-", "membrane 2×"), ("travelling", "--", "travelling"),
             ("membrane_amp4x", ":", "membrane 4×")]
    for spike, ls, name in specs:
        T = spec_fail(fail, spike)
        ax.axvline(T, color=C_GREY, lw=0.7, ls=ls)


def fail_labels(ax, fail):
    """Rotated labels next to the failure lines (2× and travelling lie 0.06 °C apart: one label)."""
    T2, Tt, T4 = (spec_fail(fail, s) for s in ("membrane_amp2x", "travelling", "membrane_amp4x"))
    y = 0.97
    ax.text(max(T2, Tt) + 0.3, y, "membrane 2× and travelling:\npeak < 0 mV beyond here", rotation=90,
            transform=ax.get_xaxis_transform(), ha="left", va="top", color=C_GREY,
            fontsize=P.FONT_SIZE_PT - 1)
    ax.text(T4 + 0.3, y, "membrane 4×:\npeak < 0 mV", rotation=90, transform=ax.get_xaxis_transform(),
            ha="left", va="top", color=C_GREY, fontsize=P.FONT_SIZE_PT - 1)


def fig_r_vs_T(d):
    sw, eps, fail = d["sweep"], d["epsilon"], d["failure_temps"]
    mem, mem4, trav = pick(sw, "membrane", "amp2x"), pick(sw, "membrane", "amp4x"), pick(sw, "travelling", "cable")
    e_m, e_c = eps_used(eps, "membrane", "r"), eps_used(eps, "travelling", "r")
    fig, (a, b) = plt.subplots(1, 2, figsize=(P.FIG_WIDTH_FULL_CM * CM, 7.0 * CM), layout="constrained")

    a.fill_between(mem.T_C, mem.r - e_m, mem.r + e_m, color=C_BAND, lw=0)
    a.fill_between(trav.T_C, trav.r - e_c, trav.r + e_c, color=C_BAND, lw=0)
    a.plot(mem4.T_C, mem4.r, color=C_GREY, lw=0.9, ls="-", marker="o", ms=2.5, label="membrane, 4× stimulus")
    a.plot(mem.T_C, mem.r, color=C_INK, **MEM, label="membrane, 2× stimulus")
    a.plot(trav.T_C, trav.r, color=C_INK, **TRAV, label="travelling")
    fail_lines(a, fail)
    a.set_xlabel("Temperature (°C)")
    a.set_ylabel("r = Q$_{Na}$ / Q$_{min}$")
    a.set_ylim(0, None)
    fail_labels(a, fail)

    b.plot(mem.T_C, mem.Q_Na, color=C_NA, **MEM)
    b.plot(trav.T_C, trav.Q_Na, color=C_NA, **TRAV)
    b.plot(mem.T_C, mem.Q_min, color=C_INK, **MEM)
    b.plot(trav.T_C, trav.Q_min, color=C_INK, **TRAV)
    b.text(8.5, 1700, "Q$_{Na}$ (sodium that entered)", ha="left", va="bottom")
    b.text(8.5, 122, "Q$_{min}$ (capacitor charge)", ha="left", va="bottom")
    fail_lines(b, fail)
    b.set_yscale("log")
    b.set_xlabel("Temperature (°C)")
    b.set_ylabel("Charge (nC/cm²)")
    b.set_ylim(50, 2500)
    fail_labels(b, fail)
    fig.legend(*a.get_legend_handles_labels(), loc="outside lower center", ncols=3)
    for ax in (a, b):
        ax.set_xlim(5, 36)
    letter(a, "A")
    letter(b, "B", dx=-40)
    save(fig, "fig_r_vs_T")


# ---------------------------------------------------------------------------
# 5. Velocity
# ---------------------------------------------------------------------------
def fig_velocity(d):
    trav = pick(d["sweep"], "travelling", "cable")
    e = eps_used(d["epsilon"], "travelling", "velocity_m_s")
    fig, ax = plt.subplots(figsize=(P.FIG_WIDTH_SINGLE_CM * CM, 6.5 * CM), layout="constrained")
    ax.fill_between(trav.T_C, trav.velocity_m_s - e, trav.velocity_m_s + e, color=C_BAND, lw=0)
    ax.plot(trav.T_C, trav.velocity_m_s, color=C_INK, **TRAV, label="model (travelling spike)")
    ax.plot([P.VALIDATION_T_C], [P.HH1952_VELOCITY_COMPUTED_M_S], ls="none", marker="*", ms=9,
            color=C_ACCENT, label="HH 1952 computed")
    ax.plot([P.VALIDATION_T_C], [P.HH1952_VELOCITY_MEASURED_M_S], ls="none", marker="*", ms=9,
            mfc="white", color=C_ACCENT, label="HH 1952 measured")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Conduction velocity θ (m/s)")
    ax.legend(loc="lower right")
    save(fig, "fig_velocity")


# ---------------------------------------------------------------------------
# 6. Δr
# ---------------------------------------------------------------------------
def fig_delta_r(d):
    c = d["comparison"].sort_values("T_C")
    two_eps = float(c.two_eps.iloc[0])
    fig, ax = plt.subplots(figsize=(P.FIG_WIDTH_SINGLE_CM * CM, 6.5 * CM), layout="constrained")
    T = c.T_C.to_numpy()
    ax.fill_between([T.min() - 1, T.max() + 1], -two_eps, two_eps, color=C_BAND, lw=0,
                    label=f"±2ε = ±{two_eps:.3f} (numerical uncertainty)")
    ax.axhline(0, color=C_GREY, lw=0.6)
    ax.plot(T, c.delta_r_amp4x, color=C_GREY, lw=0.9, marker="o", ms=2.5, label="membrane at 4× stimulus")
    ax.plot(T, c.delta_r, color=C_INK, marker="o", ms=3.5, label="membrane at 2× stimulus")
    ax.set_xlim(T.min() - 1, T.max() + 1)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Δr = r$_{travel}$ − r$_{membrane}$")
    ax.set_ylim(-0.1, 1.0)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.1))
    save(fig, "fig_delta_r")


# ---------------------------------------------------------------------------
# 7. Overlap fraction
# ---------------------------------------------------------------------------
def fig_overlap(d):
    sw = d["sweep"]
    mem, trav = pick(sw, "membrane", "amp2x"), pick(sw, "travelling", "cable")
    fig, ax = plt.subplots(figsize=(P.FIG_WIDTH_SINGLE_CM * CM, 6.5 * CM), layout="constrained")
    ax.plot(mem.T_C, 100 * mem.ov_frac_A, color=C_ACCENT, **MEM, label="A, membrane")
    ax.plot(trav.T_C, 100 * trav.ov_frac_A, color=C_ACCENT, **TRAV, label="A, travelling")
    ax.plot(mem.T_C, 100 * mem.ov_frac_B, color=C_INK, **MEM, label="B, membrane")
    ax.plot(trav.T_C, 100 * trav.ov_frac_B, color=C_INK, **TRAV, label="B, travelling")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Overlap fraction Q$_{ov}$/Q$_{Na}$ (%)")
    ax.set_ylim(30, 100)
    ax.legend(loc="lower left", ncols=1)
    save(fig, "fig_overlap")


# ---------------------------------------------------------------------------
# 8. Kymographs
# ---------------------------------------------------------------------------
def fig_kymograph(d):
    T_warm = d["config"]["derived"]["T_warm_C"]
    fig, axs = plt.subplots(1, 2, figsize=(P.FIG_WIDTH_FULL_CM * CM, 6.5 * CM), layout="constrained")
    for ax, T, s in zip(axs, (P.SWEEP_T_FIRST_C, T_warm), "AB"):
        c = d["traces"]["cab"][T]
        im = ax.pcolormesh(c["t_snap"], c["x"], c["V_snap"].T, cmap=V_CMAP, vmin=V_RANGE[0], vmax=V_RANGE[1],
                           shading="nearest", rasterized=True)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Position along axon x (cm)")
        width = float(pick(d["sweep"], "travelling", "cable").set_index("T_C").loc[T].width)
        t_far = c["t_snap"][np.argmax(c["V_snap"][:, -1])]
        ax.set_xlim(0, t_far + 4.0 * width)
        ax.text(0.97, 0.05, f"{T:g} °C", transform=ax.transAxes, ha="right", va="bottom")
        letter(ax, s, dx=-36)
    cb = fig.colorbar(im, ax=axs, shrink=0.9, extend="max")
    cb.set_label("V (mV)")
    save(fig, "fig_kymograph")


ALL = (fig_validation, fig_spike_anatomy, fig_temperature_currents, fig_r_vs_T, fig_velocity,
       fig_delta_r, fig_overlap, fig_kymograph)


def make_all(log=print):
    style()
    d = load()
    for f in ALL:
        f(d)
        log(f"  wrote figures/{f.__name__}.png/.pdf")
