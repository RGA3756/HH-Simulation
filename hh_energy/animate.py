"""
GIF of the travelling spike at 6.3 °C: V(x) as a moving line, with the axon
drawn underneath as a bar coloured by voltage. Uses the V(x, t) snapshots
(every 0.05 ms) stored by sweep.run_experiments(). Units: x cm, V mV, t ms.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from . import params as P
from .figures import CM, C_INK, C_GREY, FIGURES, RESULTS, V_CMAP, V_RANGE, style
from .sweep import load_traces


def make_gif(T_C=P.SWEEP_T_FIRST_C, path=None, log=print):
    style()
    c = load_traces(RESULTS / "traces.npz")["cab"][T_C]
    x, t_snap, V_snap = c["x"], c["t_snap"], c["V_snap"]
    rest = float(V_snap[0, 0])

    # Show the spike's transit: until the far end has peaked and recovered for ~6 ms
    t_far = t_snap[np.argmax(V_snap[:, -1])]
    t_show = min(t_snap[-1], t_far + 6.0)
    n_frames = int(round(P.GIF_DURATION_S * P.GIF_FPS))
    idx = np.unique(np.round(np.linspace(0, np.searchsorted(t_snap, t_show), n_frames)).astype(int))

    fig = plt.figure(figsize=(P.FIG_WIDTH_FULL_CM * CM, 8.5 * CM), dpi=150)
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 0.6], hspace=0.55, left=0.1, right=0.97, top=0.93, bottom=0.12)
    ax = fig.add_subplot(gs[0])
    bar = fig.add_subplot(gs[1], sharex=ax)
    ax.axhline(rest, color=C_GREY, lw=0.6)
    ax.text(5.0, rest + 2, "rest", color=C_GREY, ha="right", va="bottom")
    (line,) = ax.plot(x, V_snap[0], color=C_INK, lw=1.8)
    ax.set_xlim(0, P.LENGTH_CM)
    ax.set_ylim(-90, 75)
    ax.set_ylabel("V (mV)")
    label = ax.text(0.01, 0.97, "", transform=ax.transAxes, ha="left", va="top")
    ax.text(0.99, 0.97, f"{T_C:g} °C", transform=ax.transAxes, ha="right", va="top")
    im = bar.imshow(V_snap[:1], aspect="auto", cmap=V_CMAP, vmin=V_RANGE[0], vmax=V_RANGE[1],
                    extent=(x[0] - 0.5 * (x[1] - x[0]), x[-1] + 0.5 * (x[1] - x[0]), 0, 1))
    bar.set_yticks([])
    bar.set_xlabel("Position along the axon x (cm)   (bar: axon coloured by voltage, dark = spike)")
    for s in ("left", "right", "top"):
        bar.spines[s].set_visible(False)

    def update(k):
        line.set_ydata(V_snap[k])
        im.set_data(V_snap[k:k + 1])
        label.set_text(f"t = {t_snap[k]:5.2f} ms")
        return line, im, label

    anim = FuncAnimation(fig, update, frames=idx, blit=False)
    path = path or FIGURES / "travelling_spike.gif"
    anim.save(path, writer=PillowWriter(fps=P.GIF_FPS), dpi=150)
    plt.close(fig)
    log(f"  wrote {path.relative_to(RESULTS.parent)} ({len(idx)} frames, {len(idx) / P.GIF_FPS:.1f} s)")
    return path
