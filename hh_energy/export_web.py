"""
Builds web/explorer.html (single self-contained file, works offline) and runs
the JavaScript-vs-Python consistency check.

The page embeds a JSON block with the model constants, the sweep results and
the travelling-spike voltage maps V(x, t) (≈100 positions × ≤300 frames per
temperature, quantised to one byte for display). Units: mV, ms, cm, nC/cm².
"""
import base64
import datetime
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from . import params as P
from .figures import RESULTS, V_RANGE
from .sweep import load_traces

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = Path(__file__).with_name("explorer_template.html")
OUT = ROOT / "web" / "explorer.html"
JS_CHECK = ROOT / "tests" / "js_consistency.mjs"

MODEL_PARAMS = ("Q10", "T_REF_C", "VTRAP_SMALL", "V_REST_BRACKET_MV", "DT_MEMBRANE_MS", "STIM_START_MS",
                "STIM_DURATION_MS", "CHUNK_MS", "T_MAX_MS", "NO_SPIKE_STOP_MS", "PASSIVE_RUN_MS",
                "REST_RETURN_TOL_MV", "SPIKE_V_MIN_MV", "SPIKE_MARGIN_OVER_PASSIVE_MV",
                "ELEMENTARY_CHARGE_C", "NA_PER_ATP", "UM2_PER_CM2")
V_OFFSET, V_SPAN = -100.0, 180.0      # 1-byte quantisation range for the voltage maps (mV)


def build_data():
    sweep = pd.read_csv(RESULTS / "sweep.csv")
    fail = pd.read_csv(RESULTS / "failure_temps.csv")
    cfg = json.loads((RESULTS / "config.json").read_text())
    traces = load_traces(RESULTS / "traces.npz")
    ok = sweep[sweep.spiked.astype(bool) & (sweep.qna_def == "primary")]
    mem = ok[(ok.spike_type == "membrane") & (ok.stim_variant == "amp2x")].sort_values("T_C")
    trav = ok[ok.spike_type == "travelling"].sort_values("T_C")
    spec = fail[fail.is_spec_definition.astype(bool)]

    params = {k: (list(getattr(P, k)) if isinstance(getattr(P, k), tuple) else getattr(P, k)) for k in MODEL_PARAMS}
    params.update({k: float(v) for k, v in cfg["membrane"].items()})

    temps = []
    x_all = None
    for _, row in trav.iterrows():
        c = traces["cab"][row.T_C]
        x_all = c["x"]
        stride_x = max(1, int(round((x_all.size - 1) / P.WEB_N_POSITIONS)))
        V = c["V_snap"][:, ::stride_x]
        t_far = c["t_snap"][np.argmax(c["V_snap"][:, -1])]
        n_show = int(np.searchsorted(c["t_snap"], t_far + 4.0 * row.width)) + 1
        stride_t = int(np.ceil(n_show / P.WEB_N_FRAMES))
        V = V[:n_show:stride_t]
        q = np.clip(np.round((V - V_OFFSET) / V_SPAN * 255.0), 0, 255).astype(np.uint8)
        temps.append({"T": float(row.T_C), "dt_frame": float(c["t_snap"][1] - c["t_snap"][0]) * stride_t,
                      "n_frames": int(V.shape[0]), "velocity": float(row.velocity_m_s), "r": float(row.r),
                      "data": base64.b64encode(q.tobytes()).decode("ascii")})
    stride_x = max(1, int(round((x_all.size - 1) / P.WEB_N_POSITIONS)))
    return {
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "params": params,
        "threshold": cfg["derived"]["threshold_uA_cm2"],
        "amp_primary": cfg["derived"]["amp_primary_uA_cm2"],
        "T_min": P.SWEEP_T_FIRST_C,
        "T_max": float(spec[spec.spike == "travelling"].T_fail_C.iloc[0]),
        "failure_spec": [{"spike": s, "T": float(T)} for s, T in zip(spec.spike, spec.T_fail_C)],
        "v_range": list(V_RANGE),
        "sweep": {
            "membrane": [{"T": float(T), "r": float(r), "Q_Na": float(q), "Q_min": float(qm)}
                         for T, r, q, qm in zip(mem.T_C, mem.r, mem.Q_Na, mem.Q_min)],
            "travelling": [{"T": float(T), "r": float(r)} for T, r in zip(trav.T_C, trav.r)],
        },
        "cable": {"x": [round(float(v), 6) for v in x_all[::stride_x]], "v_offset": V_OFFSET,
                  "v_scale": V_SPAN / 255.0, "rest": cfg["derived"]["rest"]["V"], "temps": temps},
    }


def build_html(log=print):
    data = build_data()
    html = TEMPLATE.read_text().replace("__HH_DATA__", json.dumps(data, separators=(",", ":")))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(html)
    n = sum(t["n_frames"] for t in data["cable"]["temps"])
    log(f"  wrote web/explorer.html ({OUT.stat().st_size / 1e6:.2f} MB; {len(data['cable']['temps'])} temperatures, "
        f"{len(data['cable']['x'])} positions, {n} frames in total)")
    return data


def js_consistency(log=print):
    """Run the page's model code in Node at 6.3 and 18.5 °C and compare r with the Python sweep."""
    node = shutil.which("node")
    if node is None:
        raise RuntimeError("Node.js not found: the JS consistency check cannot run")
    temps = [str(T) for T in P.JS_CONSISTENCY_T_C]
    out = subprocess.run([node, str(JS_CHECK), str(OUT), *temps], capture_output=True, text=True, check=True)
    js = {round(d["T"], 6): d for d in json.loads(out.stdout)}
    sweep = pd.read_csv(RESULTS / "sweep.csv")
    py = sweep[(sweep.spike_type == "membrane") & (sweep.stim_variant == "amp2x") & (sweep.qna_def == "primary")]
    py = py.set_index("T_C")
    rows = []
    for T in P.JS_CONSISTENCY_T_C:
        j, p_ = js[round(T, 6)], py.loc[T]
        rows.append({"T_C": T, "r_python": p_.r, "r_js": j["r"], "rel_diff_r": (j["r"] - p_.r) / p_.r,
                     "Q_Na_python": p_.Q_Na, "Q_Na_js": j["Q_Na"], "Q_min_python": p_.Q_min, "Q_min_js": j["Q_min"],
                     "spiked_js": j["spiked"], "tolerance": P.JS_CONSISTENCY_TOL,
                     "pass": bool(abs(j["r"] - p_.r) / p_.r < P.JS_CONSISTENCY_TOL and j["spiked"]),
                     "node_version": subprocess.run([node, "--version"], capture_output=True, text=True).stdout.strip()})
        log(f"  JS vs Python at {T} °C: r_js = {j['r']:.6f}, r_python = {p_.r:.6f}, "
            f"relative difference = {(j['r'] - p_.r) / p_.r:.2e}")
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "js_consistency.csv", index=False)
    return df
