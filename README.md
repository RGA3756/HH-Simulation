# How much energy does a nerve spike cost at different temperatures?

This project simulates the Hodgkin–Huxley (HH) model of the squid giant axon and measures how much sodium
(Na⁺) one spike lets into the cell, compared with the smallest amount that could possibly do the job.

- **Q_min = C_m (V_peak − V_rest)** is the charge needed to "fill the membrane capacitor" from rest to the
  spike peak. It is the theoretical minimum.
- **Q_Na** is the sodium charge that actually enters. It is larger, because sodium keeps flowing in while
  potassium (K⁺) is already flowing out (the "overlap").
- **r = Q_Na / Q_min** says how many times the minimum the spike really uses. **Efficiency η = 100 % / r.**
  Every 3 Na⁺ ions that enter cost one ATP molecule when the Na⁺/K⁺ pump throws them out again.

Everything is computed for two kinds of spike, at temperatures from 6.3 °C upward:

1. **Membrane spike:** one small patch of membrane ("space-clamped"), started by a short current pulse.
2. **Travelling spike:** a spike running along a 5 cm axon, measured at the middle (x = 2.5 cm).

Every number in the paper outputs comes from code in this folder that was actually run. Nothing is
estimated or copied in by hand.

---

## Install and run

You need Python 3.11 or newer, and Node.js for the explorer check. In a terminal, inside this folder:

```bash
python3 -m venv .venv
```

```bash
.venv/bin/python -m pip install -r requirements.txt
```

```bash
.venv/bin/python run_all.py
```

`run_all.py` deletes the old outputs and rebuilds everything from scratch (about 2.5 minutes on an Apple M3 laptop).
If any test fails, it stops before running the experiments. A copy of everything it prints is saved in
`results/run_log.txt`.

To run only the tests:

```bash
.venv/bin/python -m pytest -v
```

To use the interactive explorer, open `web/explorer.html` in any web browser (it works offline).

---

## What each file is

| File | What it contains |
|---|---|
| `hh_energy/params.py` | Every constant (Table 1) and every simulation setting. Nothing is hard-coded anywhere else. |
| `hh_energy/rates.py` | The gate equations: opening/closing rates α and β, steady states, time constants, and the temperature factor φ(T). |
| `hh_energy/membrane.py` | The membrane-patch model, the RK4 integrator, and the DOP853 reference solver. |
| `hh_energy/cable.py` | The 5 cm axon model (exponential Euler for the gates, backward Euler for the voltage). |
| `hh_energy/metrics.py` | All measurements: Q_Na (two definitions), Q_min, r, efficiency, overlap (two variants), velocity, spike detection, ATP. |
| `hh_energy/sweep.py` | The experiments: threshold, temperature sweep, failure temperatures, validation, grid refinement, comparison. |
| `hh_energy/figures.py` | All paper figures. |
| `hh_energy/animate.py` | The GIF of the travelling spike. |
| `hh_energy/export_web.py` | Builds `web/explorer.html` and runs the JavaScript-vs-Python check. |
| `hh_energy/explorer_template.html` | The explorer page before the data are inserted. |
| `hh_energy/report.py` | Writes `RESULTS_SUMMARY.md`, `METHODS_FILL_IN.md` and `results/env_info.txt`. |
| `tests/` | The automatic tests (see below), plus `js_consistency.mjs`, which runs the explorer's model headless in Node.js. |
| `run_all.py` | One command that reproduces everything. |
| `results/sweep.csv` | One row per temperature × spike type × variant, with every measurement. |
| `results/comparison.csv` | r_membrane, r_travel, Δr and the "real difference" flag at each temperature, for each variant. |
| `results/refinement.csv`, `results/epsilon.csv` | Grid-refinement values and the numerical uncertainty ε. |
| `results/near_failure.csv`, `results/refinement_edge.csv`, `results/failure_vs_grid.csv` | Extra checks close to the failure temperature (not part of the planned Methods): r, Q_Na and Q_min of both impulse types at the warmest temperatures where both exist, the grid refinement repeated there, and the failure temperature re-found on every refinement grid. Used in section 7 of `RESULTS_SUMMARY.md`. |
| `results/failure_temps.csv` | Failure temperatures (the spec's rules, plus clearly labelled diagnostic rules). |
| `results/validation.csv` | Velocity against Hodgkin & Huxley (1952), and RK4 against DOP853. |
| `results/js_consistency.csv` | The explorer's JavaScript r against the Python r. |
| `results/config.json` | Every parameter and setting used, plus derived values (threshold, T_warm, ...). |
| `results/env_info.txt` | Computer, operating system, Python and package versions. |
| `results/traces.npz` | The simulated voltage and current traces that the figures are drawn from. |
| `RESULTS_SUMMARY.md` | Auto-generated: checklist, results table, soft expectations, answers to the research question, and the values for the paper sentence ([T_fail], [A], [B], [E], higher/lower, more/less). |
| `METHODS_FILL_IN.md` | Auto-generated: the exact values to paste into the Methods section. |

---

## The figures (in `figures/`, as 300 dpi PNG and as PDF)

- **fig_validation**: (A) the membrane spike at 6.3 °C from our RK4 solver and from a very precise
  reference solver, drawn on top of each other. (B) The difference between them on a log scale, far below the
  0.01 mV test limit. This shows the numerical method is not distorting the spike.
- **fig_spike_anatomy** (the key teaching figure), membrane spike at 6.3 °C. (A) Voltage, with Q_min marked as
  the rise from rest to peak. (B) The three gates: m and h belong to the sodium channel, n to the potassium
  channel. (C) Sodium entering (orange) and potassium leaving (blue). Dark orange is the first Q_min worth of
  sodium, which is enough to fill the capacitor. Grey is the overlap, where both flow at once and largely
  cancel. Note the dip in the sodium current at the peak: there V is close to E_Na, so there is almost no
  driving force.
- **fig_temperature_currents**: the same two currents at 6.3 °C, a middle temperature and T_warm, for both
  spike types. Time is divided by the spike width so the shapes can be compared. The grey overlap shrinks as it
  gets warmer.
- **fig_r_vs_T**: (A) r against temperature for both spike types, with the ±ε numerical-uncertainty band (so
  thin it is hard to see) and vertical lines where the spike peak drops below 0 mV. (B) Q_Na and Q_min against
  temperature on a log axis. Q_Na falls a lot; Q_min only a little.
- **fig_velocity**: conduction velocity of the travelling spike against temperature, with Hodgkin & Huxley's
  1952 computed (filled star) and measured (open star) values at 18.5 °C.
- **fig_delta_r**: Δr = r_travel − r_membrane. The grey band is ±2ε. Points outside it are real differences,
  not numerical noise.
- **fig_overlap**: the fraction of Q_Na that is "overlap", by both definitions (A: sodium entering after the
  peak; B: sodium entering while potassium also flows out), for both spike types.
- **fig_kymograph**: the travelling spike as a heat map, with position along the axon on the vertical axis and
  time on the horizontal axis. The spike is the dark tilted stripe, and its slope is the velocity. Left 6.3 °C,
  right T_warm.
- **travelling_spike.gif**: an 8-second animation of the voltage along the axon at 6.3 °C, with the axon
  drawn underneath and coloured by voltage.

Colours: Na⁺ orange (#D9772B), K⁺ blue (#3B6E9E), matching Figure 1 of the paper. Membrane spike = solid
line and filled markers; travelling spike = dashed line and open markers.

---

## The tests, and why each one matters

| Test | What it checks | Why it matters |
|---|---|---|
| 1 | φ(6.3) = 1 and φ(16.3) = 3 | Temperature enters the model only through φ. If this is wrong, every temperature result is wrong. |
| 2 | α_m and α_n give the right values at the "0/0" voltages (−40 and −55 mV) and change smoothly near them | A naive formula divides 0 by 0 there and returns "not a number", which would silently break a simulation. |
| 3 | The resting potential is −65 mV (±0.5) and no net current flows there | Every simulation starts from this state. A wrong start would add a fake drift to every charge. |
| 4 | Charge balance: the stimulus charge equals the ionic charge plus the change of capacitor charge | Bookkeeping check that the currents, signs and units add up (no charge appears or disappears). |
| 5 | RK4 agrees with an ultra-precise solver to within 0.01 mV at 6.3 and 25 °C | Shows the time step is small enough, including at warm temperatures where everything is faster. |
| 6 | A passive cable with a constant injected current matches the textbook cosh solution and input resistance within 1 % | Checks the ×1000 units conversion in the cable equation, using λ computed independently in ohms and centimetres. A wrong factor would make λ about 30 times too small. |
| 7 | The travelling spike reaches x = 4 cm at 6.3 and 18.5 °C | The cable model actually conducts spikes. |

`RESULTS_SUMMARY.md` lists every test with its result and the actual measured margin.

---

## Things to know when you write the paper

- **"Failure temperature" means "the spike peak drops below 0 mV".** In this model the heat block is
  gradual. Just above the spec's failure temperatures the response is still active, and the travelling wave
  still runs the whole axon with a peak below 0 mV. `results/failure_temps.csv` therefore also lists, clearly
  labelled as *diagnostic* (not part of the original plan), the temperatures at which the response really
  disappears.
- **The membrane stimulus is fixed** at 2 × (variant: 4 ×) the 6.3 °C threshold at every temperature, as
  planned. The membrane spike's failure temperature, and the size of Δr, therefore partly depend on that
  choice. The 4 × variant shows how much.
- **Near the stimulated end of the axon** (x < 0.1 cm) the voltage briefly overshoots to about +70 mV because
  the 2000 µA/cm² pulse is very strong. This is far from the measuring point at 2.5 cm and does not affect the
  measurements there.
- **The HH 1952 velocity check** is reported as PASS if within ±5 %. That 5 % is my choice; the plan only asked
  for the % difference.
- The ε bands in the figures are real but very thin, because the numerical error is small compared with the
  effects being measured.

## Units used everywhere

V in mV, t in ms, current density in µA/cm², C_m in µF/cm², conductance in mS/cm², charge in nC/cm²
(µA/cm² × ms = nC/cm²), length in cm, R_i in Ω·cm, T in °C. Velocity is computed in cm/ms and reported in m/s
(× 10).
