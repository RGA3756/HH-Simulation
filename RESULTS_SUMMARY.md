# Results summary

*Generated automatically by `run_all.py` on 2026-11-17 21:49. Do not edit by hand: every number below is read from the files in `results/`.*

## 1. Checklist

| Check | What it verifies | Result | Measured |
|---|---|---|---|
| Test 1: Temperature factor | φ(6.3) = 1 and φ(16.3) = 3 to 1e−12 (1 case) | **PASS** | phi(6.3) = 1.000000000000000, phi(16.3) = 3.000000000000000 |
| Test 2: Removable singularities | α_m(−40) = 1, α_n(−55) = 0.1; rates finite and continuous within 1e−6 mV (4 cases) | **PASS** | alpha_m(-40) = 1.000000000000000, alpha_n(-55) = 0.100000000000000 |
| Test 3: Resting state | \|V_rest − (−65)\| < 0.5 mV and total ionic current < 1e−9 µA/cm² (1 case) | **PASS** | V_rest = -64.974052 mV, total ionic current = 5.33e-15 uA/cm2 |
| Test 4: Membrane charge balance | ∫I_stim dt = ∫I_ion dt + C_m ΔV within 0.1 % (6.3 and 25 °C) (2 cases) | **PASS** | relative imbalance = 2.92e-06 (limit 1e-3); relative imbalance = 1.55e-05 (limit 1e-3) |
| Test 5: RK4 vs DOP853 | max \|ΔV\| < 0.01 mV at 6.3 and 25 °C (2 cases) | **PASS** | max \|dV\| = 6.89e-06 mV (limit 0.01 mV); max \|dV\| = 1.35e-05 mV (limit 0.01 mV) |
| Test 6: Passive cable / units | steady state matches cosh((L − x)/λ) and the input resistance within 1 % (2 cases) | **PASS** | D = 0.336158 uA/cm2 per (mV/cm2); lambda = 0.7035 cm; max shape error = 6.0e-05; V(0) - V_rest = 0.52320 mV vs analytic 0.52321 mV (limit 1%) |
| Test 7: Propagation | travelling spike reaches x = 4 cm above 0 mV at 6.3 and 18.5 °C (2 cases) | **PASS** | max V at x = 4.0 cm = 37.93 mV; max V at x = 4.0 cm = 25.27 mV |
| Validation: HH 1952 velocity | θ at 18.5 °C vs HH's computed 18.8 m/s (PASS if within ±5 %, a threshold chosen for this project) | **PASS** | θ = 18.655 m/s, difference −0.77 % (vs measured 21.2 m/s: −12.0 %) |
| Explorer: JavaScript vs Python | r from the in-browser model within 0.5 % of Python at 6.3, 18.5 °C | **PASS** | 6.3 °C: r_JS = 13.237812, r_Py = 13.237812 (rel. diff. −8.05e−16); 18.5 °C: r_JS = 4.019942, r_Py = 4.019942 (rel. diff. 4.42e−16) |

## 2. Results by temperature

Primary definitions: Q_Na = ∫(Ĩ_Na − Ĩ_Na,rest) dt over the whole run; membrane stimulus = 2 × threshold. Travelling spike measured at x = 2.5 cm. Charges in nC/cm². Overlap A = share of Q_Na that enters after the voltage peak.

| T (°C) | φ | r membrane | r travel | Δr | flag | θ (m/s) | Q_Na mem | Q_Na trav | Q_min mem | Q_min trav | overlap A mem | overlap A trav |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6.3 | 1.000 | 13.238 | 13.959 | 0.721 | real difference | 12.28 | 1390.8 | 1436.1 | 105.1 | 102.9 | 88.5 % | 87.4 % |
| 8.0 | 1.205 | 11.074 | 11.781 | 0.707 | real difference | 13.11 | 1156.5 | 1199.6 | 104.4 | 101.8 | 86.7 % | 85.6 % |
| 10.0 | 1.502 | 9.008 | 9.700 | 0.692 | real difference | 14.11 | 932.2 | 973.3 | 103.5 | 100.3 | 84.2 % | 83.2 % |
| 12.0 | 1.871 | 7.363 | 8.040 | 0.677 | real difference | 15.15 | 752.8 | 792.5 | 102.2 | 98.6 | 81.1 % | 80.4 % |
| 14.0 | 2.330 | 6.054 | 6.718 | 0.664 | real difference | 16.22 | 609.3 | 648.0 | 100.7 | 96.5 | 77.4 % | 77.3 % |
| 16.0 | 2.903 | 5.015 | 5.667 | 0.652 | real difference | 17.30 | 494.6 | 532.6 | 98.6 | 94.0 | 73.1 % | 73.8 % |
| 18.0 | 3.616 | 4.196 | 4.834 | 0.638 | real difference | 18.38 | 402.9 | 440.2 | 96.0 | 91.1 | 68.2 % | 70.1 % |
| 18.5 | 3.820 | 4.020 | 4.654 | 0.634 | real difference | 18.65 | 383.0 | 420.0 | 95.3 | 90.3 | 66.8 % | 69.1 % |
| 20.0 | 4.505 | 3.554 | 4.174 | 0.620 | real difference | 19.46 | 329.5 | 365.8 | 92.7 | 87.6 | 62.6 % | 66.2 % |
| 22.0 | 5.612 | 3.057 | 3.652 | 0.595 | real difference | 20.49 | 270.8 | 305.7 | 88.6 | 83.7 | 56.6 % | 62.1 % |
| 24.0 | 6.990 | 2.681 | 3.241 | 0.560 | real difference | 21.44 | 223.4 | 256.5 | 83.4 | 79.1 | 50.3 % | 58.1 % |
| 26.0 | 8.708 | 2.404 | 2.919 | 0.515 | real difference | 22.26 | 184.7 | 215.6 | 76.8 | 73.9 | 44.1 % | 54.2 % |
| 28.0 | 10.848 | 2.214 | 2.671 | 0.457 | real difference | 22.83 | 152.1 | 181.2 | 68.7 | 67.8 | 38.2 % | 50.6 % |

Above that, only the membrane spike at 4 × threshold still meets the main spike rule: 30.0 °C: r = 1.706; 32.0 °C: r = 1.546.

Δr flag: "real difference" if |Δr| > 2ε with ε = sqrt(ε_mem² + ε_cab²) = 0.0203, so 2ε = 0.0407. Every variant (raw Q_Na, overlap B, 4 × stimulus) is in `results/sweep.csv` and `results/comparison.csv`.

## 3. Failure temperatures

| Spike | Rule | Main rule? | T_fail (°C) |
|---|---|---|---|
| membrane_amp2x | V_peak > 0 mV and > passive peak + 20 mV (main rule) | yes | 28.797 ± 0.016 (between 28.7812 and 28.8125) |
| membrane_amp4x | V_peak > 0 mV and > passive peak + 20 mV (main rule) | yes | 33.016 ± 0.016 (between 33.0000 and 33.0312) |
| travelling | V(x = 4.0 cm) > 0 mV (main rule) | yes | 28.859 ± 0.016 (between 28.8438 and 28.8750) |
| membrane_amp2x | DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement) | no (diagnostic) | 33.109 ± 0.016 (between 33.0938 and 33.1250) |
| membrane_amp4x | DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement) | no (diagnostic) | 39.828 ± 0.016 (between 39.8125 and 39.8438) |
| travelling | DIAGNOSTIC: depolarisation at 4.0 cm >= 90% of that at 2.0 cm (>= 20 mV) | no (diagnostic) | 33.516 ± 0.016 (between 33.5000 and 33.5312) |

**How to read this.** The main rules mark the temperature at which the spike *peak* falls below 0 mV. The heat block in this model is gradual: just above those temperatures the response is still active and the travelling wave still runs the full length of the axon with a peak below 0 mV. The diagnostic rows mark where the response really disappears. The main-rule values are therefore best described as "the spike peak falls below 0 mV at …" rather than "conduction fails at …".

## 4. Numerical uncertainty ε (grid refinement)

ε = max |value(refined) − value(base)| over the two refined grids. Computed at 6.3 °C and at T_warm = 26.0 °C (highest coarse T with a travelling spike, 28.0 °C, minus 2 °C). The larger value is used at every temperature.

| Model | Quantity | ε at 6.3 °C | ε at 26 °C | ε used |
|---|---|---|---|---|
| membrane | Q_Na | 5.13e−04 | 1.96e−04 | 5.13e−04 |
| membrane | Q_min | 7.01e−05 | 1.05e−03 | 1.05e−03 |
| membrane | r | 1.37e−05 | 3.05e−05 | 3.05e−05 |
| travelling | Q_Na | 7.14e−01 | 9.92e−02 | 7.14e−01 |
| travelling | Q_min | 4.21e−02 | 4.84e−01 | 4.84e−01 |
| travelling | r | 1.26e−02 | 2.03e−02 | 2.03e−02 |
| travelling | velocity_m_s | 1.80e−02 | 1.57e−01 | 1.57e−01 |

## 5. Soft sanity expectations (matches / differs, not pass/fail)

| Expectation | Model value | Verdict |
|---|---|---|
| V_peak at 6.3 °C roughly +30 to +50 mV | membrane 40.09 mV, travelling 37.90 mV | matches |
| Q_min at 6.3 °C about 80–115 nC/cm² | membrane 105.1, travelling 102.9 nC/cm² | matches |
| r at 6.3 °C of order 10 (Sengupta et al. 2010: about 11, i.e. about 9 % efficiency, under their own definitions) | membrane r = 13.24 (η = 7.6 %), travelling r = 13.96 (η = 7.2 %) | matches (order of magnitude; numerically 20 % / 27 % above 11) |
| r decreases from 6.3 °C to about 18 °C (Sengupta: efficiency rises to about 24 % at 18 °C) | r falls monotonically: membrane 13.24 → 4.20 (η at 18 °C = 23.8 %), travelling 13.96 → 4.83 (η = 20.7 %) | matches |
| Velocity at 18.5 °C about 16–21 m/s | 18.65 m/s | matches |

Rules used for the verdicts: V_peak and Q_min inside the stated ranges; "order 10" = between 10^0.5 and 10^1.5; "decreases" = every step on the coarse grid from 6.3 to 18 °C goes down; velocity inside 16–21 m/s.

## 6. Answers to the research question

1. **r falls steeply as the axon warms.** For the membrane spike r drops from 13.24 at 6.3 °C to 2.21 at 28.0 °C (efficiency 7.6 % → 45.2 %); for the travelling spike from 13.96 to 2.67 at 28.0 °C.
2. **The fall in r comes from Q_Na, not from Q_min.** For the membrane spike Q_Na falls 9.1-fold (1391 → 152 nC/cm²) while Q_min falls only 1.53-fold (105 → 69 nC/cm², which on its own would raise r); the share of Na⁺ entering after the peak (overlap A) drops from 89 % to 38 %.
3. **The membrane patch does not reproduce the travelling spike beyond numerical uncertainty.** r_travel − r_membrane is positive at every temperature from 6.3 to 28.0 °C: Δr = 0.721 at 6.3 °C and 0.457 at 28.0 °C, against 2ε = 0.041; the patch underestimates r by 5.4 % at 6.3 °C and 20.6 % at 28.0 °C.
4. **The size of that gap depends on the patch stimulus.** With 4 × threshold instead of 2 ×, Δr is 0.762–0.820 over the same range, so the membrane estimate is not unique.
5. **Warm limit.** The spike peak falls below 0 mV at 28.80 °C (membrane, 2 ×) and 28.86 °C (travelling, main rule), but the travelling wave keeps running without decaying until 33.52 °C (diagnostic rule).

## 7. Headline values

Every value below is computed from `results/` (sweep, refinement and the extra near-failure runs). Percentage difference = (r_propagating − r_membrane) ÷ r_membrane × 100, primary Q_Na definition, membrane stimulus 2 × threshold, base grid unless stated.

**Rounding rule.** The two differences in r are rounded so that their last digit is no finer than how much they change between the base and the finest grid at that temperature. The numerical error is given to one significant figure. The failure temperature is given as "about N °C" if refining the grid moves it by more than 0.1 °C.

| Quantity | Rounded value | Unrounded | Note |
|---|---|---|---|
| Failure temperature of the travelling spike | **about 29 °C** | 28.859 °C (base grid) | Main rule: V at x = 4.0 cm no longer exceeds 0 mV. Refining the grid moves it by up to 0.20 °C (finest grid: 29.05 °C). The wave keeps travelling, with a peak below 0 mV, until 33.52 °C (diagnostic). |
| Difference in r at 6.3 °C | **5.4%** | 5.447 % | at 6.3 °C; finest grid 5.35 % |
| Difference in r near failure | **21%** | 20.637 % | at 28.0 °C, the warmest sweep temperature where both impulses exist (r cannot be measured at T_fail itself); finest grid 19.59 % |
| Travelling r compared with membrane r | **higher** | – | r_propagating vs r_membrane at all 13 sweep temperatures with both impulses |
| Numerical error in r | **0.9%** | 0.872 % | largest change in r over every grid check (travelling, 28.0 °C, finest vs base grid) |
| Real difference? | **yes** | – | the two differences are 6.2 × and 23.7 × the numerical error; "real" means more than 2 × |
| Overlap sodium in the travelling spike | **more** (definition B), **more** (definition A) | – | overlap sodium per unit of capacitor charge (Q_ov/Q_min) in the propagating impulse; how much of the difference it explains is below |

Summary list:

```
failure temperature (travelling): about 29 °C
difference in r at 6.3 °C: 5.4%
difference in r near failure (28.0 °C): 21%
travelling r vs membrane r: higher
numerical error in r: 0.9%
overlap sodium, travelling vs membrane: more
```

### Check 1: is the difference near failure an edge effect?

Near failure both impulses get smaller, which lowers Q_min = C_m ΔV and raises r even if sodium entry does not change. Since r_prop/r_mem = (Q_Na,prop/Q_Na,mem) ÷ (Q_min,prop/Q_min,mem), the gap splits exactly into a "more sodium" part and a "smaller C_m ΔV" part (shares computed on a log scale).

| T (°C) | source | r gap | Q_Na prop vs mem | C_m ΔV prop vs mem | gap from more sodium | gap from smaller C_m ΔV |
|---|---|---|---|---|---|---|
| 6.3 | sweep | 5.45 % | 3.3 % | −2.1 % | 60 % | 40 % |
| 8.0 | sweep | 6.39 % | 3.7 % | −2.5 % | 59 % | 41 % |
| 10.0 | sweep | 7.68 % | 4.4 % | −3.0 % | 58 % | 42 % |
| 12.0 | sweep | 9.20 % | 5.3 % | −3.6 % | 58 % | 42 % |
| 14.0 | sweep | 10.97 % | 6.4 % | −4.2 % | 59 % | 41 % |
| 16.0 | sweep | 12.99 % | 7.7 % | −4.7 % | 61 % | 39 % |
| 18.0 | sweep | 15.21 % | 9.2 % | −5.2 % | 62 % | 38 % |
| 18.5 | sweep | 15.77 % | 9.7 % | −5.3 % | 63 % | 37 % |
| 20.0 | sweep | 17.45 % | 11.0 % | −5.5 % | 65 % | 35 % |
| 22.0 | sweep | 19.47 % | 12.9 % | −5.5 % | 68 % | 32 % |
| 24.0 | sweep | 20.91 % | 14.8 % | −5.1 % | 73 % | 27 % |
| 26.0 | sweep | 21.43 % | 16.8 % | −3.8 % | 80 % | 20 % |
| 28.0 | sweep | 20.64 % | 19.1 % | −1.3 % | 93 % | 7 % |
| 28.781 | extra run | 19.83 % | 20.3 % | 0.4 % | 102 % | −2 % |

- **Not an edge effect.** Near T_fail the gap comes from extra sodium (93 % at 28.0 °C, 102 % at 28.781 °C), not from a smaller C_m ΔV: close to failure both impulses shrink by about the same amount. The smaller-C_m ΔV part is largest at the cold end (40 % of the gap at 6.3 °C).
- **The trend is smooth.** The gap rises from 5.45 % at 6.3 °C to a maximum of 21.43 % at 26.0 °C (largest single step 2.21 percentage points, 16.0 → 18.0 °C), then changes by −0.79 points to 20.64 % at 28.0 °C and is 19.83 % at 28.781 °C. No jump at the end. For reference, the value at 18.5 °C is 15.77 % (not grid-checked at that temperature).

### Check 2: where was the numerical error measured, and does the failure temperature move with the grid?

- The main refinement check ran at 6.3 °C and T_warm = 26.0 °C: largest change in r 0.70 % (travelling, 26.0 °C); after the first halving only 0.47 %; membrane patch at most 1.3e−03 %.
- Repeated at 28.0 °C, next to T_fail (extra check): largest change 0.87 % (travelling). The numerical error above uses the larger of the two, 0.87 %.
- Failure temperature re-found on each grid (bisection to 0.01 °C):

| Impulse | Grid | T_fail (°C) | shift from base grid (°C) |
|---|---|---|---|
| membrane_amp2x | Δt = 5 µs | 28.785 ± 0.004 | 0.000 |
| membrane_amp2x | Δt = 2.5 µs | 28.785 ± 0.004 | 0.000 |
| membrane_amp2x | Δt = 1.25 µs | 28.785 ± 0.004 | 0.000 |
| travelling | Δx = 100 µm, Δt = 5 µs | 28.848 ± 0.004 | 0.000 |
| travelling | Δx = 50 µm, Δt = 2.5 µs | 28.988 ± 0.004 | 0.141 |
| travelling | Δx = 25 µm, Δt = 1.25 µs | 29.051 ± 0.004 | 0.203 |

- The propagating impulse's failure temperature moves by up to 0.20 °C (more than 0.1 °C), so it is best given as **about 29 °C** rather than 28.9 °C. The successive shifts (0.141, 0.062 °C) get smaller, so the value is settling as the grid is refined.
- The membrane impulse's failure temperature does not move at this resolution (largest shift 0.000 °C, bisection tolerance 0.01 °C).

### Does overlap explain the difference?

Since r = Q_ov/Q_min + (Q_Na − Q_ov)/Q_min, the difference in r splits exactly into an overlap part and the rest.

| Overlap definition | T (°C) | Q_ov membrane | Q_ov propagating | share of Q_Na (membrane / propagating) | Q_ov/Q_min (membrane / propagating) | part of the r difference explained |
|---|---|---|---|---|---|---|
| B: Na⁺ entering while K⁺ flows out at the same time | 6.3 | 1295.7 nC/cm² | 1320.9 nC/cm² | 93.2 % / 92.0 % | 12.33 / 12.84 | 70 % |
| B: Na⁺ entering while K⁺ flows out at the same time | 28.0 | 95.5 nC/cm² | 103.9 nC/cm² | 62.8 % / 57.4 % | 1.39 / 1.53 | 31 % |
| A: Na⁺ entering after the voltage peak | 6.3 | 1231.4 nC/cm² | 1254.9 nC/cm² | 88.5 % / 87.4 % | 11.72 / 12.20 | 66 % |
| A: Na⁺ entering after the voltage peak | 28.0 | 58.1 nC/cm² | 91.6 nC/cm² | 38.2 % / 50.6 % | 0.85 / 1.35 | 110 % |

- Definition B: Na⁺ entering while K⁺ flows out at the same time: the propagating impulse has **more** overlap sodium per unit of capacitor charge; it explains 70 % of the r difference at 6.3 °C and 31 % at 28.0 °C; overlap explains most of the difference at one end but only 31 % at the other, so it is only a partial explanation. As a share of Q_Na, the propagating impulse's overlap is smaller at every temperature, so the difference is in the amount of overlap sodium, not in the overlap fraction.
- Definition A: Na⁺ entering after the voltage peak: the propagating impulse has **more** overlap sodium per unit of capacitor charge; it explains 66 % of the r difference at 6.3 °C and 110 % at 28.0 °C; overlap explains most of the difference at both ends. As a share of Q_Na, the propagating impulse's overlap is bigger only at some temperatures, so the difference is in the amount of overlap sodium, not in the overlap fraction.
