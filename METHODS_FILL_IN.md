# Values for the Methods section

*Generated automatically by `run_all.py` on 2026-10-22 18:27. Every value is read from `results/`.*

## Computer and software

- Computer / machine model: Mac15,3; CPU: Apple M3 (arm64)
- Operating system: macOS-26.3.1-arm64-arm-64bit-Mach-O
- Python 3.14.7; numpy 2.5.4; scipy 1.18.1; matplotlib 3.11.2; pandas 3.0.6
- Node.js (explorer consistency check only): v24.21.0

## Model and numerics

- Membrane: C_m = 1 µF/cm², ḡ_Na = 120, ḡ_K = 36, g_L = 0.3 mS/cm², E_Na = 50, E_K = -77, E_L = -54.3 mV (NEURON hh.mod rate functions).
- Temperature factor φ(T) = 3^((T − 6.3)/10), applied to the gating rates only.
- Resting potential V_rest = −64.9741 mV (root of the total steady-state ionic current).
- Membrane patch: classical RK4, Δt = 5 µs, stimulus integrated piece by piece; reference solver scipy DOP853 (rtol 1e-10, atol 1e-12).
- Cable: a = 238 µm, R_i = 35.4 Ω·cm, L = 5 cm, sealed ends, vertex-centred grid with Δx = 100 µm (501 nodes), Δt = 5 µs; gates by exponential Euler, voltage by backward Euler (tridiagonal solve). Coupling coefficient a/(2R_i) = 0.3362 µA/cm² per (mV/cm²). Explicit stability number DΔt/(C_mΔx²) = 16.8 (explicit stepping would need ≤ 0.5).
- Travelling-spike stimulus: 2000 µA/cm² for 0.3 ms from t = 0.1 ms on all nodes with x < 0.1 cm; measured at x = 2.5 cm; velocity from the 0 mV crossings at x = 2.0 and 3.0 cm.
- Run length: until V (membrane, or the cable node at 2.5 cm) is back within 0.5 mV of V_rest after the after-hyperpolarisation (cap 100 ms). Charges integrated with the trapezoid rule.

## Stimulus amplitude (membrane spike)

- Threshold at 6.3 °C (0.3 ms pulse, bisection to 0.1 µA/cm²): **21.875 µA/cm²** (bracket 21.777–21.875)
- Primary amplitude (2 × threshold): **43.750 µA/cm²**; variant (4 ×): **87.500 µA/cm²**; same amplitude at every temperature.

## RK4 vs DOP853

- 6.3 °C: max |V_RK4 − V_DOP853| = **6.89e−06 mV** (during the spike: 6.89e−06 mV); relative difference in Q_Na = **2.94e−08**
- 25.0 °C: max |V_RK4 − V_DOP853| = **1.35e−05 mV** (during the spike: 1.35e−05 mV); relative difference in Q_Na = **−7.19e−07**

## Numerical uncertainty ε (grid refinement)

- Cable grids (Δx, Δt): (100 µm, 5 µs), (50 µm, 2.5 µs), (25 µm, 1.25 µs); membrane Δt: 5 µs, 2.5 µs, 1.25 µs.
- Refinement temperatures: 6.3 °C and T_warm = 26.0 °C.
- ε for r (membrane): 1.37e−05 at 6.3 °C, 3.05e−05 at 26.0 °C → used: **3.05e−05**
- ε for r (travelling): 1.26e−02 at 6.3 °C, 2.03e−02 at 26.0 °C → used: **2.03e−02**
- ε for θ (travelling): 1.80e−02 m/s at 6.3 °C, 1.57e−01 m/s at 26.0 °C → used: **1.57e−01 m/s**
- ε for Q_Na (membrane): 5.13e−04 nC/cm² at 6.3 °C, 1.96e−04 nC/cm² at 26.0 °C → used: **5.13e−04 nC/cm²**
- ε for Q_Na (travelling): 7.14e−01 nC/cm² at 6.3 °C, 9.92e−02 nC/cm² at 26.0 °C → used: **7.14e−01 nC/cm²**
- ε for Q_min (membrane): 7.01e−05 nC/cm² at 6.3 °C, 1.05e−03 nC/cm² at 26.0 °C → used: **1.05e−03 nC/cm²**
- ε for Q_min (travelling): 4.21e−02 nC/cm² at 6.3 °C, 4.84e−01 nC/cm² at 26.0 °C → used: **4.84e−01 nC/cm²**
- Extra check, not part of the planned method: the same refinement at 28.0 °C (warmest sweep temperature with both impulses) changes r by at most 0.87 % (travelling).
- Combined ε for Δr = sqrt(ε_mem² + ε_cab²) = **0.0203**; flag threshold 2ε = 0.0407.

## Temperature grid

- Coarse sweep: 6.3, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 18.5, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0 °C (6.3 °C, then 2 °C steps from 8 °C, plus 18.5 °C, until every spike type had failed).

## Failure temperatures (bisection to 0.05 °C)

- membrane_amp2x: **28.797 ± 0.016 °C** (V_peak > 0 mV and > passive peak + 20 mV (spec))
- membrane_amp4x: **33.016 ± 0.016 °C** (V_peak > 0 mV and > passive peak + 20 mV (spec))
- travelling: **28.859 ± 0.016 °C** (V(x = 4.0 cm) > 0 mV (spec))
- membrane_amp2x: **33.109 ± 0.016 °C** (DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement)) — DIAGNOSTIC, not the spec's rule
- membrane_amp4x: **39.828 ± 0.016 °C** (DIAGNOSTIC: V_peak > passive peak + 20 mV (no 0 mV requirement)) — DIAGNOSTIC, not the spec's rule
- travelling: **33.516 ± 0.016 °C** (DIAGNOSTIC: depolarisation at 4.0 cm >= 90% of that at 2.0 cm (>= 20 mV)) — DIAGNOSTIC, not the spec's rule
- Extra check, not part of the planned method: travelling failure temperature re-found on each refinement grid (bisection to 0.01 °C): 28.848 °C, 28.988 °C, 29.051 °C (largest shift from the base grid 0.203 °C)
- Extra check, not part of the planned method: membrane_amp2x failure temperature re-found on each refinement grid (bisection to 0.01 °C): 28.785 °C, 28.785 °C, 28.785 °C (largest shift from the base grid 0.000 °C)

## Validation against Hodgkin & Huxley (1952)

- dx = 100 um, dt = 5 us: θ(18.5 °C) = **18.655 m/s**; vs HH computed 18.8 m/s: −0.77 %; vs HH measured 21.2 m/s: −12.00 % (base grid)
- dx = 50 um, dt = 2.5 us: θ(18.5 °C) = **18.692 m/s**; vs HH computed 18.8 m/s: −0.57 %; vs HH measured 21.2 m/s: −11.83 %
- dx = 25 um, dt = 1.25 us: θ(18.5 °C) = **18.710 m/s**; vs HH computed 18.8 m/s: −0.48 %; vs HH measured 21.2 m/s: −11.74 %

## Explorer consistency check

- 6.3 °C: r (JavaScript) = 13.237812, r (Python) = 13.237812, relative difference −8.05e−16 (tolerance 0.5 %)
- 18.5 °C: r (JavaScript) = 4.019942, r (Python) = 4.019942, relative difference 4.42e−16 (tolerance 0.5 %)
