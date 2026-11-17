"""
Every constant and simulation setting used in the project (nothing is
hard-coded elsewhere).

Units used throughout the code
------------------------------
voltage V, reversal potentials E ...... mV
time t ................................ ms
current density I ..................... µA/cm²
membrane capacitance C_m .............. µF/cm²
conductance g ......................... mS/cm²
charge Q .............................. nC/cm²   (µA/cm² · ms = nC/cm²; µF/cm² · mV = nC/cm²)
length x, radius a .................... cm
axial resistivity R_i ................. Ω·cm
temperature T ......................... °C
velocity .............................. computed in cm/ms, reported in m/s (× 10)
"""
from dataclasses import dataclass, asdict, replace


# ---------------------------------------------------------------------------
# Table 1: Hodgkin–Huxley membrane (identical to NEURON's hh.mod)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MembraneParams:
    C_m: float = 1.0      # µF/cm²
    g_Na: float = 120.0   # mS/cm²  (maximal sodium conductance, ḡ_Na)
    g_K: float = 36.0     # mS/cm²  (maximal potassium conductance, ḡ_K)
    g_L: float = 0.3      # mS/cm²  (leak conductance)
    E_Na: float = 50.0    # mV
    E_K: float = -77.0    # mV
    E_L: float = -54.3    # mV

    def passive(self) -> "MembraneParams":
        """Same membrane with ḡ_Na = 0 (the passive control for spike detection)."""
        return replace(self, g_Na=0.0)


HH = MembraneParams()

# Temperature scaling of the gating rates: phi(T) = Q10 ** ((T - T_REF) / 10)
T_REF_C = 6.3         # °C, temperature at which HH fitted their rate equations
Q10 = 3.0

# vtrap switches to its Taylor form when |x/y| is below this value
VTRAP_SMALL = 1e-6

# Resting potential sanity check
V_REST_EXPECTED_MV = -65.0
V_REST_TOLERANCE_MV = 0.5
V_REST_BRACKET_MV = (-80.0, -50.0)   # bracket for the root finder

# ---------------------------------------------------------------------------
# Cable (Hodgkin & Huxley 1952 squid giant axon)
# ---------------------------------------------------------------------------
RADIUS_CM = 0.0238        # a = 238 µm
R_I_OHM_CM = 35.4         # axial resistivity R_i
LENGTH_CM = 5.0           # axon length L (sealed ends)

# Units conversion for the cable term (a / 2R_i) ∂²V/∂x²:
#   a / (2 R_i)  has units cm / (Ω·cm) = 1/Ω = S.
#   S × (mV/cm²) = mA/cm².  Multiply by 1000 to get µA/cm².
MA_PER_CM2_TO_UA_PER_CM2 = 1000.0
# Coupling coefficient D in µA/cm² per (mV/cm²); ≈ 0.336 for a = 0.0238 cm, R_i = 35.4 Ω·cm.
D_COUPLING = RADIUS_CM / (2.0 * R_I_OHM_CM) * MA_PER_CM2_TO_UA_PER_CM2

# ---------------------------------------------------------------------------
# Stimuli
# ---------------------------------------------------------------------------
STIM_START_MS = 0.1        # both spike types: pulse starts at t = 0.1 ms
STIM_DURATION_MS = 0.3     # both spike types: pulse lasts 0.3 ms

# Membrane spike: amplitude = factor × threshold(6.3 °C)
THRESHOLD_T_C = 6.3
THRESHOLD_TOL_UA = 0.1             # bisection tolerance for the threshold amplitude
THRESHOLD_BRACKET_UA = (0.0, 200.0)
AMP_FACTOR_PRIMARY = 2.0
AMP_FACTOR_VARIANT = 4.0

# Travelling spike: 2000 µA/cm² on every node with x < 0.1 cm
CABLE_STIM_AMP_UA = 2000.0
CABLE_STIM_X_MAX_CM = 0.1

# ---------------------------------------------------------------------------
# Numerical methods
# ---------------------------------------------------------------------------
DT_MEMBRANE_MS = 0.005            # RK4 step (5 µs)
REF_METHOD = "DOP853"             # scipy reference integrator
REF_RTOL = 1e-10
REF_ATOL = 1e-12

DX_CABLE_CM = 0.01                # 100 µm  -> 501 nodes on 5 cm
DT_CABLE_MS = 0.005               # 5 µs
SNAPSHOT_EVERY_MS = 0.05          # full V(x, t) is stored every 0.05 ms
RECORD_X_CM = (2.0, 2.5, 3.0, 4.0)   # V, I_Na, I_K stored every step here
MEASURE_X_CM = 2.5                # travelling spike is measured here (middle of axon)
VELOCITY_X_CM = (2.0, 3.0)        # conduction velocity between these two nodes
PROPAGATION_X_CM = 4.0            # propagated if V here exceeds PROPAGATION_V_MV
PROPAGATION_V_MV = 0.0

# Run length ("run until the spike is over"):
REST_RETURN_TOL_MV = 0.5          # stop once V is back within 0.5 mV of V_rest after the AHP
T_MAX_MS = 100.0                  # hard cap on any run
CHUNK_MS = 5.0                    # runs are extended in chunks of this length
NO_SPIKE_STOP_MS = 30.0           # a run with no spike at all is stopped here
PASSIVE_RUN_MS = 10.0             # passive control (ḡ_Na = 0) only needs its peak
CABLE_ARRIVAL_MARGIN_MV = 10.0    # "spike has arrived" at a node: V > V_rest + 10 mV
CABLE_FAR_NODE_WAIT_MS = 5.0      # after x = 2.5 cm recovers, wait at most this long for x = 4 cm to peak

# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------
SPIKE_V_MIN_MV = 0.0              # membrane spike: V_peak > 0 mV ...
SPIKE_MARGIN_OVER_PASSIVE_MV = 20.0   # ... and V_peak > passive peak + 20 mV
RAW_ONSET_MARGIN_MV = 1.0         # raw window (cable) starts when V > V_rest + 1 mV
RAW_END_TOL_MV = 1.0              # raw window ends when V is within 1 mV of V_rest after the AHP
# Diagnostics (NOT part of the spec; reported separately and labelled as such).
# The spec's 0 mV rules mark where the spike PEAK falls below 0 mV. These two
# extra rules mark where the response really disappears:
#   membrane: active response = V_peak > passive peak + 20 mV (no 0 mV requirement)
#   cable:    non-decremental wave = depolarisation at 4.0 cm ≥ 90 % of that at 2.0 cm,
#             and the depolarisation at 2.0 cm is ≥ 20 mV
DIAG_WAVE_RATIO = 0.9
DIAG_WAVE_MIN_DEPOL_MV = 20.0
DIAG_WAVE_X_CM = (2.0, 4.0)
ELEMENTARY_CHARGE_C = 1.602176634e-19
NA_PER_ATP = 3.0                  # Na+/K+ pump: 3 Na+ out per ATP
UM2_PER_CM2 = 1e8

# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------
SWEEP_T_FIRST_C = 6.3
SWEEP_T_GRID_START_C = 8.0        # then 8, 10, 12, ... in 2 °C steps
SWEEP_T_STEP_C = 2.0
SWEEP_T_EXTRA_C = (18.5,)         # always included
SWEEP_T_CAP_C = 60.0              # safety stop if spikes never fail
FAILURE_BISECT_TOL_C = 0.05
FAILURE_GRID_TOL_C = 0.01        # extra check: failure temperature re-found on each refinement grid
T_FAIL_ABOUT_SHIFT_C = 0.1       # if refinement moves T_fail by more than this, quote it as "about N °C"

VALIDATION_T_C = 18.5
HH1952_VELOCITY_COMPUTED_M_S = 18.8
HH1952_VELOCITY_MEASURED_M_S = 21.2

# Grid refinement: (dx in cm, dt in ms) for the cable; dt in ms for the membrane
CABLE_REFINEMENT = ((0.01, 0.005), (0.005, 0.0025), (0.0025, 0.00125))
MEMBRANE_REFINEMENT = (0.005, 0.0025, 0.00125)
T_WARM_OFFSET_C = 2.0             # T_warm = highest propagating coarse T − 2 °C
DELTA_R_FLAG_FACTOR = 2.0         # "real difference" if |Δr| > 2 ε

# ---------------------------------------------------------------------------
# Passive-cable test (tests/test_passive_cable.py)
# ---------------------------------------------------------------------------
PASSIVE_TEST_I_UA = 50.0          # constant current density on node x = 0
PASSIVE_TEST_RUN_MS = 60.0        # long enough to reach steady state (τ_m ≈ 1.5 ms)
PASSIVE_TEST_TOL = 0.01           # 1 %

# ---------------------------------------------------------------------------
# Figures, animation, web export
# ---------------------------------------------------------------------------
COLOR_NA = "#D9772B"
COLOR_K = "#3B6E9E"
FIG_DPI = 300
FIG_WIDTH_FULL_CM = 16.0
FIG_WIDTH_SINGLE_CM = 8.0
FONT_SIZE_PT = 9
GIF_DURATION_S = 8.0
GIF_FPS = 20
WEB_N_POSITIONS = 100
WEB_N_FRAMES = 300
JS_CONSISTENCY_T_C = (6.3, 18.5)
JS_CONSISTENCY_TOL = 0.005        # 0.5 %


def config_dict() -> dict:
    """Every parameter and setting above, for results/config.json."""
    out = {"membrane": asdict(HH), "D_coupling_uA_per_mV": D_COUPLING}
    for name, value in globals().items():
        if name.isupper() and not isinstance(value, MembraneParams):
            out[name] = list(value) if isinstance(value, tuple) else value
    return out
