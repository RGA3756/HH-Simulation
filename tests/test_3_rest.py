"""Test 3: resting state."""
from hh_energy import params as P
from hh_energy.membrane import resting_state, total_steady_state_current


def test_3_resting_state(record_property):
    """V_rest is within 0.5 mV of −65 mV and the total ionic current there is < 1e−9 µA/cm²."""
    rest = resting_state()
    record_property("measured", f"V_rest = {rest.V:.6f} mV, total ionic current = {total_steady_state_current(rest.V):.2e} uA/cm2")
    assert abs(rest.V - (-65.0)) < P.V_REST_TOLERANCE_MV
    assert abs(total_steady_state_current(rest.V)) < 1e-9
    assert abs(rest.I_Na + rest.I_K + rest.I_L) < 1e-9
