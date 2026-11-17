"""Test 7: the travelling spike propagates."""
import pytest

from hh_energy import cable, metrics
from hh_energy import params as P


@pytest.mark.parametrize("T_C", [6.3, 18.5])
def test_7_cable_propagates(T_C, record_property):
    """The spike reaches x = 4.0 cm above 0 mV, and x = 2.5 cm returns to rest afterwards."""
    run = cable.simulate(T_C)
    record_property("measured", f"max V at x = 4.0 cm = {run.V_rec[P.PROPAGATION_X_CM].max():.2f} mV")
    assert metrics.cable_propagated(run.V_rec[P.PROPAGATION_X_CM])
    assert run.back_to_rest
