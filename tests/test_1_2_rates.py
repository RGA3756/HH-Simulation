"""Tests 1 and 2: temperature factor and the removable singularities of α_m, α_n."""
import numpy as np
import pytest

from hh_energy import params as P
from hh_energy.rates import alpha_m, alpha_n, phi, vtrap


def test_1_phi(record_property):
    """φ(6.3) = 1 and φ(16.3) = 3: the Q10 = 3 rule is applied correctly."""
    record_property("measured", f"phi(6.3) = {phi(6.3):.15f}, phi(16.3) = {phi(16.3):.15f}")
    assert abs(phi(6.3) - 1.0) < 1e-12
    assert abs(phi(16.3) - 3.0) < 1e-12


def test_2_vtrap_limits(record_property):
    """At the singular voltages the rates take their analytic limits."""
    record_property("measured", f"alpha_m(-40) = {alpha_m(-40.0):.15f}, alpha_n(-55) = {alpha_n(-55.0):.15f}")
    assert abs(alpha_m(-40.0) - 1.0) < 1e-12
    assert abs(alpha_n(-55.0) - 0.1) < 1e-12


def _series(x, y):
    """x / (exp(x/y) − 1) from its Taylor series (accurate for |x/y| ≤ 1e−3)."""
    u = x / y
    return y * (1 - u / 2 + u**2 / 12 - u**4 / 720)


@pytest.mark.parametrize("rate, V_s, scale", [(alpha_m, -40.0, 0.1), (alpha_n, -55.0, 0.01)])
def test_2_rates_continuous(rate, V_s, scale):
    """Within 1e−6 mV of each singular voltage the rate is finite and matches the exact series."""
    for d in np.concatenate([-np.logspace(-15, -6, 40), [0.0], np.logspace(-15, -6, 40)]):
        V = V_s + d
        val = rate(V)
        assert np.isfinite(val)
        assert abs(val - scale * _series(-(V - V_s), 10.0)) < 1e-9


def test_2_vtrap_switch_is_seamless():
    """The Taylor branch and the exact formula agree at the |x/y| = 1e−6 switch point."""
    y = 10.0
    for x in (P.VTRAP_SMALL * y * (1 - 1e-9), P.VTRAP_SMALL * y * (1 + 1e-9)):
        for s in (1, -1):
            assert abs(vtrap(s * x, y) - _series(s * x, y)) < 1e-9
