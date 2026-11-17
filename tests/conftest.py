import pytest

from hh_energy import sweep


@pytest.fixture(scope="session")
def amplitudes():
    """Primary (2×) and variant (4×) membrane stimulus amplitudes from the 6.3 °C threshold."""
    threshold, _, _ = sweep.find_threshold()
    return sweep.stimulus_amplitudes(threshold)
