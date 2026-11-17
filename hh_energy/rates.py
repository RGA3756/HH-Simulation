"""
Hodgkin–Huxley gating kinetics (NEURON hh.mod form).

V in mV, rates α and β in 1/ms, time constants τ in ms, T in °C.
All functions accept floats or numpy arrays.
"""
import numpy as np

from . import params as P


def phi(T_C):
    """Temperature factor φ(T) = 3^((T − 6.3)/10) applied to every gating rate."""
    return P.Q10 ** ((np.asarray(T_C, dtype=float) - P.T_REF_C) / 10.0)


def vtrap(x, y):
    """x / (exp(x/y) − 1), with the Taylor form y·(1 − x/(2y)) when |x/y| < 1e−6.

    This removes the 0/0 at x = 0 (V = −40 mV for α_m, V = −55 mV for α_n).
    """
    x = np.asarray(x, dtype=float)
    u = x / y
    small = np.abs(u) < P.VTRAP_SMALL
    with np.errstate(divide="ignore", invalid="ignore"):
        regular = x / (np.exp(u) - 1.0)
    out = np.where(small, y * (1.0 - u / 2.0), regular)
    return out[()] if out.ndim == 0 else out


def alpha_m(V):
    return 0.1 * vtrap(-(V + 40.0), 10.0)


def beta_m(V):
    return 4.0 * np.exp(-(V + 65.0) / 18.0)


def alpha_h(V):
    return 0.07 * np.exp(-(V + 65.0) / 20.0)


def beta_h(V):
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))


def alpha_n(V):
    return 0.01 * vtrap(-(V + 55.0), 10.0)


def beta_n(V):
    return 0.125 * np.exp(-(V + 65.0) / 80.0)


def all_rates(V):
    """(α_m, β_m, α_h, β_h, α_n, β_n) at 6.3 °C, in 1/ms."""
    return alpha_m(V), beta_m(V), alpha_h(V), beta_h(V), alpha_n(V), beta_n(V)


def steady_state(V):
    """(m∞, h∞, n∞) = α/(α+β). Independent of temperature."""
    am, bm, ah, bh, an, bn = all_rates(V)
    return am / (am + bm), ah / (ah + bh), an / (an + bn)


def time_constants(V, T_C):
    """(τ_m, τ_h, τ_n) = 1/(φ(T)(α+β)) in ms."""
    am, bm, ah, bh, an, bn = all_rates(V)
    f = phi(T_C)
    return 1.0 / (f * (am + bm)), 1.0 / (f * (ah + bh)), 1.0 / (f * (an + bn))
