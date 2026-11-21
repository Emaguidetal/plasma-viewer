"""Lumped CR radiation and a conditional optically-thin thermal continuum diagnostic.

CR coefficients and their provenance are loaded by GasModel. The free-free
formula is Rybicki & Lightman (1979), Radiative Processes in Astrophysics,
Eq. 5.15b, converted to SI. It assumes Maxwellian electrons, optically thin
emission and a constant bolometric Gaunt factor. It is not a helium line
spectrum or a calibrated camera signal and is NOT an additional kinetic sink.
"""
import numpy as np

from . import plasma3d as k3
from . import plasma3d_cr as kc

H = 6.62607015e-34
C = 299792458.0


def cr_emission(state, prm):
    if state.shape[0] != 6 or not np.isfinite(state).all() or np.any(state < 0):
        raise ValueError("emission requires six finite nonnegative plasma state fields")
    resonance = state[3] * prm[kc.P_ARES] * prm[kc.P_FRES] * prm[kc.P_ESTAR] * k3.QE
    excimer = state[4] * prm[kc.P_AX] * prm[kc.P_EX] * k3.QE
    return resonance, excimer


def free_free_estimate(state, wavelength_band_m=(400e-9, 700e-9), gaunt=1.2):
    lo, hi = wavelength_band_m
    if not 0 < lo < hi or not np.isfinite(gaunt) or gaunt <= 0:
        raise ValueError("invalid wavelength band or Gaunt factor")
    if state.shape[0] != 6 or not np.isfinite(state).all() or np.any(state < 0):
        raise ValueError("continuum requires six finite nonnegative state fields")
    ne = state[0] + state[1] + 2 * state[5]
    mean_energy_ev = np.divide(state[2], ne, out=np.zeros_like(ne), where=ne > 0)
    te_k = mean_energy_ev * (k3.QE / (1.5 * k3.KB))
    if not np.isfinite(te_k).all():
        raise ArithmeticError("non-finite electron temperature in continuum diagnostic")
    charge_squared_ions = state[0] + state[1] + 4 * state[5]
    bolometric = 1.426e-40 * gaunt * ne * charge_squared_ions * np.sqrt(te_k)
    thermal = k3.KB * np.maximum(te_k, 1e-100)
    fraction = np.exp(-H * C / (hi * thermal)) - np.exp(-H * C / (lo * thermal))
    return bolometric, bolometric * fraction
