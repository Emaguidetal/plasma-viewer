"""Conservative 1D Cartesian/cylindrical Euler gas dynamics with passive plasma states.

Cylindrical volumes are annular areas per unit laser-axis length. When six
plasma tracers are present, total energy includes electron and chemical energy.
The gas and electrons have gamma=5/3; chemical energy does not exert pressure.
This is an inviscid, local model, not a full cavity gas-flow calculation.
"""
from __future__ import annotations

import math

import numpy as np

from . import plasma3d as k3
from . import plasma3d_cr as kc

M_HE = 4.002602 * 1.66053906660e-27
GAMMA = 5.0 / 3.0


def chemical_energy(state, prm):
    ip = prm[k3.P_IP]
    return k3.QE * (ip * state[0] + (ip - prm[kc.P_DION]) * state[1]
                    + prm[kc.P_ESTAR] * state[3] + prm[kc.P_EX] * state[4]
                    + (ip + prm[kc.P_I2]) * state[5])


class Euler1D:
    def __init__(self, edges, *, cylindrical=False, gamma=GAMMA, plasma_prm=None):
        self.edges = np.asarray(edges, dtype=float)
        if (self.edges.ndim != 1 or len(self.edges) < 3
                or not np.isfinite(self.edges).all() or np.any(np.diff(self.edges) <= 0)):
            raise ValueError("Euler edges must be finite and strictly increasing")
        if cylindrical and self.edges[0] != 0:
            raise ValueError("cylindrical domain must start at the axis r=0")
        if gamma <= 1:
            raise ValueError("gamma must exceed one")
        self.cylindrical = cylindrical
        self.gamma = gamma
        self.prm = plasma_prm
        self.centres = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.volumes = math.pi * np.diff(self.edges ** 2) if cylindrical else np.diff(self.edges)
        self.faces = 2 * math.pi * self.edges if cylindrical else np.ones(len(self.edges))
        self.boundary_out = np.zeros(3)

    def primitive(self, state):
        if state.shape != (9 if self.prm is not None else 3, len(self.centres)):
            raise ValueError("invalid conserved-state shape")
        rho = state[0]
        if not np.isfinite(state).all() or np.any(rho <= 0):
            raise ArithmeticError("non-finite state or nonpositive mass density")
        velocity = state[1] / rho
        internal = state[2] - 0.5 * rho * velocity ** 2
        if self.prm is not None:
            if np.any(state[3:] < 0):
                raise ArithmeticError("negative plasma tracer")
            occupied = state[3] + 2 * state[4] + state[6] + 2 * state[7] + state[8]
            if np.any(occupied > rho / M_HE * (1 + 1e-12)):
                raise ArithmeticError("plasma tracers exceed the available helium nuclei")
            internal -= chemical_energy(state[3:], self.prm)
            gas_internal = internal - state[5] * k3.QE
        else:
            gas_internal = internal
        if np.any(gas_internal <= 0):
            raise ArithmeticError("nonpositive heavy-gas internal energy")
        pressure = (self.gamma - 1) * internal
        return rho, velocity, pressure, gas_internal

    def stable_dt(self, state, cfl=0.25):
        if not 0 < cfl <= 0.4:
            raise ValueError("CFL must lie in (0, 0.4]")
        rho, velocity, pressure, _ = self.primitive(state)
        speed = np.abs(velocity) + np.sqrt(self.gamma * pressure / rho)
        return float(cfl * np.min(np.diff(self.edges) / speed))

    def _flux(self, state):
        rho, velocity, pressure, _ = self.primitive(state)
        flux = state * velocity
        flux[1] += pressure
        flux[2] += pressure * velocity
        return flux, velocity, np.sqrt(self.gamma * pressure / rho)

    def step(self, state, dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Euler timestep must be positive and finite")
        if dt > self.stable_dt(state, 0.4) * (1 + 1e-12):
            raise ValueError("Euler timestep exceeds the supported CFL")
        _, velocity, pressure, _ = self.primitive(state)
        left = state[:, :1].copy()
        if self.cylindrical:
            left[1] *= -1
        padded = np.concatenate([left, state, state[:, -1:]], axis=1)
        # Reuse primitive/flux arithmetic without changing the physical mesh.
        ul, ur = padded[:, :-1], padded[:, 1:]
        def flux_sides(u):
            rho = u[0]
            vel = u[1] / rho
            eint = u[2] - 0.5 * rho * vel ** 2
            if self.prm is not None:
                eint -= chemical_energy(u[3:], self.prm)
            p = (self.gamma - 1) * eint
            f = u * vel
            f[1] += p
            f[2] += p * vel
            return f, vel, np.sqrt(self.gamma * p / rho)
        fl, vl, cl = flux_sides(ul)
        fr, vr, cr = flux_sides(ur)
        sl = np.minimum(vl - cl, vr - cr)
        sr = np.maximum(vl + cl, vr + cr)
        flux = (sr * fl - sl * fr + sl * sr * (ur - ul)) / (sr - sl)
        flux = np.where(sl >= 0, fl, np.where(sr <= 0, fr, flux))
        result = state - dt * np.diff(flux * self.faces, axis=1) / self.volumes
        if self.cylindrical:
            result[1] += dt * pressure * np.diff(self.faces) / self.volumes
        if self.prm is not None:
            face_velocity = 0.5 * (vl + vr)
            if self.cylindrical:
                face_velocity[0] = 0
            divergence = np.diff(self.faces * face_velocity) / self.volumes
            result[5] *= np.exp(-(self.gamma - 1) * divergence * dt)
        self.primitive(result)
        self.boundary_out += dt * (flux[:3, -1] * self.faces[-1] - flux[:3, 0] * self.faces[0])
        return result

    def integrate(self, state, t_end, samples):
        if t_end <= 0 or samples < 2:
            raise ValueError("positive duration and at least two samples required")
        targets = np.linspace(0, t_end, samples)
        frames = [state.copy()]
        t = 0.0
        for target in targets[1:]:
            while t < target * (1 - 1e-13):
                dt = min(target - t, self.stable_dt(state))
                state = self.step(state, dt)
                t += dt
            frames.append(state.copy())
        return targets, np.asarray(frames)


def sod_star(gamma=1.4):
    """Exact pressure/contact velocity for the standard Sod shock tube."""
    def f(p, rho, p0):
        a = math.sqrt(gamma * p0 / rho)
        if p > p0:
            return (p - p0) * math.sqrt(2 / ((gamma + 1) * rho)
                                        / (p + (gamma - 1) / (gamma + 1) * p0))
        return 2 * a / (gamma - 1) * ((p / p0) ** ((gamma - 1) / (2 * gamma)) - 1)
    lo, hi = 0.1, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid, 1, 1) + f(mid, 0.125, 0.1) > 0:
            hi = mid
        else:
            lo = mid
    pstar = 0.5 * (lo + hi)
    ustar = 0.5 * (f(pstar, 0.125, 0.1) - f(pstar, 1, 1))
    shock_speed = math.sqrt(gamma * 0.1 / 0.125) * math.sqrt(
        (gamma + 1) / (2 * gamma) * pstar / 0.1 + (gamma - 1) / (2 * gamma))
    return pstar, ustar, shock_speed
