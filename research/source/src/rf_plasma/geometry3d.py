"""Staircased 3D geometries and TE10 port helpers for the custom FDTD solver."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import epsilon_0 as EPS0, mu_0 as MU0

from .fdtd3d import Grid, te10_discrete_beta


@dataclass
class Layout:
    """Everything the coupled driver needs to know about one geometry."""

    name: str
    grid: Grid
    node_inside: np.ndarray
    pml_lo: int
    pml_hi: int
    rect_i0: int          # x node index of the WR-284 feed's left wall
    rect_j0: int          # y node index of the feed's bottom wall
    rect_na: int          # feed width in cells
    rect_nb: int          # feed height in cells
    source_k: int
    meas_k: int
    trans_k: int | None   # Poynting plane on the output side (None for a closed cavity)
    energy_k_lo: int      # energy is integrated from here to the end of the non-PML domain
    energy_k_hi: int
    focus_i: int
    focus_k: int
    focus_j: np.ndarray   # Ey edge indices crossed by the laser (laser along -y)
    description: str
    gas_k_lo: int = 0     # gas (and plasma) occupies node planes gas_k_lo..gas_k_hi
    gas_k_hi: int = -1

    @property
    def feed_width_m(self) -> float:
        return self.rect_na * self.grid.dx

    @property
    def feed_height_m(self) -> float:
        return self.rect_nb * self.grid.dx


def _cells(length_m: float, dx: float) -> int:
    return max(1, int(round(length_m / dx)))


def iris_cavity(dx: float, width_m: float, height_m: float, cavity_length_m: float,
                iris_aperture_m: float, iris_thickness_m: float = 2e-3,
                feed_length_m: float = 0.045, pml_cells: int = 16) -> Layout:
    """WR-284 feed, symmetric inductive iris, closed rectangular cavity (PEC everywhere)."""
    na, nb = _cells(width_m, dx), _cells(height_m, dx)
    k_src = pml_cells + 4
    k_meas = k_src + _cells(0.02, dx)
    k_iris0 = pml_cells + _cells(feed_length_m, dx)
    k_iris1 = k_iris0 + _cells(iris_thickness_m, dx)
    k_end = k_iris1 + _cells(cavity_length_m, dx)
    nz = k_end
    grid = Grid(na, nb, nz, dx)
    ins = np.ones((na + 1, nb + 1, nz + 1), bool)
    x = np.arange(na + 1) * dx - na * dx / 2
    outside_aperture = np.abs(x) > iris_aperture_m / 2 + 1e-12
    for k in range(k_iris0, k_iris1 + 1):
        ins[outside_aperture, :, k] = False
    focus_k = k_iris1 + (k_end - k_iris1) // 2
    return Layout(
        name="iris_cavity", grid=grid, node_inside=ins, pml_lo=pml_cells, pml_hi=0,
        rect_i0=0, rect_j0=0, rect_na=na, rect_nb=nb, source_k=k_src, meas_k=k_meas,
        trans_k=None, energy_k_lo=k_meas, energy_k_hi=nz, focus_i=na // 2, focus_k=focus_k,
        focus_j=np.arange(nb), description=(
            f"WR-284 ({na*dx*1e3:.0f}x{nb*dx*1e3:.0f} mm) feed, {iris_aperture_m*1e3:.1f} mm inductive iris, "
            f"{(k_end-k_iris1)*dx*1e3:.0f} mm cavity; laser along y at cavity centre"),
        gas_k_lo=k_iris1 + 1, gas_k_hi=nz - 1,
    )


def straight_guide(dx: float, width_m: float, height_m: float, length_m: float, pml_cells: int = 16) -> Layout:
    """Matched WR-284 section with CPML at both ends (port/absorber validation)."""
    na, nb = _cells(width_m, dx), _cells(height_m, dx)
    nz = 2 * pml_cells + _cells(length_m, dx)
    grid = Grid(na, nb, nz, dx)
    ins = np.ones((na + 1, nb + 1, nz + 1), bool)
    k_src = pml_cells + 4
    k_meas = k_src + _cells(0.02, dx)
    return Layout(
        name="straight_guide", grid=grid, node_inside=ins, pml_lo=pml_cells, pml_hi=pml_cells,
        rect_i0=0, rect_j0=0, rect_na=na, rect_nb=nb, source_k=k_src, meas_k=k_meas,
        trans_k=nz - pml_cells - 4, energy_k_lo=k_meas, energy_k_hi=nz - pml_cells - 4,
        focus_i=na // 2, focus_k=(k_meas + nz - pml_cells) // 2, focus_j=np.arange(nb),
        description="matched straight WR-284 guide",
        gas_k_lo=k_meas + 2, gas_k_hi=nz - pml_cells - 6,
    )


def _taper_limit(theta: np.ndarray, fraction: np.ndarray, radius: float, a: float, b: float) -> np.ndarray:
    """Same smooth-step polar interpolation as transition_rf.transition_polygon (0=circle, 1=rect)."""
    c, s = np.abs(np.cos(theta)), np.abs(np.sin(theta))
    with np.errstate(divide="ignore"):
        xl = np.where(c > 1e-15, a / 2 / np.maximum(c, 1e-300), np.inf)
        yl = np.where(s > 1e-15, b / 2 / np.maximum(s, 1e-300), np.inf)
    rect = np.minimum(xl, yl)
    smooth = fraction**2 * (3 - 2 * fraction)
    return (1 - smooth) * radius + smooth * rect


def taper_guide(dx: float, width_m: float, height_m: float, diameter_m: float, taper_length_m: float,
                feed_length_m: float = 0.045, circular_length_m: float = 0.09, pml_cells: int = 16) -> Layout:
    """WR-284 feed -> smooth taper -> 90 mm circular guide terminated by CPML."""
    radius = diameter_m / 2
    n = _cells(diameter_m, dx)
    na, nb = _cells(width_m, dx), _cells(height_m, dx)
    i0, j0 = (n - na) // 2, (n - nb) // 2
    k_t0 = pml_cells + _cells(feed_length_m, dx)
    k_t1 = k_t0 + _cells(taper_length_m, dx)
    k_out = k_t1 + _cells(circular_length_m, dx)
    nz = k_out + pml_cells
    grid = Grid(n, n, nz, dx)
    xc = n * dx / 2
    x = np.arange(n + 1) * dx - xc
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.hypot(X, Y)
    theta = np.arctan2(Y, X)
    ins = np.zeros((n + 1, n + 1, nz + 1), bool)
    xr = (np.arange(n + 1) - i0) * dx
    yr = (np.arange(n + 1) - j0) * dx
    rect_mask = ((xr >= -1e-12) & (xr <= na * dx + 1e-12))[:, None] & ((yr >= -1e-12) & (yr <= nb * dx + 1e-12))[None, :]
    a_eff, b_eff = na * dx, nb * dx
    for k in range(nz + 1):
        if k <= k_t0:
            ins[:, :, k] = rect_mask
        elif k >= k_t1:
            ins[:, :, k] = rho <= radius + 1e-9
        else:
            frac = (k_t1 - k) / (k_t1 - k_t0)
            limit = _taper_limit(theta, np.full_like(theta, frac), radius, a_eff, b_eff)
            ins[:, :, k] = rho <= limit + 1e-9
    focus_k = k_t1 + (k_out - k_t1) // 2
    focus_i = n // 2
    js = np.array([j for j in range(n) if ins[focus_i, j, focus_k] and ins[focus_i, j + 1, focus_k]])
    return Layout(
        name="taper_guide", grid=grid, node_inside=ins, pml_lo=pml_cells, pml_hi=pml_cells,
        rect_i0=i0, rect_j0=j0, rect_na=na, rect_nb=nb, source_k=pml_cells + 4,
        meas_k=pml_cells + 4 + _cells(0.02, dx), trans_k=k_out - 4,
        energy_k_lo=pml_cells + 4 + _cells(0.02, dx), energy_k_hi=k_out - 4,
        focus_i=focus_i, focus_k=focus_k, focus_j=js, description=(
            f"WR-284 feed, {taper_length_m*1e3:.0f} mm smooth taper, {diameter_m*1e3:.0f} mm circular guide "
            "with CPML termination; laser along y at circular-guide centre"),
        gas_k_lo=k_t0, gas_k_hi=k_out - 6,
    )


class TE10Port:
    """Soft TE10 current sheet plus modal (forward/backward) projection at a reference plane."""

    def __init__(self, layout: Layout, frequency_hz: float):
        self.layout = layout
        g = layout.grid
        self.f = frequency_hz
        nxp1, ny = g.nx + 1, g.ny
        prof = np.zeros((nxp1, ny))
        i = np.arange(nxp1)
        xi = (i - layout.rect_i0) * g.dx
        sx = np.where((i > layout.rect_i0) & (i < layout.rect_i0 + layout.rect_na),
                      np.sin(math.pi * xi / layout.feed_width_m), 0.0)
        jj = np.arange(ny)
        sy = ((jj >= layout.rect_j0) & (jj < layout.rect_j0 + layout.rect_nb)).astype(float)
        prof[:] = sx[:, None] * sy[None, :]
        self.profile = prof
        self.norm = float(np.sum(prof**2))   # discrete int sin^2 dA / dx^2
        self.beta = te10_discrete_beta(frequency_hz, layout.feed_width_m, g.dx, g.dt)
        omega = 2 * math.pi * frequency_hz
        # discrete TE10 wave impedance from the staggered Faraday update
        self.z_te = MU0 * (math.sin(omega * g.dt / 2) / (g.dt / 2)) / (
            math.sin(self.beta * g.dx / 2) / (g.dx / 2))
        self.amplitude = 0.0
        self.envelope = lambda t: 1.0

    def inject(self, fdtd):
        g = self.layout.grid
        t = (fdtd.step_count + 1) * g.dt
        val = self.amplitude * self.envelope(t) * math.sin(2 * math.pi * self.f * t)
        if val != 0.0:
            fdtd.Ey[:, :, self.layout.source_k] -= (g.dt / EPS0) * val * self.profile

    def modal_amplitudes(self, v: np.ndarray, cur: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Complex forward/backward TE10 E-amplitudes from per-step projections.

        v[n]   : sum Ey(t_n)  * profile        (E at t_n = m dt after step m)
        cur[n] : sum Hx(t_n - dt/2) * profile  (H from the same step, averaged to the E plane)
        Returns sliding one-period phasors (len = len(v) - period + 1).
        """
        g = self.layout.grid
        omega = 2 * math.pi * self.f
        period = int(round(1 / (self.f * dt)))
        n = np.arange(len(v))
        ph_e = v * np.exp(-1j * omega * n * dt)
        ph_h = cur * np.exp(-1j * omega * (n - 0.5) * dt)
        kern = np.ones(period) * 2 / period
        ve = np.convolve(ph_e, kern, mode="valid") / self.norm
        vh = np.convolve(ph_h, kern, mode="valid") / self.norm
        vh = vh / math.cos(self.beta * g.dx / 2)  # spatial average of Hx(k-1/2), Hx(k+1/2)
        forward = 0.5 * (ve - self.z_te * vh)
        backward = 0.5 * (ve + self.z_te * vh)
        return forward, backward

    def power(self, amplitude: np.ndarray) -> np.ndarray:
        """Cycle-averaged power (W) carried by a TE10 wave of peak E amplitude |A|."""
        g = self.layout.grid
        return np.abs(amplitude) ** 2 * self.norm * g.dx**2 / (2 * self.z_te)


def raised_cosine_ramp(t_ramp: float):
    def env(t: float) -> float:
        if t >= t_ramp:
            return 1.0
        return 0.5 * (1 - math.cos(math.pi * t / t_ramp))
    return env
