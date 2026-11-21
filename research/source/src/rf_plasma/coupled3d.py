"""Cycle-resolved 3D RF-plasma coupling: FDTD + bulk local-field kinetics + laser filament.

Per block of `block` FDTD steps:
  FDTD steps (TE10 drive, probes) -> bulk RF work to nodes -> bulk kinetics (Boltzmann tables)
  -> new bulk Drude coefficients -> filament kinetics/hydro/transport with laser (Maxwell tables)
  -> new multi-ring Drude coefficients on the filament edges.
The plasma current memory is retained on every edge and ring (no envelope approximation).

Subcommands: calibrate (cavity ringdown vs iris aperture), prefill (cold steady state, scaled to
the 10 MV/m focus target), run (one case), matrix (case matrix + controls).
"""

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import time

import numba
import numpy as np
from scipy.signal import hilbert

from . import plasma3d as k3
from . import plasma3d_cr as kc
from .fdtd3d import FDTD, _plane_projection, _poynting_z
from .gas_model import load_gas_model
from .geometry3d import TE10Port, iris_cavity, raised_cosine_ramp, straight_guide, taper_guide

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "coupled3d"
REF = json.loads((ROOT / "cases" / "reference.json").read_text(encoding="utf-8"))
KB = 1.380649e-23
QE = 1.602176634e-19
C0 = 299792458.0
F_NOMINAL = REF["rf"]["frequency_hz"]
GUIDE_LENGTH_M = 0.20  # matched straight WR-284 section between the CPML absorbers


# ---------------------------------------------------------------- geometry / calibration
def make_layout(geometry: str, dx: float, aperture_m: float | None = None):
    r = REF["rectangular"]
    if geometry == "cavity":
        if aperture_m is None:
            raise ValueError("cavity needs an iris aperture (run calibrate)")
        return iris_cavity(dx, r["width_m"], r["height_m"], r["closed_cavity_length_m"], aperture_m)
    if geometry == "taper":
        c = REF["circular"]
        return taper_guide(dx, r["width_m"], r["height_m"], c["diameter_m"], c["taper_length_m"])
    if geometry == "guide":
        return straight_guide(dx, r["width_m"], r["height_m"], GUIDE_LENGTH_M)
    raise ValueError(geometry)


def _cal_key(geometry, dx):
    return f"{geometry}_dx{dx * 1e3:g}mm"


def _load_cal() -> dict:
    p = RESULTS / "calibration.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _save_cal(entries: dict):
    """Merge top-level entries into calibration.json (re-read at save time; safe for parallel jobs)."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    cal = _load_cal()
    cal.update(entries)
    (RESULTS / "calibration.json").write_text(json.dumps(cal, indent=2), encoding="utf-8")


def phasor(sig: np.ndarray, dt: float, f: float, periods: float = 3.0) -> complex:
    """Least-squares complex amplitude P of sig ~ Re(P exp(i w t)) + c over the last periods."""
    n = min(len(sig), int(round(periods / (f * dt))))
    t = (np.arange(len(sig) - n, len(sig))) * dt
    w = 2 * math.pi * f
    a = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones(n)])
    c = np.linalg.lstsq(a, sig[-n:], rcond=None)[0]
    return complex(c[0], -c[1])


def _focus_edge(lay):
    js = lay.focus_j
    jc = int(js[np.argmin(np.abs((js + 0.5) - lay.grid.ny / 2))])
    return lay.focus_i, jc, lay.focus_k


def ringdown(dx: float, aperture_m: float, t_ring_s: float = 30e-9) -> dict:
    lay = make_layout("cavity", dx, aperture_m)
    g = lay.grid
    fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, F_NOMINAL)
    port = TE10Port(lay, F_NOMINAL)
    port.amplitude = 1.0
    t_c, sig = 2.5e-9, 0.6e-9
    port.envelope = lambda t: math.exp(-0.5 * ((t - t_c) / sig) ** 2)
    i, j, k = _focus_edge(lay)
    n = int((2 * t_c + t_ring_s) / g.dt)
    ey = np.empty(n)
    for s in range(n):
        fd.step(port.inject)
        ey[s] = fd.Ey[i, j, k]
    t = np.arange(n) * g.dt
    sel = t > 2 * t_c + 3e-9
    an = hilbert(ey[sel])
    env, ph = np.log(np.abs(an)), np.unwrap(np.angle(an))
    tt = t[sel]
    trim = slice(len(tt) // 10, -len(tt) // 10)
    slope = np.polyfit(tt[trim], env[trim], 1)[0]
    f0 = np.polyfit(tt[trim], ph[trim], 1)[0] / (2 * math.pi)
    tau = -1.0 / slope
    return dict(aperture_m=aperture_m, f0_hz=float(f0), tau_amp_s=float(tau), q_loaded=float(math.pi * f0 * tau),
                cells=g.cells, dt_s=g.dt)


def calibrate(dx: float, apertures_mm, q_target: float):
    rows = []
    for a in apertures_mm:
        t0 = time.time()
        r = ringdown(dx, a * 1e-3)
        r["wall_s"] = time.time() - t0
        rows.append(r)
        print(f"aperture {a:5.1f} mm: f0 {r['f0_hz'] / 1e9:.5f} GHz  Q_L {r['q_loaded']:.0f}  ({r['wall_s']:.0f} s)", flush=True)
    best = min(rows, key=lambda r: abs(math.log(r["q_loaded"] / q_target)))
    key = _cal_key("cavity", dx)
    entry = _load_cal().get(key, {})
    entry.update(scan=rows, chosen=best, q_target=q_target,
                 note="PEC walls: Q_L is external (iris) only; illustrative choice of coupling")
    _save_cal({key: entry})
    return best


def prefill(geometry: str, dx: float, threads: int):
    """Drive the empty device to steady state at unit source amplitude; store state and focus field."""
    numba.set_num_threads(threads)
    cal = _load_cal()
    key = _cal_key(geometry, dx)
    if geometry == "cavity":
        ch = cal[key]["chosen"]
        lay = make_layout(geometry, dx, ch["aperture_m"])
        f0, tau = ch["f0_hz"], ch["tau_amp_s"]
        t_boost = tau * math.log(2.0)          # single-mode exact fill: 2x drive until a = a_ss
        t_total = t_boost + 2.5 * tau
        ramp = raised_cosine_ramp(1e-9)
        env = lambda t: ramp(t) * (2.0 if t < t_boost else 1.0)
    else:
        lay = make_layout(geometry, dx)
        f0, t_total = F_NOMINAL, 10e-9
        env = raised_cosine_ramp(2e-9)
    g = lay.grid
    fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, f0)
    port = TE10Port(lay, f0)
    port.amplitude = 1.0
    port.envelope = env
    i, j, k = _focus_edge(lay)
    n = int(t_total / g.dt)
    ey = np.empty(n)
    t0 = time.time()
    for s in range(n):
        fd.step(port.inject)
        ey[s] = fd.Ey[i, j, k]
    per = int(round(1 / (f0 * g.dt)))
    p_end = abs(phasor(ey, g.dt, f0))
    p_prev = abs(phasor(ey[: -int(0.25 * (tau if geometry == "cavity" else 2e-9) / g.dt)], g.dt, f0))
    settle = abs(p_end / p_prev - 1.0)
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"prefill_{key}.npz"
    np.savez_compressed(path, **fd.state(), ey_focus=ey[-20 * per:])
    entry = _load_cal().get(key, {})
    entry["prefill"] = dict(file=path.name, f_drive_hz=f0, steps=n, focus_peak_per_unit=p_end,
                            amplitude_for_target=REF["rf"]["empty_focus_field_peak_v_m"] / p_end,
                            settle_rel_change=settle, wall_s=time.time() - t0)
    _save_cal({key: entry})
    print(json.dumps(entry["prefill"], indent=2))
    return entry["prefill"]


# ---------------------------------------------------------------- case definition
def _impurity_seed(case, rb, beam_w, p_peak, nt=241):
    """Photo-electrons from a trace impurity (PPT rate of the impurity species) in the vacuum focal beam.

    Returns (count per ring per metre of laser path, excess energy per electron in eV), both (nE, R).
    Saturation (1 - exp(-int w dt)) caps the yield at the impurity density."""
    from .gas_model import load_gas_model
    gi = load_gas_model(case.impurity, case.wavelength_m, allow_unverified=True)
    ln_i = gi.ln_i0 + gi.hi * np.arange(len(gi.ppt_tab))
    n_imp = case.impurity_frac * case.pressure_pa / (KB * case.temperature_k)
    rc = 0.5 * (rb[:, 1:] + rb[:, :-1])
    area = math.pi * np.diff(rb ** 2, axis=1)
    tr = np.linspace(-2.0, 2.0, nt) * case.fwhm_s
    pt = p_peak * np.exp(-4 * math.log(2) * (tr / case.fwhm_s) ** 2)
    i_rt = (2 * pt[None, None, :] / (math.pi * beam_w[:, None, None] ** 2)
            * np.exp(-2 * rc[:, :, None] ** 2 / beam_w[:, None, None] ** 2))
    ln_w = np.interp(np.log(np.maximum(i_rt, 1e-300)), ln_i, gi.ppt_tab, left=-690.0)
    expo = np.trapezoid(np.exp(ln_w), tr, axis=-1)
    count = n_imp * (1.0 - np.exp(-expo)) * area
    # prescribed laser-localized seed (stands in for unresolved seeding: aerosols, surfaces, metastables).
    # k > 0: effective multiphoton order, n ~ int I^k dt, normalized to laser_seed_m3 at the reference pulse focus
    if case.laser_seed_m3 > 0:
        k = case.laser_seed_k
        if k > 0:
            e_r, t_r, w_r = case.laser_seed_ref
            i_pk = 2 * p_peak / (math.pi * beam_w ** 2)
            i_ref = 2 * e_r * math.sqrt(4 * math.log(2) / math.pi) / t_r / (math.pi * w_r ** 2)
            amp = (i_pk / i_ref) ** k * (case.fwhm_s / t_r)
        else:
            k, amp = 1.0, np.ones_like(beam_w)
        count += case.laser_seed_m3 * amp[:, None] * np.exp(-2 * k * rc ** 2 / beam_w[:, None] ** 2) * area
    ph = gi.prm[k3.P_PHOTON_EV]
    excess = max(math.ceil(gi.prm[k3.P_IP] / ph) * ph - gi.prm[k3.P_IP], 0.05)
    return count, np.full_like(count, excess)


@dataclass
class Case:
    geometry: str = "cavity"
    gas: str = "Ar"
    laser: bool = True
    rf: bool = True
    fwhm_s: float = 300e-12
    t_laser_s: float = 1.5e-9
    t_end_s: float = 12e-9
    energy_j: float = REF["laser"]["energy_at_focus_j"]
    waist_m: float = REF["laser"]["waist_radius_m"]
    wavelength_m: float = REF["laser"]["wavelength_m"]
    field_target_v_m: float = REF["rf"]["empty_focus_field_peak_v_m"]
    pressure_pa: float = REF["gas"]["gauge_pressure_pa"] + REF["gas"]["ambient_pressure_pa"]
    temperature_k: float = REF["gas"]["temperature_k"]
    seed_m3: float = 1e9
    ee_scale: float = 1.0
    dx_m: float = 1e-3
    rings: int = 48
    block: int = 8
    hydro: bool = True
    allow_unverified: bool = True
    snapshot_s: list = field(default_factory=list)
    frame_every_blocks: int = 0  # >0: store animation frames of computed fields/plasma every N blocks
    tag: str = ""
    # extended physics (plasma3d_cr); all False reproduces the original kernels exactly
    cr: bool = False            # lumped He* / He2* / He2+(alpha) collisional-radiative chemistry
    transport: bool = False     # bulk free/ambipolar diffusion, electron energy convection/conduction
    gas_heating: bool = False   # bulk isochoric gas heating + gas heat conduction
    propagation: bool = False   # radial paraxial laser field: diffraction + absorption (always, in the
                                # extended path) and, with this flag, plasma/neutral/Kerr refraction
    laser_dr_m: float = 1e-6
    laser_rmax_m: float = 0.8e-3
    ee_tables: bool = False     # solved e-e Boltzmann closure (swarm_ee) instead of the w-blend
    fil_free: bool = False      # filament Allis-Rose free -> ambipolar radial transport (ionization-front growth)
    impurity: str = "N2"        # trace species photoionized by the laser (PPT, multiphoton)
    impurity_frac: float = 0.0  # impurity mole fraction (UHP He 99.999 %: <= 1e-5 total impurities)
    laser_seed_m3: float = 0.0  # prescribed peak electron density deposited in the beam at t_laser (Gaussian)
    laser_seed_k: float = 0.0   # effective order: seed ~ int I^k dt relative to laser_seed_ref (0: fixed profile)
    laser_seed_ref: tuple = (1e-3, 1e-9, 20e-6)  # reference pulse (energy J, FWHM s, waist m) for laser_seed_m3

    @property
    def extended(self) -> bool:
        return self.cr or self.transport or self.gas_heating or self.propagation or self.ee_tables or self.fil_free

    @property
    def name(self) -> str:
        mode = ("rf+laser" if self.rf and self.laser else "rf-only" if self.rf else "laser-only")
        pulse = f"_{self.fwhm_s * 1e12:.0f}ps" if self.laser else ""
        return f"{self.geometry}_{self.gas}_{mode}{pulse}{('_' + self.tag) if self.tag else ''}"


def default_case(geometry, gas, fwhm_s, laser=True, rf=True, **kw) -> Case:
    t_l = 1.5e-9 if fwhm_s < 1e-9 else 6e-9
    t_end = t_l + 10e-9
    snaps = [t_l - fwhm_s, t_l, t_l + fwhm_s, t_l + 1e-9, t_l + 3e-9, t_end]
    if not laser:
        snaps = [2e-9, 5e-9, 10e-9, t_end]
    return Case(geometry=geometry, gas=gas, laser=laser, rf=rf, fwhm_s=fwhm_s, t_laser_s=t_l,
                t_end_s=t_end, snapshot_s=snaps, **kw)


def _ring_edges(rings: int, r_out: float, r_uniform: float = 60e-6) -> np.ndarray:
    nu = rings // 2
    inner = np.linspace(0.0, r_uniform, nu + 1)
    outer = np.geomspace(r_uniform, r_out, rings - nu + 1)[1:]
    return np.concatenate([inner, outer])


# ---------------------------------------------------------------- run
def run_case(case: Case, threads: int = 8, out_root: Path = RESULTS) -> dict:
    numba.set_num_threads(threads)
    t_wall = time.time()
    dx = case.dx_m
    cal = _load_cal()
    key = _cal_key(case.geometry, dx)
    aperture = cal[key]["chosen"]["aperture_m"] if case.geometry == "cavity" else None
    lay = make_layout(case.geometry, dx, aperture)
    g = lay.grid
    dt = g.dt
    f_drive = cal[key]["prefill"]["f_drive_hz"] if case.rf else F_NOMINAL
    fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, f_drive)
    port = TE10Port(lay, f_drive)
    if case.rf:
        pre = cal[key]["prefill"]
        scale = case.field_target_v_m / pre["focus_peak_per_unit"]
        st = dict(np.load(RESULTS / pre["file"]))
        st.pop("ey_focus")
        for n_, v in st.items():
            if n_ != "step_count":
                st[n_] = v * scale
        fd.load_state(st)
        port.amplitude = scale
        port.envelope = lambda t: 1.0
    step0 = fd.step_count

    gm = load_gas_model(case.gas, case.wavelength_m, allow_unverified=case.allow_unverified,
                        ee_scale=case.ee_scale, extended=case.extended, cr=case.cr, transport=case.transport,
                        ee_n_gas=(case.pressure_pa / (KB * case.temperature_k)) if case.ee_tables else None)
    tabs = gm.tabs
    xt = case.extended
    if case.fil_free:
        gm.prm[kc.P_FIL_FREE] = 1.0
    n0 = case.pressure_pa / (KB * case.temperature_k)
    tg0 = case.temperature_k
    p0 = case.pressure_pa
    eps_g = 1.5 * KB * tg0 / QE
    active = lay.node_inside.copy()
    active[:, :, : lay.gas_k_lo] = False
    active[:, :, lay.gas_k_hi + 1:] = False
    shape = lay.node_inside.shape
    st6 = np.zeros((6,) + shape)
    nxp, nx2p, wb, n_star, n_exc, n_alpha = st6  # views: the kernels update st6 in place
    nxp[active] = case.seed_m3
    wb[:] = nxp * eps_g
    tg_b = np.full(shape, tg0)
    wp2 = np.zeros(shape)
    nu = np.zeros(shape)
    node_work = np.zeros(shape)
    bulk_led = np.zeros((shape[0], 9))
    if xt:
        coef_b = np.zeros((6,) + shape)
        delta_b = np.zeros((5,) + shape)
        d_max = 1.0
        bulk_e0 = float(wb.sum()) * QE * dx ** 3
        heat_b0 = float(tg_b[active].sum())

    fi, fj, fk = _focus_edge(lay)
    js = np.asarray(lay.focus_j, dtype=np.int64)
    ne_ = len(js)
    fd.attach_multipole(fi, fk, js, case.rings)
    mp = fd.multipole
    r_out = dx / math.sqrt(math.pi)
    rb = np.tile(_ring_edges(case.rings, r_out, min(max(60e-6, 1.5 * case.waist_m), 0.5 * r_out)), (ne_, 1))
    ub = np.zeros_like(rb)
    area = math.pi * np.diff(rb ** 2, axis=1)
    heavy = n0 * area
    eg = 0.5 * gm.prm[k3.P_DOF] * KB * tg0 * heavy
    cs6 = np.zeros((6,) + area.shape)
    cxp, cx2p, ue = cs6[0], cs6[1], cs6[2]
    cxp[:] = case.seed_m3 * area
    ue[:] = cxp * eps_g
    fil_led = np.zeros((ne_, 12))
    mass_r = heavy * gm.prm[k3.P_MASS]

    def fil_energy():
        mb = np.zeros_like(ub); mb[:, 1:] += 0.5 * mass_r; mb[:, 1:-1] += 0.5 * mass_r[:, 1:]
        return float((ue.sum() * QE + eg.sum() + 0.5 * np.sum(mb * ub ** 2)) * dx)

    e_fil0 = fil_energy()
    laser_in = np.zeros(ne_)
    trans = np.ones(1)
    k3.filament_drude(rb, cxp, cx2p, ue, heavy, dt, dx, gm.prm, tabs, gm.ln_e0, gm.h,
                      mp["w"], mp["a"], mp["b"])

    # laser geometry (propagation along -y, focus at the guide centre)
    y_e = (js + 0.5) * dx
    y_f = g.ny * dx / 2
    z_r = math.pi * case.waist_m ** 2 / case.wavelength_m
    beam_w = case.waist_m * np.sqrt(1 + ((y_e - y_f) / z_r) ** 2)
    order = np.argsort(-y_e).astype(np.int64)
    p_peak = case.energy_j * math.sqrt(4 * math.log(2) / math.pi) / case.fwhm_s

    def beam_power(t):  # W, unattenuated; quasi-static across the guide (transit ~0.1 ns << RF/cavity scales)
        tr = t - case.t_laser_s
        return np.full(ne_, p_peak * math.exp(-4 * math.log(2) * (tr / case.fwhm_s) ** 2)) if case.laser \
            else np.zeros(ne_)

    if xt:
        # radial paraxial field entering edge order[0], vacuum Gaussian normalized to 1 W
        npg = int(round(case.laser_rmax_m / case.laser_dr_m))
        r_g = (np.arange(npg) + 0.5) * case.laser_dr_m
        dA_g = math.pi * case.laser_dr_m ** 2 * (2 * np.arange(npg) + 1)
        k_l = 2 * math.pi / case.wavelength_m
        q0 = -(y_e[order[0]] - y_f) - 1j * z_r
        field0 = np.exp(1j * k_l * r_g ** 2 / (2 * q0))
        field0 /= math.sqrt(float(np.sum(np.abs(field0) ** 2 * dA_g)))
        n_abs = max(4, int(0.2 * npg))
        lat_mask = np.ones(npg)
        lat_mask[-n_abs:] = np.cos(0.5 * math.pi * np.linspace(0, 1, n_abs) ** 2) ** 0.02
        field = np.zeros(npg, np.complex128)
        omega_l = 2 * math.pi * C0 / case.wavelength_m
        n_crit = k3.EPS0 * k3.ME * omega_l ** 2 / QE ** 2 if case.propagation else 1e300
        refr = gm.optics.get("refractivity_stp", 0.0) * n0 / k3.N_LOSCHMIDT if case.propagation else 0.0
        n2 = gm.optics.get("n2_m2_W_per_atm", 0.0) * p0 / 101325.0 if case.propagation else 0.0
        beam_out = np.zeros(2 + 2 * ne_)
        laser_lat = 0.0
    fil_drude = (lambda: kc.filament_drude_cr(rb, cs6, heavy, dt, dx, gm.prm, tabs, gm.ln_e0, gm.h,
                                              mp["w"], mp["a"], mp["b"])) if xt else \
        (lambda: k3.filament_drude(rb, cxp, cx2p, ue, heavy, dt, dx, gm.prm, tabs, gm.ln_e0, gm.h,
                                   mp["w"], mp["a"], mp["b"]))
    imp_seed = _impurity_seed(case, rb, beam_w, p_peak) \
        if (case.laser and (case.impurity_frac > 0 or case.laser_seed_m3 > 0)) else None
    imp_done = imp_seed is None

    dtb = case.block * dt
    nblk = int(math.ceil(case.t_end_s / dtb))
    nst = nblk * case.block
    ey_f = np.zeros(nst); v_in = np.zeros(nst); i_in = np.zeros(nst); s_out = np.zeros(nst)
    ts = {k_: np.zeros(nblk) for k_ in (
        "t", "em_energy", "bulk_electrons", "bulk_ne_max", "bulk_te_at_max", "fil_electrons", "rf_to_bulk",
        "rf_to_fil", "rf_to_nodes", "laser_delivered", "laser_absorbed", "laser_exit",
        "fil_ne_axis_focus", "fil_te_axis_focus", "fil_tg_axis_focus", "fil_r50_focus", "fil_rb_out_focus")}
    sp = {k_: np.zeros((nblk, ne_), np.float32) for k_ in ("fil_ne_axis", "fil_te_axis", "fil_line_e")}
    bulk_line = np.zeros((nblk, shape[1]), np.float32)
    tx = {k_: np.zeros(nblk) for k_ in (
        "bulk_tg_max", "bulk_hestar_max", "bulk_he2star_max", "bulk_alpha_max", "fil_hestar_axis_focus",
        "fil_he2star_axis_focus", "fil_alpha_axis_focus", "laser_lateral", "beam_axis_I_focus",
        "beam_r2nd_focus", "bulk_transport_dmax")} if xt else {}
    spx = {k_: np.zeros((nblk, ne_), np.float32) for k_ in ("beam_axis_I", "beam_r2nd", "fil_hestar_axis")} \
        if xt else {}
    snaps = {}
    frames = {k_: [] for k_ in ("t", "ey_xz", "ey_yz", "ne_xz", "ne_yz", "fil_ne_axis", "fil_te_axis",
                                "fil_r50", "fil_line_e", "fil_x2p_frac", "laser_power", "em_energy",
                                "fil_rb_focus", "fil_ne_r_focus", "fil_te_r_focus", "fil_tg_r_focus",
                                "fil_x2p_r_focus")}
    snap_t = sorted(case.snapshot_s)
    snap_i = 0
    e_focus = int(np.where(js == fj)[0][0])
    vol = dx ** 3
    laser_deliv = laser_abs = laser_exit = 0.0
    rf_nodes = 0.0
    src = port.inject if case.rf else None
    status = "ok"
    try:
        for b in range(nblk):
            if case.rf:
                for s in range(case.block):
                    fd.step(src)
                    q = b * case.block + s
                    ey_f[q] = fd.Ey[fi, fj, fk]
                    v_in[q], i_in[q] = _plane_projection(fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
                    if lay.trans_k is not None:
                        s_out[q] = _poynting_z(fd.Ex, fd.Ey, fd.Hx, fd.Hy, lay.trans_k) * dx * dx
                fd.collect_bulk_work(node_work)
                rf_nodes += float(node_work.sum()) * vol
                if xt:
                    kc.bulk_update_cr(st6, tg_b, active, node_work, dtb, n0, gm.prm, tabs, gm.ln_e0, gm.h,
                                      gm.ln_i0, gm.hi, gm.ppt_tab, wp2, nu, bulk_led, case.gas_heating)
                    if case.transport:
                        t_s = 0.0
                        while t_s < dtb * (1 - 1e-12):
                            dts = min(dtb - t_s, 0.15 * dx * dx / max(d_max, 1e-30))
                            d_max = kc.bulk_transport(st6, tg_b, active, dts, dx, n0, gm.prm, tabs, gm.ln_e0,
                                                      gm.h, coef_b, delta_b, case.gas_heating)
                            t_s += dts
                else:
                    k3.bulk_update(nxp, nx2p, wb, active, node_work, dtb, n0, tg0, gm.prm, tabs, gm.ln_e0,
                                   gm.h, gm.ln_i0, gm.hi, gm.ppt_tab, wp2, nu, bulk_led)
                fd.set_bulk_plasma(wp2, nu)
            t_now = (b + 1) * dtb
            t_mid = t_now - 0.5 * dtb
            if not imp_done and t_now >= case.t_laser_s:
                # impurity photo-electrons (count per ring per length), injected at the pulse centre with
                # the photon excess energy; booked as laser input in the filament ledger
                cxp += imp_seed[0]
                ue += imp_seed[0] * imp_seed[1]
                fil_led[:, 1] += (imp_seed[0] * imp_seed[1]).sum(axis=1) * QE
                imp_done = True
            pw = beam_power(t_mid)
            if xt:
                field[:] = field0 * math.sqrt(pw[order[0]])
                beam_out[:] = 0.0
                kc.filament_update_cr(rb, ub, heavy, eg, cs6, mp["P"], order, laser_in, dtb, n0, tg0, p0,
                                      gm.prm, tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab, field,
                                      case.laser_dr_m, k_l, n2, refr, n_crit, fil_led, case.hydro, dx,
                                      bool(case.laser), lat_mask, beam_out)
                exit_now = beam_out[1] * dtb
                laser_lat += beam_out[0]
            else:
                i_axis = 2 * pw / (math.pi * beam_w ** 2)
                trans[0] = 1.0
                k3.filament_update(rb, ub, heavy, eg, cxp, cx2p, ue, mp["P"], trans, order, laser_in, dtb, n0,
                                   tg0, p0, gm.prm, tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab,
                                   i_axis, beam_w, fil_led, case.hydro, dx)
                exit_now = float(pw[order[-1]] * trans[0]) * dtb
            absorbed = laser_in * dx
            laser_deliv += float(pw[order[0]]) * dtb
            laser_abs += float(absorbed.sum())
            laser_exit += exit_now
            laser_in[:] = 0.0
            fil_drude()
            # ---- diagnostics
            ne_b = nxp + nx2p + 2 * n_alpha
            area = math.pi * np.diff(rb ** 2, axis=1)
            cnt = cxp + cx2p + 2 * cs6[5]
            ne_r = cnt / area
            te_r = np.where(cnt > 0, ue / (1.5 * np.maximum(cnt, 1e-300)), eps_g / 1.5)
            npart = heavy - cx2p - cs6[4]
            tg_r = eg / (0.5 * gm.prm[k3.P_DOF] * KB * np.maximum(npart, 1e-300))
            imax = np.unravel_index(np.argmax(ne_b), shape)
            ts["t"][b] = t_now
            ts["em_energy"][b] = fd.energy(lay.energy_k_lo, lay.energy_k_hi) if case.rf else 0.0
            ts["bulk_electrons"][b] = float(ne_b[active].sum()) * vol
            ts["bulk_ne_max"][b] = float(ne_b[imax])
            ts["bulk_te_at_max"][b] = float(wb[imax] / max(ne_b[imax], 1e-30) / 1.5)
            ts["fil_electrons"][b] = float(cnt.sum()) * dx
            ts["rf_to_bulk"][b] = float(bulk_led[:, 0].sum()) * vol
            ts["rf_to_fil"][b] = float(fil_led[:, 0].sum()) * dx
            ts["rf_to_nodes"][b] = rf_nodes
            ts["laser_delivered"][b] = laser_deliv
            ts["laser_absorbed"][b] = laser_abs
            ts["laser_exit"][b] = laser_exit
            ts["fil_ne_axis_focus"][b] = ne_r[e_focus, 0]
            ts["fil_te_axis_focus"][b] = te_r[e_focus, 0]
            ts["fil_tg_axis_focus"][b] = tg_r[e_focus, 0]
            cum = np.cumsum(cnt[e_focus])
            ts["fil_r50_focus"][b] = rb[e_focus, 1 + int(np.searchsorted(cum, 0.5 * cum[-1]))]
            ts["fil_rb_out_focus"][b] = rb[e_focus, -1]
            sp["fil_ne_axis"][b] = ne_r[:, 0]
            sp["fil_te_axis"][b] = te_r[:, 0]
            sp["fil_line_e"][b] = cnt.sum(axis=1)
            bulk_line[b] = ne_b[fi, :, fk]
            if xt:
                tx["bulk_tg_max"][b] = float(tg_b[active].max())
                tx["bulk_hestar_max"][b] = float(n_star.max())
                tx["bulk_he2star_max"][b] = float(n_exc.max())
                tx["bulk_alpha_max"][b] = float(n_alpha.max())
                tx["fil_hestar_axis_focus"][b] = cs6[3, e_focus, 0] / area[e_focus, 0]
                tx["fil_he2star_axis_focus"][b] = cs6[4, e_focus, 0] / area[e_focus, 0]
                tx["fil_alpha_axis_focus"][b] = cs6[5, e_focus, 0] / area[e_focus, 0]
                tx["laser_lateral"][b] = laser_lat
                tx["beam_axis_I_focus"][b] = beam_out[2 + e_focus]
                tx["beam_r2nd_focus"][b] = beam_out[2 + ne_ + e_focus]
                tx["bulk_transport_dmax"][b] = d_max if case.transport else 0.0
                spx["beam_axis_I"][b] = beam_out[2:2 + ne_]
                spx["beam_r2nd"][b] = beam_out[2 + ne_:]
                spx["fil_hestar_axis"][b] = cs6[3, :, 0] / area[:, 0]
            if case.frame_every_blocks > 0 and b % case.frame_every_blocks == 0:
                cum_all = np.cumsum(cnt, axis=1)
                i50 = np.array([np.searchsorted(cum_all[e_], 0.5 * cum_all[e_, -1]) for e_ in range(ne_)])
                frames["t"].append(t_now)
                frames["ey_xz"].append(fd.Ey[:, fj, :].astype(np.float32))
                frames["ey_yz"].append(fd.Ey[fi, :, :].astype(np.float32))
                frames["ne_xz"].append(ne_b[:, fj, :].astype(np.float32))
                frames["ne_yz"].append(ne_b[fi, :, :].astype(np.float32))
                frames["fil_ne_axis"].append(ne_r[:, 0].astype(np.float32))
                frames["fil_te_axis"].append(te_r[:, 0].astype(np.float32))
                frames["fil_r50"].append(rb[np.arange(ne_), 1 + np.minimum(i50, rb.shape[1] - 2)].astype(np.float32))
                frames["fil_line_e"].append(cnt.sum(axis=1).astype(np.float32))
                frames["fil_x2p_frac"].append((cx2p.sum(axis=1) / np.maximum(cnt.sum(axis=1), 1e-300)).astype(np.float32))
                frames["laser_power"].append(pw.astype(np.float32))
                frames["em_energy"].append(ts["em_energy"][b])
                frames["fil_rb_focus"].append(rb[e_focus].astype(np.float32))
                frames["fil_ne_r_focus"].append(ne_r[e_focus].astype(np.float32))
                frames["fil_te_r_focus"].append(te_r[e_focus].astype(np.float32))
                frames["fil_tg_r_focus"].append(tg_r[e_focus].astype(np.float32))
                frames["fil_x2p_r_focus"].append((cx2p[e_focus] / area[e_focus]).astype(np.float32))
                if xt:
                    frames.setdefault("tg_yz", []).append(tg_b[fi, :, :].astype(np.float32))
                    frames.setdefault("hestar_yz", []).append(n_star[fi, :, :].astype(np.float32))
                    frames.setdefault("fil_hestar_r_focus", []).append((cs6[3, e_focus] / area[e_focus]).astype(np.float32))
                    frames.setdefault("beam_axis_I", []).append(beam_out[2:2 + ne_].astype(np.float32))
                    frames.setdefault("beam_r2nd", []).append(beam_out[2 + ne_:].astype(np.float32))
            while snap_i < len(snap_t) and t_now >= snap_t[snap_i] - 0.5 * dtb:
                tag = f"s{snap_i}"
                snaps[f"{tag}_t"] = t_now
                snaps[f"{tag}_ne_xz"] = ne_b[:, fj, :].astype(np.float32)
                snaps[f"{tag}_ne_yz"] = ne_b[fi, :, :].astype(np.float32)
                snaps[f"{tag}_ey_yz"] = fd.Ey[fi, :, :].astype(np.float32)
                snaps[f"{tag}_ey_xz"] = fd.Ey[:, fj, :].astype(np.float32)
                snaps[f"{tag}_fil_rb"] = rb.astype(np.float32)
                snaps[f"{tag}_fil_ne"] = ne_r.astype(np.float32)
                snaps[f"{tag}_fil_te"] = te_r.astype(np.float32)
                snaps[f"{tag}_fil_tg"] = tg_r.astype(np.float32)
                snap_i += 1
            if b % max(1, nblk // 20) == 0:
                print(f"  {case.name}: t={t_now * 1e9:6.2f} ns  bulk ne_max={ts['bulk_ne_max'][b]:.2e}  "
                      f"fil ne0={ts['fil_ne_axis_focus'][b]:.2e} Te0={ts['fil_te_axis_focus'][b]:.2f} eV  "
                      f"wall {time.time() - t_wall:.0f} s", flush=True)
            if not np.isfinite(ts["bulk_ne_max"][b] + ts["fil_ne_axis_focus"][b]):
                status = f"non-finite state at block {b}"
                break
    except KeyboardInterrupt:
        status = f"interrupted at block {b}"
    nb_done = b + 1

    # ---------------- post-processing: modal powers and RF energy balance
    out = out_root / case.name
    out.mkdir(parents=True, exist_ok=True)
    summary = dict(case=asdict(case), name=case.name, status=status, layout=lay.description,
                   grid=dict(nx=g.nx, ny=g.ny, nz=g.nz, dx=dx, dt=dt, cells=g.cells),
                   f_drive_hz=f_drive, n_gas_m3=n0, filament_edges=int(ne_), wall_s=time.time() - t_wall,
                   provenance=gm.provenance, threads=threads)
    rf = {}
    if case.rf:
        nq = nb_done * case.block
        fwd, bwd = port.modal_amplitudes(v_in[:nq], i_in[:nq], dt)
        per = int(round(1 / (f_drive * dt)))
        p_f, p_b = port.power(fwd), port.power(bwd)
        p_out = np.convolve(s_out[:nq], np.ones(per) / per, mode="valid") if lay.trans_k is not None else np.zeros_like(p_f)
        t_c = (np.arange(len(p_f)) + per / 2) * dt
        e_in = float(np.sum(p_f - p_b - p_out) * dt)
        du = float(ts["em_energy"][nb_done - 1] - ts["em_energy"][0])
        absorbed = float(ts["rf_to_nodes"][nb_done - 1] + ts["rf_to_fil"][nb_done - 1])
        rf = dict(P_forward_W=float(np.median(p_f)), gamma_start=float(abs(bwd[0] / fwd[0])),
                  gamma_end=float(abs(bwd[-1] / fwd[-1])), gamma_min=float(np.min(np.abs(bwd / fwd))),
                  gamma_max=float(np.max(np.abs(bwd / fwd))),
                  energy_in_net_J=e_in, delta_em_J=du, rf_absorbed_J=absorbed,
                  # normalized by the largest term incl. the initial stored energy (absorption may be ~0)
                  balance_residual=(e_in - du - absorbed) / max(abs(e_in), abs(du), abs(absorbed),
                                                                 float(ts["em_energy"][0]), 1e-30),
                  # raw peak over one RF period (a Hilbert envelope has end artefacts)
                  focus_peak_start_V_m=float(np.abs(ey_f[:per]).max()),
                  focus_peak_min_V_m=float(min(np.abs(ey_f[i:i + per]).max() for i in range(0, nq - per + 1, per))),
                  focus_peak_end_V_m=float(np.abs(ey_f[nq - per:nq]).max()))
        np.savez_compressed(out / "rf_series.npz", t_phasor=t_c, fwd=fwd.astype(np.complex64),
                            bwd=bwd.astype(np.complex64), p_out=p_out.astype(np.float32),
                            ey_focus=ey_f[:nq].astype(np.float32), dt=dt)
    summary["rf"] = rf
    summary["laser"] = dict(delivered_J=laser_deliv, absorbed_J=laser_abs, exit_J=laser_exit,
                            fraction_absorbed=laser_abs / laser_deliv if laser_deliv > 0 else 0.0)
    if imp_seed is not None:
        summary["laser"]["seed_electrons"] = float(imp_seed[0].sum() * dx)
        summary["laser"]["seed_electrons_focus_edge"] = float(imp_seed[0][e_focus].sum() * dx)
    bl = bulk_led.sum(axis=0) * vol
    fl = fil_led.sum(axis=0) * dx
    names = ["rf_in", "laser_in", "excitation", "elastic", "ionization", "recomb_carry", "gas_heat", "sink", "chem"]
    summary["ledger_bulk_J"] = dict(zip(names, map(float, bl)))
    summary["ledger_filament_J"] = dict(zip(names + ["hydro_ke", "p0_work", "laser_absorbed"], map(float, fl)))
    # closure: inputs = change of (electron + gas internal + kinetic) energy + sinks + chem + p0 work
    fin = float(fl[0] + fl[1])
    fres = fin - (fil_energy() - e_fil0 + fl[7] + fl[8] + fl[10])
    # floor: double-precision roundoff of the filament thermal energy
    summary["ledger_filament_J"]["closure_residual_rel"] = fres / max(fin, abs(fl[7]), 1e-12 * e_fil0)
    i_pk = int(np.argmax(ts["fil_ne_axis_focus"][:nb_done]))
    summary["peaks"] = dict(fil_ne_axis_max=float(ts["fil_ne_axis_focus"][i_pk]), t_fil_ne_max=float(ts["t"][i_pk]),
                            fil_te_axis_max=float(ts["fil_te_axis_focus"][:nb_done].max()),
                            fil_tg_axis_max=float(ts["fil_tg_axis_focus"][:nb_done].max()),
                            fil_electrons_end=float(ts["fil_electrons"][nb_done - 1]),
                            bulk_ne_max_end=float(ts["bulk_ne_max"][nb_done - 1]),
                            bulk_electrons_end=float(ts["bulk_electrons"][nb_done - 1]))
    if xt:
        summary["physics"] = dict(cr=case.cr, transport=case.transport, gas_heating=case.gas_heating,
                                  propagation_refraction=case.propagation, propagation_field=True,
                                  optics=gm.optics)
        summary["laser"]["lateral_J"] = laser_lat
        summary["laser"]["closure_residual_rel"] = ((laser_deliv - laser_abs - laser_exit - laser_lat)
                                                    / max(laser_deliv, 1e-30))
        b_in = float(bl[0] + bl[1])
        dw_b = float(wb.sum()) * QE * vol - bulk_e0
        summary["ledger_bulk_J"]["closure_residual_rel"] = ((b_in - dw_b - bl[6] - bl[7] - bl[8])
                                                            / max(abs(b_in), abs(bl[7]), 1e-30))
        heat_tg = float(0.5 * gm.prm[k3.P_DOF] * KB * n0 * (tg_b[active].sum() - heat_b0)) * vol
        summary["ledger_bulk_J"]["gas_heat_from_tg_J"] = heat_tg
        summary["peaks"].update({f"{k_}_max": float(v[:nb_done].max()) for k_, v in tx.items()})
        # isochoric bulk check: sound travel over the run vs cell size
        c_s = math.sqrt((1 + 2 / gm.prm[k3.P_DOF]) * KB * summary["peaks"]["bulk_tg_max_max"] / gm.prm[k3.P_MASS])
        summary["peaks"]["bulk_acoustic_length_over_dx"] = c_s * ts["t"][nb_done - 1] / dx
    np.savez_compressed(out / "timeseries.npz", **{k_: v[:nb_done] for k_, v in ts.items()},
                        **{k_: v[:nb_done] for k_, v in sp.items()}, bulk_ne_line_y=bulk_line[:nb_done],
                        **{k_: v[:nb_done] for k_, v in tx.items()}, **{k_: v[:nb_done] for k_, v in spx.items()},
                        fil_y=y_e, fil_js=js, beam_w=beam_w)
    np.savez_compressed(out / "snapshots.npz", **snaps, dx=dx, focus=np.array([fi, fj, fk]),
                        final_ne=(nxp + nx2p + 2 * n_alpha).astype(np.float32),
                        gas_k=np.array([lay.gas_k_lo, lay.gas_k_hi]))
    if frames["t"]:
        np.savez_compressed(out / "frames.npz", **{k_: np.asarray(v) for k_, v in frames.items()},
                            dx=dx, focus=np.array([fi, fj, fk]), fil_y=y_e, beam_w=beam_w,
                            gas_k=np.array([lay.gas_k_lo, lay.gas_k_hi]),
                            node_inside_xz=lay.node_inside[:, fj, :], meas_k=lay.meas_k,
                            trans_k=-1 if lay.trans_k is None else lay.trans_k)
    summary["wall_s"] = time.time() - t_wall
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    print(f"done {case.name}: {status}, {summary['wall_s']:.0f} s", flush=True)
    return summary


# (FWHM, energy at focus): reference 1 mJ is below the 300 ps breakdown threshold of all three gases
# in this model; 30 mJ / 300 ps and 100 mJ / 3 ns are above it (laser-only scan, taper_*laser-only*).
MATRIX_PULSES = ((300e-12, 1e-3), (300e-12, 30e-3), (3e-9, 100e-3))


def matrix_cases(geometries=("cavity", "taper"), gases=("Ar", "He", "N2")):
    cases = []
    for geo in geometries:
        for gas in gases:
            for fw, en in MATRIX_PULSES:
                cases.append(default_case(geo, gas, fw, energy_j=en, tag=f"E{en * 1e3:g}mJ"))
            cases.append(default_case(geo, gas, 3e-9, laser=False))           # RF-only control (16 ns)
    return cases


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m rf_plasma.coupled3d")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("calibrate"); c.add_argument("--apertures", type=float, nargs="+", default=[20, 25, 30, 35])
    c.add_argument("--q-target", type=float, default=400.0)
    p = sub.add_parser("prefill"); p.add_argument("--geometry", choices=["cavity", "taper", "guide"], required=True)
    r = sub.add_parser("run")
    r.add_argument("--geometry", choices=["cavity", "taper", "guide"], required=True)
    r.add_argument("--gas", choices=["Ar", "He", "N2"], required=True)
    r.add_argument("--fwhm", type=float, default=300e-12)
    r.add_argument("--no-laser", action="store_true"); r.add_argument("--no-rf", action="store_true")
    r.add_argument("--field", type=float, default=None); r.add_argument("--t-end", type=float, default=None)
    r.add_argument("--seed", type=float, default=None); r.add_argument("--tag", default="")
    r.add_argument("--energy", type=float, default=None, help="laser energy at focus (J)")
    r.add_argument("--waist", type=float, default=None, help="1/e^2 waist radius (m)")
    r.add_argument("--ee-scale", type=float, default=None,
                   help="multiplier on nu_ee in the EEDF closure criterion (default 1)")
    r.add_argument("--block", type=int, default=None, help="FDTD steps per kinetics/coupling block")
    r.add_argument("--rings", type=int, default=None, help="filament Lagrangian rings")
    m = sub.add_parser("matrix"); m.add_argument("--geometry", nargs="+", default=["cavity", "taper"])
    m.add_argument("--gas", nargs="+", default=["Ar", "He", "N2"]); m.add_argument("--skip-existing", action="store_true")
    for s in (c, p, r, m):
        s.add_argument("--dx", type=float, default=1e-3); s.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    if a.cmd == "calibrate":
        numba.set_num_threads(a.threads)
        print(calibrate(a.dx, a.apertures, a.q_target))
    elif a.cmd == "prefill":
        prefill(a.geometry, a.dx, a.threads)
    elif a.cmd == "run":
        kw = dict(dx_m=a.dx, tag=a.tag)
        if a.field: kw["field_target_v_m"] = a.field
        if a.seed: kw["seed_m3"] = a.seed
        if a.energy: kw["energy_j"] = a.energy
        if a.waist: kw["waist_m"] = a.waist
        if a.ee_scale is not None: kw["ee_scale"] = a.ee_scale
        if a.block: kw["block"] = a.block
        if a.rings: kw["rings"] = a.rings
        cs = default_case(a.geometry, a.gas, a.fwhm, laser=not a.no_laser, rf=not a.no_rf, **kw)
        if a.t_end: cs.t_end_s = a.t_end
        run_case(cs, a.threads)
    elif a.cmd == "matrix":
        for cs in matrix_cases(a.geometry, a.gas):
            cs.dx_m = a.dx
            if a.skip_existing and (RESULTS / cs.name / "summary.json").exists():
                continue
            run_case(cs, a.threads)


if __name__ == "__main__":
    main()
