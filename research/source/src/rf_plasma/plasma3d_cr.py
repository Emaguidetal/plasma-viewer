"""Extended 3D kinetics: excited-state CR chemistry, bulk transport + gas heating, laser propagation.

Stage 2 (collisional-radiative helium, lumped):
  state = [nX+, nX2+, W, nHe*, nHe2*, nHe2+(alpha)]; ne = nX+ + nX2+ + 2 nalpha (X2+ is the dimer ion).
  He* is the n=2 manifold (2^3S, 2^1S, 2^3P, 2^1P) with statistical internal populations, energy
  E* = sum g_i E_i / 16. Processes: excitation (electron excitation energy is stored as He*),
  superelastic e + He* -> e + He (Klein-Rosseland from the LXCat cross sections), stepwise ionization
  (Thomson), He* pooling (Penning), He* + 2He -> He2* + He, He* + He2* and He2* + He2* pooling, electron
  quenching and ionization of He2*, trapped resonance radiation (Holstein escape factor on the 2^1P
  share), He2* radiation, He+ + e -> He2+(alpha) + 2e (Voronov) and its recombination.
  Energy bookkeeping is exact per substep: chem = change of the stored potential energy
  U = Ip nX+ + (Ip - D) nX2+ + (Ip + I2) nalpha + E* nHe* + Ex nHe2*.
Stage 1: bulk_transport (effective free/ambipolar diffusion, electron energy convection + conduction,
  gas heat conduction; conservative, zero-flux at inactive nodes) and isochoric bulk gas heating.
Stage 3: filament_update_cr marches a radially resolved paraxial field A(r) (|A|^2 = I) through the
  filament edges: Crank-Nicolson diffraction, refraction by plasma (-ne / 2nc), neutral-density
  depletion and Kerr (n2 I), absorption removed ring by ring with exact energy accounting.
"""

import math

import numpy as np
from numba import njit, prange

from .plasma3d import (C0, EPS0, EPS_LO, EV_K, KB, ME, NE_FLOOR, N_LOSCHMIDT, QE, P_B1_0, P_B1_EXP, P_B1N_0,
                       P_B1N_EXP, P_B2_0, P_B2_C, P_B2_EXP, P_B2_THETA, P_B2N_0, P_B2N_EXP, P_DISS_A,
                       P_DISS_T1, P_DOF, P_EDEX, P_EDISS, P_EE_SCALE, P_ETA, P_IP, P_K3B0, P_K3B_EXP,
                       P_KCONV0, P_KCONV_EXP, P_KDEX0, P_KDISS0, P_MASS, P_MU_X2P, P_MU_XP, P_NPHOT,
                       P_OMEGA_L, P_PHOTON_EV, _eloss, _fil_hydro, _lookup, _rate, _thomas,
                       _wmax, drude_params, nu_ei)
from . import plasma3d as _k3

FIL_CNT_MIN = 1e-9  # charges per metre of filament ring below which a ring is treated as empty

# CR parameters appended to the base parameter vector
P_CR_ON = _k3.N_PARAMS
P_ESTAR, P_EX, P_DION, P_I2 = P_CR_ON + 1, P_CR_ON + 2, P_CR_ON + 3, P_CR_ON + 4
P_KC3 = P_CR_ON + 5                       # He* + 2He -> He2* + He (m^6/s)
P_KMM_A, P_KMM_M = P_CR_ON + 6, P_CR_ON + 7   # He* + He* -> He+ / He2+ + e (m^3/s)
P_KMX_A, P_KMX_M = P_CR_ON + 8, P_CR_ON + 9   # He* + He2* -> He+ / He2+ + e
P_KXX = P_CR_ON + 10                      # He2* + He2* -> He2+ + 2He + e
P_KXQ = P_CR_ON + 11                      # e + He2* -> e + 2He
P_KXI0, P_KXIE = P_CR_ON + 12, P_CR_ON + 13   # e + He2* -> He2+ + 2e : k0 exp(-E/Te)
P_ARES, P_FRES = P_CR_ON + 14, P_CR_ON + 15   # trapped resonance rate (s^-1), 2^1P share of He*
P_AX = P_CR_ON + 16                       # He2* radiative rate (s^-1)
P_KAREC0, P_KAREC_EXP = P_CR_ON + 17, P_CR_ON + 18  # e + alpha + He -> He+ (m^6/s) (Te/300K)^x
# EEDF closure with e-e collisions solved (swarm_ee): 4-row Boltzmann blocks per ionization degree
# x_j = 0, 1e-7, 1e-6, ... appended to tabs at row P_EE_ROW0; P_EE_NLEV = 0 keeps the w-blend
P_EE_ROW0, P_EE_NLEV = P_CR_ON + 19, P_CR_ON + 20
# filament radial transport: 0 = ambipolar everywhere (original); 1 = Allis-Rose free -> ambipolar
# transition on the local density scale length (lets the low-density leading edge diffuse freely)
P_FIL_FREE = P_CR_ON + 21
N_PARAMS_CR = P_CR_ON + 22
ROW_KSUP, ROW_KSTEP = 8, 9                # extra table rows (Maxwellian, per He* atom)
EE_X1 = 1e-7                              # first nonzero ionization-degree level (decades above)


@njit(cache=True, inline="always")
def _wsel(ne, n_neu, eps, tabs, ln_e0, h, prm):
    """Blend weight w >= 0 (base closure), or -(1 + fractional ionization-degree level) for e-e tables."""
    nlev = int(prm[P_EE_NLEV])
    if nlev < 2:
        return _wmax(ne, n_neu, eps, tabs, ln_e0, h, prm[P_IP], prm[P_EE_SCALE])
    x = ne / max(n_neu, 1.0)
    if x <= EE_X1:
        f = max(x, 0.0) / EE_X1
    else:
        f = min(1.0 + math.log10(x / EE_X1), nlev - 1.0)
    return -1.0 - f


@njit(cache=True, inline="always")
def _rx(tabs, row, w, ln_e0, h, eps, prm):
    if w >= 0.0:
        return _rate(tabs, row, w, ln_e0, h, eps)
    f = -w - 1.0
    j = int(f)
    r0 = int(prm[P_EE_ROW0])
    v = _lookup(tabs[r0 + 4 * j + row], ln_e0, h, eps)
    fr = f - j
    if fr > 0.0 and j + 1 < int(prm[P_EE_NLEV]):
        v = (1.0 - fr) * v + fr * _lookup(tabs[r0 + 4 * (j + 1) + row], ln_e0, h, eps)
    return v


@njit(cache=True, inline="always")
def _elx(e, ne, n_gas, tabs, ln_e0, h, eps_g, two_m_over_m, w, prm):
    """Electron energy loss per electron (eV/s), as plasma3d._eloss with the selected closure."""
    exc = n_gas * _rx(tabs, 2, w, ln_e0, h, e, prm)
    el = (n_gas * _rx(tabs, 3, w, ln_e0, h, e, prm) + two_m_over_m * nu_ei(ne, e / 1.5)) * (e - eps_g)
    return exc + el


def extend_params(prm: np.ndarray) -> np.ndarray:
    out = np.zeros(N_PARAMS_CR)
    out[: len(prm)] = prm
    return out


@njit(cache=True, inline="always")
def _voronov_he1(te):
    """He+ -> He2+ ionization (Voronov 1997 fit, m^3/s); exponentially small below its 3 eV range."""
    u = 54.4 / max(te, 0.05)
    if u > 600.0:
        return 0.0
    return 2.05e-15 * (1.0 + math.sqrt(u)) / (0.265 + u) * u ** 0.25 * math.exp(-u)


@njit(cache=True)
def kinetic_cell_cr(state, dt, n_gas, tg, p_ext, intensity, prm, tabs, ln_e0, h,
                    ln_i0, hi, ppt_tab, out, path_len=0.0, max_substeps=200000):
    """Advance one cell; state (6,) = [nXp, nX2p, W, nHe*, nHe2*, nalpha]; out as kinetic_cell
    [rf, laser, excitation(net of superelastic), elastic, ionization(info), recomb_carry, gas, sink, chem].
    n_gas is ground-state neutral density, not total nuclei. Optional out[9:11]
    separates resonance/excimer radiative energy from the total sink.
    Raise on substep exhaustion rather than returning a partially advanced state.
    With prm[P_CR_ON] = 0 and prm[P_DION] = 0 this reproduces plasma3d.kinetic_cell."""
    km_t = tabs[0]
    eps_max = math.exp(ln_e0 + h * (km_t.shape[0] - 1))
    ip = prm[P_IP]
    eta = prm[P_ETA]
    cr = prm[P_CR_ON] > 0.0 and tabs.shape[0] > ROW_KSTEP
    estar = prm[P_ESTAR]; ex = prm[P_EX]; dion = prm[P_DION]; i2 = prm[P_I2]
    ipd = ip - dion
    eps_g = 1.5 * tg / EV_K
    two_m_over_m = 2.0 * ME / prm[P_MASS]
    omega_l = prm[P_OMEGA_L]
    w_ppt = 0.0
    if intensity > 0.0:
        x = (math.log(intensity) - ln_i0) / hi
        if x >= 0.0:
            n = ppt_tab.shape[0]
            if x >= n - 1:
                w_ppt = math.exp(ppt_tab[n - 1])
            else:
                i = int(x)
                f = x - i
                w_ppt = math.exp(ppt_tab[i] * (1 - f) + ppt_tab[i + 1] * f)
    nphot = prm[P_NPHOT]
    excess = max(nphot * prm[P_PHOTON_EV] - ip, 0.0)
    kconv_c = prm[P_KCONV0] * (tg / 300.0) ** prm[P_KCONV_EXP]
    tg_ev = tg / EV_K
    kdiss_c = 0.0
    if prm[P_KDISS0] > 0.0:
        kdiss_c = prm[P_KDISS0] * tg_ev ** prm[P_DISS_A] * math.exp(-prm[P_EDISS] / tg_ev)
        if prm[P_DISS_T1] > 0.0:
            kdiss_c *= math.exp(min(tg / prm[P_DISS_T1], 50.0))
    theta = prm[P_B2_THETA]
    b2_tg = 1.0
    if theta > 0.0:
        exn = math.exp(-theta / tg)
        b2_tg = (1.0 - exn) / (1.0 - prm[P_B2_C] * exn)
    t = 0.0
    nxp = state[0]; nx2p = state[1]; w = state[2]
    nm = state[3]; nx = state[4]; na = state[5]
    n_tot = n_gas + nxp + 2.0 * nx2p + na + nm + 2.0 * nx
    iters = 0
    while t < dt:
        if iters >= max_substeps:
            raise ArithmeticError("CR kinetics exhausted its substep budget before completing dt")
        iters += 1
        ne = nxp + nx2p + 2.0 * na
        ne_true = ne
        if ne < NE_FLOOR:
            ne = NE_FLOOR
        eps = w / ne
        if eps < EPS_LO:
            eps = EPS_LO
        if eps > eps_max:
            eps = eps_max
        te = eps / 1.5
        te_k = te * EV_K
        n_neu = max(n_tot - nxp - 2.0 * nx2p - na - nm - 2.0 * nx, 0.0)
        kconv = kconv_c * n_neu * n_neu
        kdiss = kdiss_c * n_neu
        wmx = _wsel(ne, n_neu, eps, tabs, ln_e0, h, prm)
        km = _rx(tabs, 0, wmx, ln_e0, h, eps, prm)
        kiz = _rx(tabs, 1, wmx, ln_e0, h, eps, prm)
        nu = n_neu * km + nu_ei(ne, te)
        r_iz = kiz * n_neu * ne_true  # the floor regularizes eps only; it never seeds an avalanche
        r_ppt = w_ppt * n_neu
        b1 = prm[P_B1_0] * (te_k / 300.0) ** prm[P_B1_EXP] if prm[P_B1_0] > 0 else 0.0
        b1 += prm[P_B1N_0] * (te_k / 300.0) ** prm[P_B1N_EXP] * n_neu if prm[P_B1N_0] > 0 else 0.0
        b2 = prm[P_B2_0] * (te_k / 300.0) ** prm[P_B2_EXP] * b2_tg
        b2 += prm[P_B2N_0] * (te_k / 300.0) ** prm[P_B2N_EXP] * n_neu if prm[P_B2N_0] > 0 else 0.0
        k3 = prm[P_K3B0] * te ** prm[P_K3B_EXP]
        kdex = 0.0
        if prm[P_KDEX0] > 0.0:
            kdex = prm[P_KDEX0] * math.exp(-(prm[P_EDEX] - 3.0 * (tg_ev - 0.026)) / te)
        l1 = b1 * ne + k3 * ne * ne
        l2 = b2 * ne + k3 * ne * ne
        d1 = kdiss + kdex * ne
        # ---- CR rates (zero when CR is off)
        ksup = kst = kxi = kxq = kaz = conv3 = arad = 0.0
        kmm = kmx = kxx = karec = ax = 0.0
        if cr:
            ksup = _lookup(tabs[ROW_KSUP], ln_e0, h, eps)
            kst = _lookup(tabs[ROW_KSTEP], ln_e0, h, eps)
            kxi = prm[P_KXI0] * math.exp(-prm[P_KXIE] / te)
            kxq = prm[P_KXQ]
            kaz = _voronov_he1(te) * ne
            conv3 = prm[P_KC3] * n_neu * n_neu
            arad = prm[P_ARES] * prm[P_FRES]
            kmm = prm[P_KMM_A] + prm[P_KMM_M]
            kmx = prm[P_KMX_A] + prm[P_KMX_M]
            kxx = prm[P_KXX]
            ax = prm[P_AX]
            karec = prm[P_KAREC0] * (te_k / 300.0) ** prm[P_KAREC_EXP] * n_neu * ne
        p_ib = 0.0
        if intensity > 0.0:
            p_ib = QE * nu * intensity / (ME * EPS0 * C0 * omega_l * omega_l)
            if path_len > 0.0:
                a_l = (ne * p_ib + r_ppt * nphot * prm[P_PHOTON_EV]) * QE / intensity * path_len
                if a_l > 1e-6:
                    fa = (1.0 - math.exp(-a_l)) / a_l
                    p_ib *= fa
                    r_ppt *= fa
        gain = (estar * ksup * ne * nm + ex * kxq * ne * nx
                + max(2 * estar - ip, 0.0) * kmm * nm * nm + max(2 * ex - ipd, 0.0) * kxx * nx * nx)
        cost = (ip - estar) * kst * ne * nm + i2 * kaz * nxp + max(ipd - ex, 0.0) * kxi * ne * nx
        heat_per_e = (p_ext + gain) / ne
        rate = max(kiz * n_neu, kiz * ne, w_ppt, l1, l2, 1e3)
        relax = _elx(eps, ne, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx, prm) / eps
        dts = min(dt - t, 0.05 / rate)
        if r_iz > 0.0:
            dts = min(dts, 0.3 * max(w, ne * EPS_LO) / (ip * r_iz))
        if cost > 0.0:
            dts = min(dts, 0.3 * max(w, ne * EPS_LO) / cost)
        if heat_per_e + p_ib > 0.0:
            dts = min(dts, max(0.3 * eps / (heat_per_e + p_ib), 0.05 / max(relax, 1.0)))
        if dts < 1e-17:
            dts = 1e-17
        # ---- He* (implicit loss; pooling lagged in nm)
        lm = (ksup + kst) * ne + conv3 + arad + 2.0 * kmm * nm + kmx * nx
        nm1 = nm / (1.0 + dts * lm)
        fm = dts * nm1
        f_sup = fm * ksup * ne; f_st = fm * kst * ne; f_cv = fm * conv3; f_rad = fm * arad
        p_mma = fm * prm[P_KMM_A] * nm if cr else 0.0
        p_mmm = fm * prm[P_KMM_M] * nm if cr else 0.0
        p_mxa = fm * prm[P_KMX_A] * nx if cr else 0.0
        p_mxm = fm * prm[P_KMX_M] * nx if cr else 0.0
        pmx = p_mxa + p_mxm
        if pmx > 0.5 * nx and pmx > 0.0:      # cannot consume more excimers than present
            s = 0.5 * nx / pmx
            nm1 += (1.0 - s) * pmx
            p_mxa *= s; p_mxm *= s; pmx *= s
        # ---- He2* (implicit, source from conversion)
        nxa = nx - pmx
        lx = (kxi + kxq) * ne + ax + 2.0 * kxx * nx
        nx1 = (nxa + f_cv) / (1.0 + dts * lx)
        fx = dts * nx1
        x_ion = fx * kxi * ne; x_q = fx * kxq * ne; x_rad = fx * ax; p_xx = fx * kxx * nx
        # ---- alpha (He2+ bare nucleus)
        f_az = dts * kaz * nxp / (1.0 + dts * kaz)
        na1 = (na + f_az) / (1.0 + dts * karec)
        r_a = dts * karec * na1
        nxp_pre = nxp - f_az + r_a
        s1 = f_st + p_mma + p_mxa
        s2 = x_ion + p_mmm + p_mxm + p_xx
        # ---- ions (implicit 2x2, charge conserving)
        a11 = 1.0 + dts * (kconv + l1)
        a12 = -dts * d1
        a21 = -dts * kconv
        a22 = 1.0 + dts * (d1 + l2)
        r1 = nxp_pre + dts * (r_iz + r_ppt) + s1
        r2 = nx2p + s2
        det = a11 * a22 - a12 * a21
        nxp_n = (r1 * a22 - a12 * r2) / det
        nx2p_n = (a11 * r2 - a21 * r1) / det
        rec1 = dts * l1 * nxp_n
        rec2 = dts * l2 * nx2p_n
        conv_net = dts * (kconv * nxp_n - d1 * nx2p_n)
        ne_n = nxp_n + nx2p_n + 2.0 * na1
        r_rec = (rec1 + rec2 + r_a) / dts
        q_e = (estar * f_sup - (ip - estar) * f_st + (2 * estar - ip) * p_mma + (2 * estar - ipd) * p_mmm
               + (estar + ex - ip) * p_mxa + (estar + ex - ipd) * p_mxm + (2 * ex - ipd) * p_xx
               - (ipd - ex) * x_ion + ex * x_q - i2 * f_az)
        src = p_ext + ne * p_ib + excess * r_ppt - ip * r_iz
        rhs = w + dts * src + q_e
        lo = EPS_LO
        hi_e = eps_max
        e1 = min(max(eps, lo), hi_e)
        ok = False
        for _ in range(8):
            g = ne_n * e1 + dts * (ne_n * _elx(e1, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx, prm) + e1 * r_rec) - rhs
            de = e1 * 1e-4 + 1e-12
            e2 = e1 + de
            g2 = ne_n * e2 + dts * (ne_n * _elx(e2, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx, prm) + e2 * r_rec) - rhs
            dg = (g2 - g) / de
            if dg <= 0.0:
                break
            e_new = e1 - g / dg
            if e_new < lo:
                e_new = 0.5 * (e1 + lo)
            if e_new > hi_e:
                e_new = 0.5 * (e1 + hi_e)
            if abs(e_new - e1) < 1e-7 * e1:
                e1 = e_new
                ok = True
                break
            e1 = e_new
        if not ok:
            a, b = lo, hi_e
            for _ in range(60):
                mid = math.sqrt(a * b)
                g = ne_n * mid + dts * (ne_n * _elx(mid, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx, prm) + mid * r_rec) - rhs
                if g > 0:
                    b = mid
                else:
                    a = mid
            e1 = math.sqrt(a * b)
        e_exc = dts * ne_n * n_neu * _rx(tabs, 2, wmx, ln_e0, h, e1, prm)
        e_el = dts * ne_n * (_elx(e1, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx, prm)
                             - n_neu * _rx(tabs, 2, wmx, ln_e0, h, e1, prm))
        e_rc = dts * e1 * r_rec
        w_new = rhs - e_exc - e_el - e_rc
        rel = ip * rec1 + ipd * rec2 + i2 * r_a
        u_old = ip * nxp + ipd * nx2p + (ip + i2) * na + estar * nm + ex * nx
        if cr:
            nm1 += e_exc / estar
            gas = e_el + dion * conv_net + (estar - ex) * f_cv + eta * (rel + e_rc)
            sink = estar * f_rad + ex * x_rad + (1.0 - eta) * (rel + e_rc)
        else:
            gas = e_el + dion * conv_net + eta * (e_exc + rel + e_rc)
            sink = (1.0 - eta) * (e_exc + rel + e_rc)
        nxp = nxp_n; nx2p = nx2p_n; w = w_new; nm = nm1; nx = nx1; na = na1
        u_new = ip * nxp + ipd * nx2p + (ip + i2) * na + estar * nm + ex * nx
        out[0] += dts * p_ext * QE
        out[1] += dts * (ne * p_ib + r_ppt * nphot * prm[P_PHOTON_EV]) * QE
        out[2] += (e_exc - estar * f_sup - ex * x_q) * QE
        out[3] += e_el * QE
        out[4] += (dts * ip * r_iz - q_e + estar * f_sup + ex * x_q) * QE
        out[5] += e_rc * QE
        out[6] += gas * QE
        out[7] += (sink + dts * r_ppt * (nphot * prm[P_PHOTON_EV] - ip - excess)) * QE
        out[8] += (u_new - u_old) * QE
        if out.shape[0] > 9:
            out[9] += estar * f_rad * QE
        if out.shape[0] > 10:
            out[10] += ex * x_rad * QE
        t += dts
    state[0] = nxp; state[1] = nx2p; state[2] = w
    state[3] = nm; state[4] = nx; state[5] = na
    return iters


# ---------------------------------------------------------------- bulk (3D node grid)
@njit(parallel=True, cache=True)
def bulk_update_cr(st6, tg, active, node_work, dt, n_gas, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab,
                   wp2_out, nu_out, ledger, heat_gas):
    """st6: (6, nx, ny, nz) state; tg: node gas temperature (K), updated isochorically when heat_gas."""
    _, nxp1, nyp1, nzp1 = st6.shape
    dof = prm[P_DOF]
    for i in prange(nxp1):
        st = np.empty(6)
        cell = np.zeros(9)
        for j in range(nyp1):
            for k in range(nzp1):
                if not active[i, j, k]:
                    wp2_out[i, j, k] = 0.0
                    nu_out[i, j, k] = 0.0
                    continue
                for q in range(6):
                    st[q] = st6[q, i, j, k]
                p_ext = node_work[i, j, k] / (QE * dt)
                n_neu = max(n_gas - st[0] - 2.0 * st[1] - st[5] - st[3] - 2.0 * st[4], 1e10)
                g0 = cell[6]
                kinetic_cell_cr(st, dt, n_neu, tg[i, j, k], p_ext, 0.0, prm, tabs, ln_e0, h,
                                ln_i0, hi, ppt_tab, cell)
                if heat_gas:
                    npart = max(n_gas - st[1] - st[4], 1e10)
                    tg[i, j, k] = max(tg[i, j, k] + (cell[6] - g0) / (0.5 * dof * KB * npart), 1.0)
                for q in range(6):
                    st6[q, i, j, k] = st[q]
                ne = st[0] + st[1] + 2.0 * st[5]
                n_neu = max(n_gas - st[0] - 2.0 * st[1] - st[5] - st[3] - 2.0 * st[4], 1e10)
                wp2, nu = drude_params(ne, st[2], n_neu, tabs, ln_e0, h, prm)
                wp2_out[i, j, k] = wp2
                nu_out[i, j, k] = nu
        for q in range(9):
            ledger[i, q] += cell[q]


@njit(cache=True, inline="always")
def _node_coef(st6, tg, i, j, k, n_gas, prm, tabs, ln_e0, h, dx2):
    """(D_eff X+, D_eff X2+, D_e, electron conductivity factor, gas conductivity, eps)."""
    ne = st6[0, i, j, k] + st6[1, i, j, k] + 2.0 * st6[5, i, j, k]
    eps = st6[2, i, j, k] / ne if ne > 0.0 else EPS_LO
    eps = max(eps, EPS_LO)
    te = eps / 1.5
    wmx = _wsel(ne, n_gas, eps, tabs, ln_e0, h, prm)
    nu = n_gas * _rx(tabs, 0, wmx, ln_e0, h, eps, prm) + nu_ei(ne, te)
    d_e = QE * te / (ME * nu) if nu > 0.0 else 0.0
    ti = tg[i, j, k] / EV_K
    da1 = prm[P_MU_XP] * N_LOSCHMIDT / n_gas * (te + ti)
    da2 = prm[P_MU_X2P] * N_LOSCHMIDT / n_gas * (te + ti)
    # free -> ambipolar transition (Allis & Rose 1954 limits), interpolated in (L / lambda_D)^2, L = dx
    x = dx2 * ne * QE / (EPS0 * te) if ne > 0.0 else 0.0
    d1 = (d_e + x * da1) / (1.0 + x)
    d2 = (d_e + x * da2) / (1.0 + x)
    kap = 2.5 * ne * d_e  # eV-flux per (eV/m): q = -kap grad Te
    kg = 0.155 * (tg[i, j, k] / 300.0) ** 0.7  # He gas conductivity, W/m/K
    return d1, d2, d_e, kap, kg, eps


@njit(parallel=True, cache=True)
def bulk_transport(st6, tg, active, dt, dx, n_gas, prm, tabs, ln_e0, h, coef, delta, gas_cond):
    """One explicit conservative step (caller sub-cycles for stability). coef: (6, nx, ny, nz) scratch,
    delta: (5, nx, ny, nz) scratch for [nXp, nX2p, W, nalpha, gas energy]. Returns max diffusivity."""
    _, nx_, ny_, nz_ = st6.shape
    dx2 = dx * dx
    dmax = np.zeros(nx_)
    for i in prange(nx_):
        for j in range(ny_):
            for k in range(nz_):
                if active[i, j, k]:
                    c = _node_coef(st6, tg, i, j, k, n_gas, prm, tabs, ln_e0, h, dx2)
                    for q in range(6):
                        coef[q, i, j, k] = c[q]
                    dm = max(c[0], c[1], 5.0 / 3.0 * c[2])
                    if dm > dmax[i]:
                        dmax[i] = dm
    cv = 1.5 * KB * n_gas
    for i in prange(nx_):
        for j in range(ny_):
            for k in range(nz_):
                for q in range(5):
                    delta[q, i, j, k] = 0.0
                if not active[i, j, k]:
                    continue
                for nb in range(6):
                    ii = np.int64(i); jj = np.int64(j); kk = np.int64(k)
                    if nb == 0:
                        ii -= 1
                    elif nb == 1:
                        ii += 1
                    elif nb == 2:
                        jj -= 1
                    elif nb == 3:
                        jj += 1
                    elif nb == 4:
                        kk -= 1
                    else:
                        kk += 1
                    if ii < 0 or jj < 0 or kk < 0 or ii >= nx_ or jj >= ny_ or kk >= nz_:
                        continue
                    if not active[ii, jj, kk]:
                        continue
                    # inflow into (i,j,k) from the neighbour; antisymmetric by construction
                    d1f = 0.5 * (coef[0, i, j, k] + coef[0, ii, jj, kk])
                    d2f = 0.5 * (coef[1, i, j, k] + coef[1, ii, jj, kk])
                    g1 = d1f * (st6[0, ii, jj, kk] - st6[0, i, j, k]) / dx2
                    g2 = d2f * (st6[1, ii, jj, kk] - st6[1, i, j, k]) / dx2
                    ga = d1f * (st6[5, ii, jj, kk] - st6[5, i, j, k]) / dx2
                    ge = g1 + g2 + 2.0 * ga  # electron inflow (quasi-neutral)
                    eps_f = 0.5 * (coef[5, i, j, k] + coef[5, ii, jj, kk])
                    kap_f = 0.5 * (coef[3, i, j, k] + coef[3, ii, jj, kk])
                    te_n = coef[5, ii, jj, kk] / 1.5
                    te_c = coef[5, i, j, k] / 1.5
                    gw = 5.0 / 3.0 * eps_f * ge + kap_f * (te_n - te_c) / dx2
                    delta[0, i, j, k] += g1
                    delta[1, i, j, k] += g2
                    delta[3, i, j, k] += ga
                    delta[2, i, j, k] += gw
                    if gas_cond:
                        kgf = 0.5 * (coef[4, i, j, k] + coef[4, ii, jj, kk])
                        delta[4, i, j, k] += kgf * (tg[ii, jj, kk] - tg[i, j, k]) / dx2
    for i in prange(nx_):
        for j in range(ny_):
            for k in range(nz_):
                if active[i, j, k]:
                    st6[0, i, j, k] = max(st6[0, i, j, k] + dt * delta[0, i, j, k], 0.0)
                    st6[1, i, j, k] = max(st6[1, i, j, k] + dt * delta[1, i, j, k], 0.0)
                    st6[5, i, j, k] = max(st6[5, i, j, k] + dt * delta[3, i, j, k], 0.0)
                    st6[2, i, j, k] = max(st6[2, i, j, k] + dt * delta[2, i, j, k], 0.0)
                    if gas_cond:
                        tg[i, j, k] += dt * delta[4, i, j, k] / cv
    return dmax.max()


# ---------------------------------------------------------------- laser propagation (radial paraxial)
@njit(cache=True)
def _cn_radial(a, dr, dz, kk, work_a, work_b, work_c, work_d, work_x):
    """Crank-Nicolson step of dA/dz = (i / 2k) (1/r) d/dr (r dA/dr), finite volume on r_j = (j+1/2) dr,
    zero flux at r = 0 and at the outer edge: conserves sum |A_j|^2 r_j exactly (Hermitian operator)."""
    n = a.shape[0]
    c = 1j * dz / (4.0 * kk * dr * dr)
    for j in range(n):
        rj = (j + 0.5)
        rm = j * 1.0  # face r_{j-1/2} / dr
        rp = j + 1.0
        lo = rm / rj if j > 0 else 0.0
        up = rp / rj if j < n - 1 else 0.0
        diag = -(lo + up)
        lap = diag * a[j]
        if j > 0:
            lap += lo * a[j - 1]
        if j < n - 1:
            lap += up * a[j + 1]
        work_d[j] = a[j] + c * lap
        work_b[j] = 1.0 - c * diag
        work_a[j] = -c * lo
        work_c[j] = -c * up
    # complex Thomas
    cp = np.empty(n, np.complex128)
    dp = np.empty(n, np.complex128)
    cp[0] = work_c[0] / work_b[0]
    dp[0] = work_d[0] / work_b[0]
    for j in range(1, n):
        m = work_b[j] - work_a[j] * cp[j - 1]
        cp[j] = work_c[j] / m
        dp[j] = (work_d[j] - work_a[j] * dp[j - 1]) / m
    a[n - 1] = dp[n - 1]
    for j in range(n - 2, -1, -1):
        a[j] = dp[j] - cp[j] * a[j + 1]


@njit(cache=True)
def _overlap(rb_e, rg_edges, ov):
    """ov[p, z] = area of grid annulus p inside ring z (exact partition)."""
    np_ = rg_edges.shape[0] - 1
    nz = rb_e.shape[0] - 1
    for p in range(np_):
        for z in range(nz):
            ov[p, z] = 0.0
    z = 0
    for p in range(np_):
        a0 = rg_edges[p]; a1 = rg_edges[p + 1]
        while z < nz and rb_e[z + 1] <= a0:
            z += 1
        zz = z
        while zz < nz and rb_e[zz] < a1:
            lo = max(a0, rb_e[zz]); hi = min(a1, rb_e[zz + 1])
            if hi > lo:
                ov[p, zz] = math.pi * (hi * hi - lo * lo)
            zz += 1


@njit(cache=True)
def _fil_diffuse_cr(e, rb, heavy, cs6, dt, prm, tabs, ln_e0, h):
    """Energy-conservative ambipolar diffusion + electron heat conduction for the CR filament.

    Unlike plasma3d._fil_diffuse, electrons (cxp + cx2p + 2 alpha) carry their upwind mean energy
    across each face (face fluxes rebuilt from the zone count changes), so sum(ue) is conserved
    to round-off. Alpha diffuses with the X+ mobility; He*/He2* are left in place (their neutral
    diffusion length over the run is < 1 um at 1 MPa)."""
    n_edges, nb = rb.shape
    nr = nb - 1
    cxp = cs6[0, e]; cx2p = cs6[1, e]; ue = cs6[2, e]; alp = cs6[5, e]
    area = np.empty(nr)
    for z in range(nr):
        area[z] = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
    te_floor = EPS_LO / 1.5
    # near-empty rings (counts ~1e-40 /m) give round-off eps = ue/cnt up to 1e44 eV. Cap only the temperature
    # entering the transport coefficients; the energy fluxes use the exact eps so sum(ue) is conserved.
    te_cap = math.exp(ln_e0 + h * (tabs.shape[1] - 1)) / 1.5
    cnt0 = np.empty(nr); eps0 = np.empty(nr); te_r = np.empty(nr)
    for z in range(nr):
        # < 1e-9 charges per metre of ring is statistically empty; its ue/cnt is pure round-off
        if cxp[z] + cx2p[z] + 2.0 * alp[z] < FIL_CNT_MIN:
            cxp[z] = 0.0; cx2p[z] = 0.0; alp[z] = 0.0; ue[z] = 0.0
        cnt0[z] = cxp[z] + cx2p[z] + 2.0 * alp[z]
        eps0[z] = ue[z] / cnt0[z] if cnt0[z] > 0.0 else 0.0
        te_r[z] = min(max(eps0[z] / 1.5, te_floor), te_cap)
    dxp = np.zeros(nr + 1); dx2p = np.zeros(nr + 1); kap = np.zeros(nr + 1)
    free = prm[P_FIL_FREE] > 0.0
    r_out = rb[e, nr]
    for f in range(1, nr):
        nel = cnt0[f - 1] / area[f - 1]
        ner = cnt0[f] / area[f]
        te_f = 0.5 * (te_r[f - 1] + te_r[f])
        ng_f = 0.5 * (heavy[e, f - 1] / area[f - 1] + heavy[e, f] / area[f])
        dxp[f] = prm[P_MU_XP] * N_LOSCHMIDT / ng_f * te_f
        dx2p[f] = prm[P_MU_X2P] * N_LOSCHMIDT / ng_f * te_f
        nmin = min(nel, ner)
        wmx_f = _wsel(nmin, ng_f, 1.5 * te_f, tabs, ln_e0, h, prm)
        nu_f = ng_f * _rx(tabs, 0, wmx_f, ln_e0, h, 1.5 * te_f, prm) + nu_ei(nmin, te_f)
        if nu_f > 0.0 and nmin > 0.0:
            kap[f] = 2.5 * nmin * te_f / (ME * nu_f) * QE
        if free and nu_f > 0.0:
            # Allis & Rose (1954) limits interpolated in (L / lambda_D)^2 with L the local density scale
            # length, clamped to [face spacing, filament radius]
            d_e = QE * te_f / (ME * nu_f)
            dr_c = 0.5 * (rb[e, f + 1] - rb[e, f - 1])
            grad = abs(ner - nel) / dr_c
            nbar = 0.5 * (nel + ner)
            l_s = nbar / grad if grad > 0.0 else r_out
            l_s = min(max(l_s, dr_c), r_out)
            x = l_s * l_s * nmin * QE / (EPS0 * te_f)
            dxp[f] = (d_e + x * dxp[f]) / (1.0 + x)
            dx2p[f] = (d_e + x * dx2p[f]) / (1.0 + x)
    _k3._diffuse_counts(cxp, rb[e], area, dxp, dt)
    _k3._diffuse_counts(cx2p, rb[e], area, dx2p, dt)
    _k3._diffuse_counts(alp, rb[e], area, dxp, dt)
    cnt1 = np.empty(nr)
    for z in range(nr):
        cnt1[z] = cxp[z] + cx2p[z] + 2.0 * alp[z]
    # face electron flux F_f (zone f-1 -> f) = -(sum of count changes inside face f); upwind energy
    flux = 0.0
    for f in range(1, nr):
        flux -= cnt1[f - 1] - cnt0[f - 1]
        # a face cannot carry more than the donor ring holds (the running sum keeps round-off from dense rings)
        if flux > 0.0:
            q = min(flux, cnt0[f - 1]) * eps0[f - 1]
        else:
            q = max(flux, -cnt0[f]) * eps0[f]
        ue[f - 1] -= q
        ue[f] += q
    for z in range(nr):
        if ue[z] < 0.0:
            ue[z] = 0.0
        te_r[z] = ue[z] / (1.5 * cnt1[z]) if cnt1[z] > 0.0 else 0.0
    _k3._conduct_te(te_r, cnt1, rb[e], kap, dt)
    for z in range(nr):
        ue[z] = 1.5 * te_r[z] * cnt1[z]


@njit(parallel=True, cache=True)
def filament_update_cr(rb, ub, heavy, eg, cs6, ring_work, order, laser_in, dt, n0, tg0, p0, prm, tabs,
                       ln_e0, h, ln_i0, hi, ppt_tab, field, dr, k_l, n2, refr, n_crit, ledger, hydro_on,
                       path_len, prop_on, lat_mask, beam_out):
    """Filament with CR kinetics and (if prop_on) radial paraxial propagation of `field` (complex,
    |A|^2 in W/m^2 on r_j = (j+1/2) dr), entering at edge order[0] and marched through every edge.
    cs6: (6, nE, R) counts per length of [X+, X2+, electron energy (eV/m), He*, He2*, alpha].
    refr: neutral refractivity (n-1) at density n0; n2: Kerr index (m^2/W) at n0; n_crit: critical density.
    lat_mask: per-substep amplitude mask (absorbing outer layer). beam_out[0] += laterally absorbed J,
    beam_out[1] = transmitted power after the last edge (W); if len >= 2 + 2 nE, beam_out[2 + e] = on-axis
    intensity and beam_out[2 + nE + e] = second-moment radius at the entry of edge e.
    ledger: (nE, 12) [9 kinetic channels, hydro_ke, p0_work, laser_absorbed]."""
    n_edges, nb = rb.shape
    nr = nb - 1
    mass = prm[P_MASS]
    dof = prm[P_DOF]
    gam = 1.0 + 2.0 / dof
    npg = field.shape[0]
    rg_edges = np.arange(npg + 1) * dr
    dA = np.empty(npg)
    for p in range(npg):
        dA[p] = math.pi * (rg_edges[p + 1] ** 2 - rg_edges[p] ** 2)
    cxp = cs6[0]; cx2p = cs6[1]; ue = cs6[2]
    for e in prange(n_edges):
        _fil_hydro(e, rb, ub, heavy, eg, cx2p, ue, dt, p0, mass, dof, gam, ledger, hydro_on)
    acc_r = np.zeros((nr, 9))
    ring_i = np.zeros(nr)
    frac = np.zeros(nr)
    ov = np.zeros((npg, nr))
    wa = np.empty(npg, np.complex128); wb = np.empty(npg, np.complex128)
    wc = np.empty(npg, np.complex128); wd = np.empty(npg, np.complex128); wx = np.empty(npg, np.complex128)
    dn = np.zeros(npg)
    for ii in range(n_edges):
        e = order[ii]
        _overlap(rb[e], rg_edges, ov)
        if prop_on and beam_out.shape[0] >= 2 + 2 * n_edges:
            p_e = 0.0
            m2 = 0.0
            for p in range(npg):
                ip_ = (field[p].real ** 2 + field[p].imag ** 2) * dA[p]
                p_e += ip_
                m2 += ip_ * ((p + 0.5) * dr) ** 2
            beam_out[2 + e] = field[0].real ** 2 + field[0].imag ** 2
            beam_out[2 + n_edges + e] = math.sqrt(2.0 * m2 / p_e) if p_e > 0.0 else 0.0
        for z in range(nr):
            area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            s = 0.0
            if prop_on:
                for p in range(npg):
                    if ov[p, z] > 0.0:
                        s += (field[p].real ** 2 + field[p].imag ** 2) * ov[p, z]
            ring_i[z] = s / area
        for z in prange(nr):
            st = np.empty(6)
            acc = np.zeros(9)
            area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            npart = heavy[e, z] - cx2p[e, z] - cs6[4, e, z]
            tg = eg[e, z] / (0.5 * dof * KB * npart) if npart > 0 else tg0
            n_gas = max((heavy[e, z] - cxp[e, z] - 2.0 * cx2p[e, z] - cs6[5, e, z] - cs6[3, e, z]
                         - 2.0 * cs6[4, e, z]) / area, 1e10)
            inten = ring_i[z] if ring_i[z] >= 1.0 else 0.0
            for q in range(6):
                st[q] = cs6[q, e, z] / area
            p_ext = ring_work[e, z] / (QE * dt)
            ring_work[e, z] = 0.0
            kinetic_cell_cr(st, dt, n_gas, tg, p_ext, inten, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab, acc,
                            path_len)
            for q in range(6):
                cs6[q, e, z] = st[q] * area
            eg[e, z] += acc[6] * area
            for qq in range(9):
                acc_r[z, qq] = acc[qq] * area
            frac[z] = 0.0
            if inten > 0.0:
                frac[z] = min(acc[1] * path_len / (inten * dt), 1.0)
        lsum = 0.0
        for z in range(nr):
            lsum += acc_r[z, 1]
            for qq in range(9):
                ledger[e, qq] += acc_r[z, qq]
        laser_in[e] += lsum
        ledger[e, 11] += lsum
        if not prop_on:
            continue
        # absorption: remove exactly the energy deposited in each ring
        for p in range(npg):
            loss = 0.0
            for z in range(nr):
                if ov[p, z] > 0.0:
                    loss += frac[z] * ov[p, z]
            g = 1.0 - loss / dA[p]
            field[p] *= math.sqrt(max(g, 0.0))
        # refraction profile after the kinetics: plasma, neutral depletion, Kerr
        for p in range(npg):
            ne_p = 0.0
            ng_p = 0.0
            cov = 0.0
            for z in range(nr):
                if ov[p, z] > 0.0:
                    area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
                    ne_p += (cxp[e, z] + cx2p[e, z] + 2.0 * cs6[5, e, z]) / area * ov[p, z]
                    ng_p += (heavy[e, z] - cxp[e, z] - 2 * cx2p[e, z]) / area * ov[p, z]
                    cov += ov[p, z]
            ng = (ng_p + (dA[p] - cov) * n0) / dA[p]
            dn[p] = -0.5 * (ne_p / dA[p]) / n_crit + refr * (ng / n0 - 1.0)
        dmax = 0.0
        for p in range(npg):
            ip_ = field[p].real ** 2 + field[p].imag ** 2
            dmax = max(dmax, abs(dn[p]) + n2 * ip_)
        nsub = max(8, int(math.ceil(k_l * dmax * path_len / 0.2)))
        dz = path_len / nsub
        for _ in range(nsub):
            for p in range(npg):
                ip_ = field[p].real ** 2 + field[p].imag ** 2
                ph = 0.5 * k_l * (dn[p] + n2 * ip_) * dz
                field[p] *= complex(math.cos(ph), math.sin(ph))
            _cn_radial(field, dr, dz, k_l, wa, wb, wc, wd, wx)
            for p in range(npg):
                ip_ = field[p].real ** 2 + field[p].imag ** 2
                ph = 0.5 * k_l * (dn[p] + n2 * ip_) * dz
                field[p] *= complex(math.cos(ph), math.sin(ph))
            lost = 0.0
            for p in range(npg):
                if lat_mask[p] < 1.0:
                    i0 = field[p].real ** 2 + field[p].imag ** 2
                    field[p] *= lat_mask[p]
                    lost += i0 * (1.0 - lat_mask[p] ** 2) * dA[p]
            beam_out[0] += lost * dt
    if prop_on:
        pw = 0.0
        for p in range(npg):
            pw += (field[p].real ** 2 + field[p].imag ** 2) * dA[p]
        beam_out[1] = pw
    for e in prange(n_edges):
        _fil_diffuse_cr(e, rb, heavy, cs6, dt, prm, tabs, ln_e0, h)


@njit(parallel=True, cache=True)
def filament_drude_cr(rb, cs6, heavy, dt_fdtd, dx, prm, tabs, ln_e0, h, w_out, a_out, b_out):
    n_edges, nb = rb.shape
    nr = nb - 1
    for e in prange(n_edges):
        for z in range(nr):
            area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            ne = (cs6[0, e, z] + cs6[1, e, z] + 2.0 * cs6[5, e, z]) / area
            n_gas = max((heavy[e, z] - cs6[0, e, z] - 2 * cs6[1, e, z] - cs6[5, e, z] - cs6[3, e, z]
                         - 2 * cs6[4, e, z]) / area, 1e10)
            wp2, nu = drude_params(ne, cs6[2, e, z] / area, n_gas, tabs, ln_e0, h, prm)
            d = 1.0 + 0.5 * nu * dt_fdtd
            a_out[e, z] = (1.0 - 0.5 * nu * dt_fdtd) / d
            b_out[e, z] = EPS0 * wp2 * dt_fdtd / d
            w_out[e, z] = area / (dx * dx)
