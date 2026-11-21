"""Numba kernels for local electron/ion kinetics, sub-cell laser filament, and radial hydro.

Species per cell: X+ (atomic or parent molecular ion), X2+ (dimer ion), electrons by
quasi-neutrality (ne = nX+ + nX2+), electron energy density W (eV m^-3).

Rates come from log-uniform tables in mean electron energy eps (eV). tabs has 4 rows, or 8 rows
(rows 0-3 two-term Boltzmann, rows 4-7 Maxwellian); with 8 rows and prm[P_EE_SCALE] > 0 each
rate is blended linearly, k = (1-w) kB + w kM, with the local, parameter-free weight
  w = nu_ee / (nu_ee + nu_eps)
  nu_ee  = NRL Plasma Formulary electron collision rate 2.91e-12 ne lnL Te^-1.5 (SI)
  nu_eps = N [kexc/eps + kel + kiz Ip/eps]_Boltzmann   (electron energy-relaxation frequency)
i.e. the EEDF is Maxwellized where e-e energy exchange outpaces energy loss to neutrals
(w = 1/2 at nu_ee = nu_eps). prm[P_EE_SCALE] (default 1) multiplies nu_ee and is only used to
bound the O(1) uncertainty of this criterion in sensitivity runs. Tables:
  km   momentum-transfer rate coefficient (m^3/s)
  kiz  ground-state ionization rate coefficient (m^3/s)
  kexc excitation/vibration energy-loss coefficient, ionization excluded (eV m^3/s)
  kel  elastic energy-loss coefficient: loss = N ne kel (eps - 1.5 k Tg)  (m^3/s)

Energy bookkeeping per cell (exact for each substep, by construction):
  electron input = dW + E_exc + E_el + Ip*R_iz*dt + eps*R_rec*dt
  laser input    = IB heating (in electron input) + R_ppt*K*h nu*dt
Gas receives E_el + eta_fast*(E_exc + released recombination energy); the rest is
recorded as a radiative/metastable sink.
"""

import math

import numpy as np
from numba import njit, prange

QE = 1.602176634e-19
ME = 9.1093837015e-31
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
C0 = 299792458.0
EV_K = QE / KB
N_LOSCHMIDT = 2.6867811e25

# parameter vector layout
P_MASS, P_IP, P_DOF = 0, 1, 2
P_KCONV0, P_KCONV_EXP = 3, 4
P_B1_0, P_B1_EXP, P_B2_0, P_B2_EXP, P_B2_THETA, P_B2_C = 5, 6, 7, 8, 9, 10
P_K3B0, P_K3B_EXP = 11, 12
P_MU_XP, P_MU_X2P = 13, 14
P_ETA = 15
P_KDISS0, P_EDISS = 16, 17
P_KDEX0, P_EDEX = 18, 19
P_PHOTON_EV, P_OMEGA_L, P_NPHOT = 20, 21, 22
P_B1N_0, P_B1N_EXP, P_B2N_0, P_B2N_EXP = 23, 24, 25, 26   # neutral-assisted: b += b0 (Te_K/300)^x N
P_DISS_A, P_DISS_T1 = 27, 28                              # kdiss = k0 Tg_eV^a exp(-E/Tg_eV) exp(Tg/T1)
P_EE_SCALE = 29                                           # EEDF closure: nu_ee multiplier (1 = criterion)
N_PARAMS = 30


@njit(cache=True, inline="always")
def _lookup(table, ln_e0, h, eps):
    x = (math.log(eps) - ln_e0) / h
    n = table.shape[0]
    if x <= 0.0:
        return table[0]
    if x >= n - 1:
        return table[n - 1]
    i = int(x)
    f = x - i
    return table[i] * (1.0 - f) + table[i + 1] * f


EPS_LO = 1e-3  # eV, lower bound of the electron mean energy in the implicit solve
# m^-3; numerical floor for empty cells. Electrons created by the floor avalanche like real ones, so it is kept
# ~1e-39 electrons per mm^3 cell: no spurious volume breakdown within ~70 ns even at nu_iz ~ 1.3e9 s^-1
NE_FLOOR = 1e-30


@njit(cache=True, inline="always")
def _rate(tabs, row, wmx, ln_e0, h, eps):
    v = _lookup(tabs[row], ln_e0, h, eps)
    if wmx > 0.0 and tabs.shape[0] >= 8:
        v = (1.0 - wmx) * v + wmx * _lookup(tabs[row + 4], ln_e0, h, eps)
    return v


@njit(cache=True, inline="always")
def nu_ei(ne, te_ev):
    """NRL formulary electron collision rate (SI): 2.91e-12 ne lnL Te^-1.5."""
    if ne <= 0.0:
        return 0.0
    te = max(te_ev, 0.02)
    ln_l = 23.0 - math.log(math.sqrt(ne * 1e-6) * te**-1.5)
    if ln_l < 2.0:
        ln_l = 2.0
    return 2.91e-12 * ne * ln_l * te**-1.5


@njit(cache=True, inline="always")
def _wmax(ne, n_neu, eps, tabs, ln_e0, h, ip, scale):
    """Maxwellian weight nu_ee / (nu_ee + nu_eps); nu_eps from the Boltzmann rows (no recursion)."""
    if scale <= 0.0 or tabs.shape[0] < 8 or ne <= 0.0:
        return 0.0
    e = max(eps, EPS_LO)
    nee = scale * nu_ei(ne, e / 1.5)
    nue = 0.0
    if n_neu > 0.0:
        nue = n_neu * (_lookup(tabs[2], ln_e0, h, e) / e + _lookup(tabs[3], ln_e0, h, e)
                       + _lookup(tabs[1], ln_e0, h, e) * ip / e)
    if nee + nue <= 0.0:
        return 0.0
    return nee / (nee + nue)


@njit(cache=True, inline="always")
def _eloss(e, ne, n_gas, tabs, ln_e0, h, eps_g, two_m_over_m, wmx):
    """Electron energy loss per electron (eV/s): inelastic + elastic (neutrals and ions)."""
    exc = n_gas * _rate(tabs, 2, wmx, ln_e0, h, e)
    el = (n_gas * _rate(tabs, 3, wmx, ln_e0, h, e) + two_m_over_m * nu_ei(ne, e / 1.5)) * (e - eps_g)
    return exc + el


@njit(cache=True)
def kinetic_cell(state, dt, n_gas, tg, p_ext, intensity, prm, tabs, ln_e0, h,
                 ln_i0, hi, ppt_tab, out, path_len=0.0):
    """Advance one cell over dt. state = [nXp, nX2p, W]; out accumulates
    [E_in_rf, E_laser, E_exc, E_el, E_iz, E_rec_carry, E_gas, E_sink, chem] in J/m^3.

    intensity is the value entering the cell; if path_len > 0 the laser terms are reduced by the
    Beer-Lambert cell average (1 - exp(-alpha L)) / (alpha L), so absorption never exceeds the
    incident flux (intensity / path_len per unit volume).

    n_gas is the neutral density at entry; neutrals are depleted inside the step (nuclei
    conserved: n_neutral + nXp + 2 nX2p = const). Electrons may cool below the gas temperature
    (inelastic losses, no superelastic return); the elastic term then reheats them from the gas,
    implicitly and inside the ledger (negative elastic / gas_heat)."""
    km_t = tabs[0]; kiz_t = tabs[1]; kexc_t = tabs[2]
    eps_max = math.exp(ln_e0 + h * (km_t.shape[0] - 1))
    mass = prm[P_MASS]
    ip = prm[P_IP]
    eta = prm[P_ETA]
    eps_g = 1.5 * tg / EV_K
    two_m_over_m = 2.0 * ME / mass
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
        ex = math.exp(-theta / tg)
        b2_tg = (1.0 - ex) / (1.0 - prm[P_B2_C] * ex)
    t = 0.0
    nxp = state[0]; nx2p = state[1]; w = state[2]
    n_tot = n_gas + nxp + 2.0 * nx2p
    iters = 0
    while t < dt and iters < 200000:
        iters += 1
        ne = nxp + nx2p
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
        n_neu = max(n_tot - nxp - 2.0 * nx2p, 0.0)
        kconv = kconv_c * n_neu * n_neu
        kdiss = kdiss_c * n_neu
        wmx = _wmax(ne, n_neu, eps, tabs, ln_e0, h, ip, prm[P_EE_SCALE])
        km = _rate(tabs, 0, wmx, ln_e0, h, eps)
        kiz = _rate(tabs, 1, wmx, ln_e0, h, eps)
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
        p_ib = 0.0
        if intensity > 0.0:
            p_ib = QE * nu * intensity / (ME * EPS0 * C0 * omega_l * omega_l)  # eV/s per electron
            if path_len > 0.0:
                a_l = (ne * p_ib + r_ppt * nphot * prm[P_PHOTON_EV]) * QE / intensity * path_len
                if a_l > 1e-6:
                    fa = (1.0 - math.exp(-a_l)) / a_l
                    p_ib *= fa
                    r_ppt *= fa
        heat_per_e = p_ext / ne
        # X+ <-> X2+ exchange is solved fully implicitly below (charge-conserving for any dts),
        # so its (often ps-scale) rates do not limit the substep
        # kiz*ne, w_ppt: neutral depletion; ip*r_iz/w: ionization cost may not overdraw W
        rate = max(kiz * n_neu, kiz * ne, w_ppt, l1, l2, 1e3)
        relax = _eloss(eps, ne, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx) / eps
        dts = min(dt - t, 0.05 / rate)
        if r_iz > 0.0:
            dts = min(dts, 0.3 * max(w, ne * EPS_LO) / (ip * r_iz))
        if heat_per_e + p_ib > 0.0:
            dts = min(dts, max(0.3 * eps / (heat_per_e + p_ib), 0.05 / max(relax, 1.0)))
        if dts < 1e-17:
            dts = 1e-17
        a11 = 1.0 + dts * (kconv + l1)
        a12 = -dts * d1
        a21 = -dts * kconv
        a22 = 1.0 + dts * (d1 + l2)
        r1 = nxp + dts * (r_iz + r_ppt)
        det = a11 * a22 - a12 * a21
        nxp_n = (r1 * a22 - a12 * nx2p) / det
        nx2p_n = (a11 * nx2p - a21 * r1) / det
        ne_n = nxp_n + nx2p_n
        r_rec = (ne + dts * (r_iz + r_ppt) - ne_n) / dts
        if r_rec < 0.0:
            r_rec = 0.0
        src = p_ext + ne * p_ib + excess * r_ppt - ip * r_iz  # explicit sources, eV m^-3 s^-1
        rhs = w + dts * src
        # implicit loss: g(e) = ne_n e + dts[ne_n (exc(e) + el(e)) + e r_rec] - rhs = 0
        lo = EPS_LO
        hi_e = eps_max
        e1 = min(max(eps, lo), hi_e)
        ok = False
        for _ in range(8):
            g = ne_n * e1 + dts * (ne_n * _eloss(e1, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx) + e1 * r_rec) - rhs
            de = e1 * 1e-4 + 1e-12
            e2 = e1 + de
            g2 = ne_n * e2 + dts * (ne_n * _eloss(e2, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx) + e2 * r_rec) - rhs
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
                g = ne_n * mid + dts * (ne_n * _eloss(mid, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx) + mid * r_rec) - rhs
                if g > 0:
                    b = mid
                else:
                    a = mid
            e1 = math.sqrt(a * b)
        e_exc = dts * ne_n * n_neu * _rate(tabs, 2, wmx, ln_e0, h, e1)
        e_el = dts * ne_n * (_eloss(e1, ne_n, n_neu, tabs, ln_e0, h, eps_g, two_m_over_m, wmx)
                             - n_neu * _rate(tabs, 2, wmx, ln_e0, h, e1))
        e_rc = dts * e1 * r_rec
        w_new = rhs - e_exc - e_el - e_rc
        e_iz = dts * ip * r_iz
        rel = dts * ip * r_rec  # chemical energy released by recombination
        out[0] += dts * p_ext * QE
        out[1] += dts * (ne * p_ib + r_ppt * nphot * prm[P_PHOTON_EV]) * QE
        out[2] += e_exc * QE
        out[3] += e_el * QE
        out[4] += e_iz * QE
        out[5] += e_rc * QE
        out[6] += (e_el + eta * (e_exc + rel + e_rc)) * QE
        out[7] += ((1.0 - eta) * (e_exc + rel + e_rc) + dts * r_ppt * (nphot * prm[P_PHOTON_EV] - ip - excess)) * QE
        out[8] += dts * ip * (r_iz + r_ppt - r_rec) * QE
        nxp = nxp_n; nx2p = nx2p_n; w = w_new
        t += dts
    state[0] = nxp; state[1] = nx2p; state[2] = w
    return iters


@njit(cache=True, inline="always")
def drude_params(ne, w, n_gas, tabs, ln_e0, h, prm):
    eps = w / ne if ne > 0 else 0.04
    eps = max(eps, 1e-3)
    wmx = _wmax(ne, n_gas, eps, tabs, ln_e0, h, prm[P_IP], prm[P_EE_SCALE])
    nu = n_gas * _rate(tabs, 0, wmx, ln_e0, h, eps) + nu_ei(ne, eps / 1.5)
    wp2 = ne * QE * QE / (EPS0 * ME)
    return wp2, nu


@njit(parallel=True, cache=True)
def bulk_update(nxp, nx2p, w, active, node_work, dt, n_gas, tg, prm, tabs, ln_e0, h,
                ln_i0, hi, ppt_tab, wp2_out, nu_out, ledger):
    """Advance every active node; node_work holds the RF work (J/m^3) of the last block.

    ledger (nthreads-independent): per-x-slab accumulators summed by the caller, shape (nx+1, 9)."""
    nxp1, nyp1, nzp1 = nxp.shape
    for i in prange(nxp1):
        st = np.empty(3)
        acc = np.zeros(9)
        for j in range(nyp1):
            for k in range(nzp1):
                if not active[i, j, k]:
                    wp2_out[i, j, k] = 0.0
                    nu_out[i, j, k] = 0.0
                    continue
                st[0] = nxp[i, j, k]; st[1] = nx2p[i, j, k]; st[2] = w[i, j, k]
                p_ext = node_work[i, j, k] / (QE * dt)
                n_neu = max(n_gas - st[0] - 2.0 * st[1], 1e10)
                kinetic_cell(st, dt, n_neu, tg, p_ext, 0.0, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab, acc)
                nxp[i, j, k] = st[0]; nx2p[i, j, k] = st[1]; w[i, j, k] = st[2]
                ne = st[0] + st[1]
                wp2, nu = drude_params(ne, st[2], max(n_gas - st[0] - 2.0 * st[1], 1e10), tabs, ln_e0, h, prm)
                wp2_out[i, j, k] = wp2
                nu_out[i, j, k] = nu
        for q in range(9):
            ledger[i, q] += acc[q]


@njit(cache=True)
def _thomas(a, b, c, d, x):
    n = b.shape[0]
    cp = np.empty(n); dp = np.empty(n)
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]


@njit(cache=True)
def _diffuse_counts(cnt, rb, area, coef_face, dt):
    """Implicit radial diffusion of a per-length count with face coefficient D (m^2/s).

    Flux through face f (between zones f-1, f) = -D 2 pi r_f (n_f - n_{f-1}) / dr_c. Conservative."""
    nz = cnt.shape[0]
    a = np.zeros(nz); b = np.zeros(nz); c = np.zeros(nz); d = np.zeros(nz)
    rc = np.empty(nz)
    for z in range(nz):
        rc[z] = 0.5 * (rb[z] + rb[z + 1])
    for z in range(nz):
        b[z] = 1.0
        d[z] = cnt[z]
    for f in range(1, nz):
        g = coef_face[f] * 2.0 * math.pi * rb[f] / (rc[f] - rc[f - 1]) * dt
        # count change in zone f-1: + g (n_f - n_{f-1}) with n = cnt/area
        b[f - 1] += g / area[f - 1]
        c[f - 1] -= g / area[f]
        b[f] += g / area[f]
        a[f] -= g / area[f - 1]
    out = np.empty(nz)
    _thomas(a, b, c, d, out)
    for z in range(nz):
        cnt[z] = max(out[z], 0.0)


@njit(cache=True)
def _conduct_te(te, cnt, rb, kap_face, dt):
    """Implicit conduction 1.5 N_z dTe_z/dt = sum_f G_f (Te_nb - Te_z); conserves sum 1.5 N Te."""
    nz = te.shape[0]
    a = np.zeros(nz); b = np.zeros(nz); c = np.zeros(nz); d = np.zeros(nz)
    for z in range(nz):
        cap = 1.5 * max(cnt[z], 1e-30)
        b[z] = cap
        d[z] = cap * te[z]
    for f in range(1, nz):
        rc_l = 0.5 * (rb[f - 1] + rb[f]); rc_r = 0.5 * (rb[f] + rb[f + 1])
        g = kap_face[f] * 2.0 * math.pi * rb[f] / (rc_r - rc_l) * dt
        b[f - 1] += g; c[f - 1] -= g
        b[f] += g; a[f] -= g
    _thomas(a, b, c, d, te)


@njit(cache=True)
def _fil_hydro(e, rb, ub, heavy, eg, cx2p, ue, dt, p0, mass, dof, gam, ledger, hydro_on):
    n_edges, nb = rb.shape
    nr = nb - 1
    area = np.empty(nr)
    for z in range(nr):
        area[z] = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
    if hydro_on:
        t = 0.0
        while t < dt:
            p_tot = np.empty(nr); q = np.empty(nr); ph = np.empty(nr); pe = np.empty(nr)
            dtc = dt - t
            for z in range(nr):
                npart = heavy[e, z] - cx2p[e, z]
                ph[z] = eg[e, z] / (0.5 * dof * area[z])
                pe[z] = (2.0 / 3.0) * ue[e, z] * QE / area[z]
                rho = heavy[e, z] * mass / area[z]
                du = ub[e, z + 1] - ub[e, z]
                cs = math.sqrt(gam * (ph[z] + pe[z]) / rho)
                q[z] = rho * (2.0 * du * du + 0.5 * cs * abs(du)) if du < 0.0 else 0.0
                p_tot[z] = ph[z] + pe[z] + q[z]
                dr = rb[e, z + 1] - rb[e, z]
                dtc = min(dtc, 0.25 * dr / (cs + abs(du) + 1e-9))
            if dtc < 1e-16:
                dtc = 1e-16
            # Energy-compatible update: boundaries move with the time-centred velocity ubar and
            # the zone work is p * 2 pi (r_old ubar)|_faces dt, which equals the kinetic-energy
            # change of the nodes (sum of mb ubar du) exactly; the outer face works against p0.
            ke_old = 0.0
            ke_new = 0.0
            ubar = np.zeros(nb)
            for b in range(1, nb):
                pr = p_tot[b] if b < nr else p0
                mb = 0.5 * heavy[e, b - 1] * mass + (0.5 * heavy[e, b] * mass if b < nr else 0.0)
                u_old = ub[e, b]
                ub[e, b] = u_old - dtc * 2.0 * math.pi * rb[e, b] * (pr - p_tot[b - 1]) / mb
                ubar[b] = 0.5 * (u_old + ub[e, b])
                ke_old += 0.5 * mb * u_old ** 2
                ke_new += 0.5 * mb * ub[e, b] ** 2
            for z in range(nr):
                d_w = 2.0 * math.pi * (rb[e, z + 1] * ubar[z + 1] - rb[e, z] * ubar[z]) * dtc
                eg[e, z] -= (ph[z] + q[z]) * d_w
                ue[e, z] -= pe[z] * d_w / QE
                if eg[e, z] < 1e-30:
                    eg[e, z] = 1e-30
                if ue[e, z] < 0.0:
                    ue[e, z] = 0.0
            ledger[e, 10] += p0 * 2.0 * math.pi * rb[e, nr] * ubar[nr] * dtc
            for b in range(1, nb):
                rb[e, b] += ubar[b] * dtc
            for z in range(nr):
                area[z] = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            ledger[e, 9] += ke_new - ke_old
            t += dtc


@njit(cache=True)
def _fil_diffuse(e, rb, heavy, cxp, cx2p, ue, dt, tg0, prm, tabs, ln_e0, h):
    n_edges, nb = rb.shape
    nr = nb - 1
    area = np.empty(nr)
    for z in range(nr):
        area[z] = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
    # ---- 3. ambipolar diffusion + electron heat conduction (implicit, conservative) ----
    te_floor = EPS_LO / 1.5  # same bound as kinetic_cell: no energy created here
    te_r = np.empty(nr)
    for z in range(nr):
        cnt = cxp[e, z] + cx2p[e, z]
        te_r[z] = max(ue[e, z] / (1.5 * cnt), te_floor) if cnt > 0.0 else te_floor
    dxp = np.zeros(nr + 1); dx2p = np.zeros(nr + 1); kap = np.zeros(nr + 1)
    for f in range(1, nr):
        nel = (cxp[e, f - 1] + cx2p[e, f - 1]) / area[f - 1]
        ner = (cxp[e, f] + cx2p[e, f]) / area[f]
        te_f = 0.5 * (te_r[f - 1] + te_r[f])
        ng_f = 0.5 * (heavy[e, f - 1] / area[f - 1] + heavy[e, f] / area[f])
        dxp[f] = prm[P_MU_XP] * N_LOSCHMIDT / ng_f * te_f
        dx2p[f] = prm[P_MU_X2P] * N_LOSCHMIDT / ng_f * te_f
        nmin = min(nel, ner)
        wmx_f = _wmax(nmin, ng_f, 1.5 * te_f, tabs, ln_e0, h, prm[P_IP], prm[P_EE_SCALE])
        nu_f = ng_f * _rate(tabs, 0, wmx_f, ln_e0, h, 1.5 * te_f) + nu_ei(nmin, te_f)
        if nu_f > 0.0 and nmin > 0.0:
            # electron heat flux q = -(5/2) n D_e grad(Te), D_e = e Te / (m nu)  [eV m^-2 s^-1]
            kap[f] = 2.5 * nmin * te_f / (ME * nu_f) * QE
    _diffuse_counts(cxp[e], rb[e], area, dxp, dt)
    _diffuse_counts(cx2p[e], rb[e], area, dx2p, dt)
    # electrons carry their mean energy with them, then conduct heat on Te
    cnt_new = np.empty(nr)
    for z in range(nr):
        cnt_new[z] = cxp[e, z] + cx2p[e, z]
    _conduct_te(te_r, cnt_new, rb[e], kap, dt)
    for z in range(nr):
        ue[e, z] = 1.5 * te_r[z] * cnt_new[z]


@njit(parallel=True, cache=True)
def filament_update(rb, ub, heavy, eg, cxp, cx2p, ue, ring_work, trans, order, laser_in,
                    dt, n0, tg0, p0, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab,
                    intensity_axis, beam_w, ledger, hydro_on, path_len):
    """Advance all filament edges (laser on axis, sub-cell Lagrangian rings).

    rb, ub: (nE, R+1) ring boundaries (m) and velocities; heavy: heavy-particle count per length
    per ring (constant); eg: gas internal energy (J/m); cxp, cx2p: ion counts per length;
    ue: electron energy (eV/m). ring_work: RF work density (J/m^3) of the last block per ring.
    intensity_axis[e]: unattenuated on-axis laser intensity (W/m^2) at the block midpoint;
    beam_w[e]: 1/e^2 beam radius at the edge.
    The laser is marched through the edges in propagation order (order[0] first) within the
    block: trans[0] (in/out, set to 1 by the caller) is the beam power transmission, reduced
    after each edge by its absorbed energy / (unattenuated beam power x dt). The beam keeps its
    Gaussian radial shape; inside each ring the Beer-Lambert cell average (path_len = edge length)
    keeps absorption below the incident flux, so depletion is same-block and energy-conserving.
    laser_in[e] accumulates the laser energy per length absorbed during the block.
    ledger: (nE, 12) accumulators [9 kinetic channels, hydro_ke, p0_work, laser_absorbed].
    """
    n_edges, nb = rb.shape
    nr = nb - 1
    mass = prm[P_MASS]
    dof = prm[P_DOF]
    gam = 1.0 + 2.0 / dof
    # ---- 1. Lagrangian hydro (VNR artificial viscosity), electron + heavy pressure ----
    for e in prange(n_edges):
        _fil_hydro(e, rb, ub, heavy, eg, cx2p, ue, dt, p0, mass, dof, gam, ledger, hydro_on)
    # ---- 2. local kinetics with laser and RF heating; laser marched edge by edge ----
    acc_r = np.zeros((nr, 9))
    for ii in range(n_edges):
        e = order[ii]
        w2 = beam_w[e] ** 2
        for z in prange(nr):
            st = np.empty(3)
            acc = np.zeros(9)
            area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            npart = heavy[e, z] - cx2p[e, z]
            tg = eg[e, z] / (0.5 * dof * KB * npart) if npart > 0 else tg0
            n_gas = max((heavy[e, z] - cxp[e, z] - 2.0 * cx2p[e, z]) / area, 1e10)
            inten = 0.0
            if intensity_axis[e] > 0.0 and trans[0] > 0.0:
                r1 = rb[e, z]; r2 = rb[e, z + 1]
                frac = 0.5 * math.pi * w2 * (math.exp(-2 * r1 * r1 / w2) - math.exp(-2 * r2 * r2 / w2)) / area
                inten = intensity_axis[e] * trans[0] * frac
                if inten < 1.0:  # W/m^2; avoids denormal 0/0 in the far Gaussian tail
                    inten = 0.0
            st[0] = cxp[e, z] / area; st[1] = cx2p[e, z] / area; st[2] = ue[e, z] / area
            p_ext = ring_work[e, z] / (QE * dt)
            ring_work[e, z] = 0.0
            kinetic_cell(st, dt, n_gas, tg, p_ext, inten, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab, acc,
                         path_len)
            cxp[e, z] = st[0] * area; cx2p[e, z] = st[1] * area; ue[e, z] = st[2] * area
            eg[e, z] += acc[6] * area
            for qq in range(9):
                acc_r[z, qq] = acc[qq] * area   # J per metre of filament
        lsum = 0.0
        for z in range(nr):
            lsum += acc_r[z, 1]
            for qq in range(9):
                ledger[e, qq] += acc_r[z, qq]
        laser_in[e] += lsum
        ledger[e, 11] += lsum
        p_beam = intensity_axis[e] * 0.5 * math.pi * w2 * dt
        if p_beam > 1e-200 and path_len > 0.0:
            trans[0] = max(trans[0] - lsum * path_len / p_beam, 0.0)
    # ---- 3. ambipolar diffusion + electron heat conduction (implicit, conservative) ----
    for e in prange(n_edges):
        _fil_diffuse(e, rb, heavy, cxp, cx2p, ue, dt, tg0, prm, tabs, ln_e0, h)


@njit(parallel=True, cache=True)
def filament_drude(rb, cxp, cx2p, ue, heavy, dt_fdtd, dx, prm, tabs, ln_e0, h, w_out, a_out, b_out):
    n_edges, nb = rb.shape
    nr = nb - 1
    for e in prange(n_edges):
        for z in range(nr):
            area = math.pi * (rb[e, z + 1] ** 2 - rb[e, z] ** 2)
            ne = (cxp[e, z] + cx2p[e, z]) / area
            n_gas = max((heavy[e, z] - cxp[e, z] - 2 * cx2p[e, z]) / area, 1e10)
            wp2, nu = drude_params(ne, ue[e, z] / area, n_gas, tabs, ln_e0, h, prm)
            d = 1.0 + 0.5 * nu * dt_fdtd
            a_out[e, z] = (1.0 - 0.5 * nu * dt_fdtd) / d
            b_out[e, z] = EPS0 * wp2 * dt_fdtd / d
            w_out[e, z] = area / (dx * dx)
