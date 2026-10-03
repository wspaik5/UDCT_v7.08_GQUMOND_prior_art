#!/usr/bin/env python3
"""
UDCT v7.08 -- reproducible checks of the analytic statements in the
"prior-art correction" note (GQUMOND Hessian term / compactness-type screening).

IMPORTANT SCOPE NOTE
--------------------
UDCT v7.08 is a correction note with NO new calculation. This script therefore
does NOT reproduce any new result. It independently re-derives / numerically
checks the *elementary statements* the note relies on:

  C1  Point-mass Hessian identity: u = |grad psi|^2 / (l0^2 psi_ij psi_ij) = r^2 / (6 l0^2)
      and the screening function f(u) = u/(1+u) limits.
  C2  Chain rule for Milgrom's example P = f(u) Q[(grad psi)^2/(f a0^2)]:
      dP/dpsi_ij = f'(u) [Q - s Q'] du/dpsi_ij, and the Hessian-term factor
      u f'(u)/f(u) = 1/(1+u)  -> small where f ~ 1 ("scales with the local
      logarithmic derivative of f").
  C3  QUMOND phantom density for two equal point masses: is it negative somewhere
      near the midpoint (null point)? What is the sign of the phantom mass enclosed
      in a small sphere around the null point?
  C4  Gauss's law: an outward angle-averaged radial force <=> negative enclosed
      effective mass; arithmetic for the "10-483 x" range quoted in the note.
  C5  Dimensionless-constant remark: l0 = X * c^2/a0 requires X << 1.

NOT checked here (needs the original sources / the v7.06 released script):
  * Equation numbers / exact form of Milgrom (2023) Eqs. 7 and 14 (read arXiv:2305.01589).
  * v7.06 Table 1 values (10-483 x), which come from the v7.06 script and data.
  * The "~6000 AU" solar balance scale quoted in the note.
  * Whether QUMOND/GQUMOND can make the enclosed effective mass negative in a strong
    external field (this is the v7.09 diagnostic, not part of v7.08).

Usage:
    python udct_v7_08_checks.py                 # run all checks
    python udct_v7_08_checks.py --nu standard   # use the "standard" nu instead of "simple"
    python udct_v7_08_checks.py --json out.json # also write machine-readable results

Only numpy is required.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------
# Physical constants (SI) -- only used in C5
# --------------------------------------------------------------------------
C_LIGHT = 299_792_458.0  # m/s
A0_MOND = 1.2e-10  # m/s^2, conventional MOND acceleration scale
AU = 1.495978707e11  # m

RESULTS: dict = {}
FAILED: list[str] = []


def report(name: str, ok: bool, detail: str = "") -> None:
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


# --------------------------------------------------------------------------
# C1: point-mass Hessian identity and screening function
# --------------------------------------------------------------------------
def screening_f(u):
    """Illustrative screening function used in the note: f(u) = u / (1 + u)."""
    u = np.asarray(u, dtype=float)
    return u / (1.0 + u)


def point_mass_potential(x, GM=1.0):
    r = np.linalg.norm(x)
    return -GM / r


def numeric_gradient_hessian(phi, x, h=1e-4):
    """Central-difference gradient and Hessian of a scalar function phi at x."""
    x = np.asarray(x, dtype=float)
    n = x.size
    grad = np.zeros(n)
    hess = np.zeros((n, n))
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h
        grad[i] = (phi(x + ei) - phi(x - ei)) / (2 * h)
        for j in range(n):
            ej = np.zeros(n)
            ej[j] = h
            hess[i, j] = (
                phi(x + ei + ej) - phi(x + ei - ej) - phi(x - ei + ej) + phi(x - ei - ej)
            ) / (4 * h * h)
    return grad, hess


def check_c1():
    print("\nC1  Point-mass Hessian identity and screening-function limits")
    GM, l0 = 1.0, 1.0
    rows = []
    max_rel_err = 0.0
    for r in [0.3, 1.0, 3.0]:
        x = np.array([r * 0.6, r * 0.0, r * 0.8])  # arbitrary direction, |x| = r
        grad, hess = numeric_gradient_hessian(lambda y: point_mass_potential(y, GM), x, h=1e-3 * r)
        g2 = float(grad @ grad)
        hh = float(np.sum(hess * hess))
        u_num = g2 / (l0**2 * hh)
        u_ana = r**2 / (6.0 * l0**2)
        rel = abs(u_num - u_ana) / u_ana
        max_rel_err = max(max_rel_err, rel)
        rows.append(dict(r=r, u_numeric=u_num, u_analytic=u_ana, rel_err=rel, f=float(screening_f(u_ana))))
    for row in rows:
        print(
            f"      r/l0={row['r']:<4}  u_num={row['u_numeric']:.6f}  u_ana={row['u_analytic']:.6f}"
            f"  f(u)={row['f']:.4f}"
        )
    report("u = r^2/(6 l0^2) for a point mass (|psi_ij psi_ij| = 6 G^2M^2/r^6)", max_rel_err < 1e-3,
           f"max rel. err = {max_rel_err:.2e}")
    report("f(u) -> 0 for u << 1 (screening ON in small systems)", float(screening_f(1e-6)) < 1e-5)
    report("f(u) -> 1 for u >> 1 (standard QUMOND recovered)", abs(float(screening_f(1e6)) - 1) < 1e-5)
    # point-mass screening radius: u = 1  <=>  r = sqrt(6) l0
    report("u = 1 at r = sqrt(6) l0 for a point mass", abs(screening_f(1.0) - 0.5) < 1e-12)
    RESULTS["C1"] = rows


# --------------------------------------------------------------------------
# C2: chain rule and logarithmic-derivative factor
# --------------------------------------------------------------------------
def P_function(grad, hess, a0, l0, Q):
    """Milgrom-type example: P = f(u) Q[(grad psi)^2 / (f a0^2)], u = |grad|^2/(l0^2 H), H = psi_ij psi_ij."""
    g2 = float(grad @ grad)
    H = float(np.sum(hess * hess))
    u = g2 / (l0**2 * H)
    f = u / (1.0 + u)
    s = g2 / (f * a0**2)
    return f * Q(s)


def check_c2():
    print("\nC2  Chain rule for P = f(u) Q[(grad psi)^2/(f a0^2)] and the factor u f'/f")
    rng = np.random.default_rng(12345)

    def Q(s):  # arbitrary smooth test function (the identity does not depend on its form)
        return np.sqrt(1.0 + s) - 1.0 + 0.3 * s

    def dQ(s):
        return 0.5 / np.sqrt(1.0 + s) + 0.3

    a0, l0 = 1.0, 1.0
    max_err = 0.0
    for _ in range(5):
        grad = rng.normal(size=3)
        A = rng.normal(size=(3, 3))
        hess = 0.5 * (A + A.T)  # random symmetric "Hessian"
        g2 = float(grad @ grad)
        H = float(np.sum(hess * hess))
        u = g2 / (l0**2 * H)
        f = u / (1 + u)
        fp = 1.0 / (1 + u) ** 2
        s = g2 / (f * a0**2)
        # analytic: dP/dpsi_ij = f'(Q - s Q') du/dpsi_ij, du/dpsi_ij = -2 u psi_ij / H
        dP_ana = fp * (Q(s) - s * dQ(s)) * (-2.0 * u * hess / H)
        # numeric: treat the 9 matrix entries as independent (P depends on them only via H)
        dP_num = np.zeros((3, 3))
        eps = 1e-6
        for i in range(3):
            for j in range(3):
                Hp = hess.copy()
                Hp[i, j] += eps
                Hm = hess.copy()
                Hm[i, j] -= eps
                dP_num[i, j] = (P_function(grad, Hp, a0, l0, Q) - P_function(grad, Hm, a0, l0, Q)) / (2 * eps)
        err = np.max(np.abs(dP_ana - dP_num)) / max(np.max(np.abs(dP_num)), 1e-12)
        max_err = max(max_err, err)
    report("dP/dpsi_ij = f'(u)[Q - sQ'] * (-2u psi_ij/H)  (analytic vs finite difference)", max_err < 1e-5,
           f"max rel. err = {max_err:.2e}")

    # u f'(u)/f(u) = 1/(1+u)
    us = np.logspace(-4, 4, 9)
    fac = us * (1.0 / (1 + us) ** 2) / screening_f(us)
    report("u f'(u)/f(u) = 1/(1+u)", np.allclose(fac, 1.0 / (1 + us), rtol=1e-12))
    small = float(1.0 / (1 + 1e4))
    report("Hessian-term factor -> 0 where f ~ 1 (u >> 1)", small < 1e-3, f"1/(1+u) at u=1e4 is {small:.1e}")
    report("Hessian-term factor -> 1 where f << 1 (u << 1)", abs(1.0 / (1 + 1e-4) - 1) < 1e-3)
    RESULTS["C2"] = dict(u=us.tolist(), factor=(1.0 / (1 + us)).tolist(), max_err=float(max_err))
    print("      u        f(u)     u f'/f = 1/(1+u)")
    for uu in us:
        print(f"      {uu:8.1e}  {float(screening_f(uu)):.4f}   {1.0/(1+uu):.4e}")


# --------------------------------------------------------------------------
# C3: QUMOND phantom density, two equal point masses
# --------------------------------------------------------------------------
def nu_simple(y):
    """'Simple' interpolating function nu(y) = 1/2 + sqrt(1/4 + 1/y)."""
    y = np.asarray(y, dtype=float)
    return 0.5 + np.sqrt(0.25 + 1.0 / y)


def nu_standard(y):
    """'Standard' interpolating function nu(y) = sqrt((1 + sqrt(1 + 4/y^2))/2)."""
    y = np.asarray(y, dtype=float)
    return np.sqrt((1.0 + np.sqrt(1.0 + 4.0 / y**2)) / 2.0)


def gN_two_masses(pos, GM=1.0, d=1.0):
    """Newtonian acceleration of two equal point masses at z = +/- d (vector, shape (...,3))."""
    pos = np.asarray(pos, dtype=float)
    out = np.zeros_like(pos)
    for sgn in (+1.0, -1.0):
        src = np.array([0.0, 0.0, sgn * d])
        dvec = src - pos
        r = np.linalg.norm(dvec, axis=-1, keepdims=True)
        out += GM * dvec / r**3
    return out


def phantom_field(pos, nu, a0, GM=1.0, d=1.0):
    """F = (nu - 1) g_N ; phantom density rho_ph = -(1/4 pi G) div F."""
    g = gN_two_masses(pos, GM, d)
    y = np.linalg.norm(g, axis=-1, keepdims=True) / a0
    return (nu(y) - 1.0) * g


def divergence_central(F, pos, h):
    """Central-difference divergence of the vector field F at points pos (shape (...,3))."""
    div = 0.0
    for k in range(3):
        e = np.zeros(3)
        e[k] = h
        div = div + (F(pos + e)[..., k] - F(pos - e)[..., k]) / (2 * h)
    return div


def check_c3(nu_name: str):
    print(f"\nC3  QUMOND phantom density, two equal masses at z=+/-d (nu = '{nu_name}', G=M=d=1, a0=0.1)")
    nu = {"simple": nu_simple, "standard": nu_standard}[nu_name]
    GM, d, a0 = 1.0, 1.0, 0.1
    F = lambda p: phantom_field(p, nu, a0, GM, d)
    four_pi_G = 4.0 * np.pi  # G = 1

    # --- (a) local sign of rho_ph on a small sphere around the null point -----
    # g_N vanishes at the midpoint, so (nu-1) ~ |g|^(-1/2) is singular there; the divergence is
    # integrable.  Sample at several radii, using h << eps.
    thetas = np.linspace(0.0, np.pi, 181)[1:-1]
    sample = {}
    for eps in (0.02, 0.05, 0.10):
        h = eps * 1e-3
        pts = np.stack([eps * np.sin(thetas), np.zeros_like(thetas), eps * np.cos(thetas)], axis=-1)
        rho_ph = -divergence_central(F, pts, h) / four_pi_G
        frac_neg_angle = float(np.mean(rho_ph < 0))
        # solid-angle weighted fraction (weight sin(theta))
        w = np.sin(thetas)
        frac_neg_solid = float(np.sum(w[rho_ph < 0]) / np.sum(w))
        sample[eps] = dict(
            rho_min=float(rho_ph.min()),
            rho_max=float(rho_ph.max()),
            frac_negative_solid_angle=frac_neg_solid,
            rho_on_axis=float(rho_ph[0]),
            rho_equator=float(rho_ph[len(rho_ph) // 2]),
        )
        print(
            f"      eps={eps:<5} rho_ph(axis)={sample[eps]['rho_on_axis']:+.3e}"
            f"  rho_ph(equator)={sample[eps]['rho_equator']:+.3e}"
            f"  solid-angle fraction with rho_ph<0: {frac_neg_solid:.2f}"
        )
    any_negative = any(v["frac_negative_solid_angle"] > 0 for v in sample.values())
    report("phantom density is NEGATIVE in part of the neighbourhood of the midpoint", any_negative)

    # --- (b) sign of the phantom mass ENCLOSED in a small sphere (flux form) ---
    # M_ph(<eps) = -(1/4 pi G) * oint (nu-1) g_N . n dA   (axisymmetric about z)
    nodes, weights = np.polynomial.legendre.leggauss(400)  # in mu = cos(theta)
    enclosed = {}
    for eps in (0.02, 0.05, 0.10):
        mu = nodes
        st = np.sqrt(1 - mu**2)
        pts = np.stack([eps * st, np.zeros_like(mu), eps * mu], axis=-1)
        n_hat = pts / eps
        Fv = F(pts)
        flux_density = np.sum(Fv * n_hat, axis=-1)
        flux = 2.0 * np.pi * eps**2 * float(np.sum(weights * flux_density))
        M_ph = -flux / four_pi_G
        # Newtonian check: oint g_N.n dA must vanish in vacuum (no mass inside)
        gN = gN_two_masses(pts, GM, d)
        flux_newton = 2.0 * np.pi * eps**2 * float(np.sum(weights * np.sum(gN * n_hat, axis=-1)))
        enclosed[eps] = dict(M_ph=M_ph, newtonian_flux=flux_newton)
        print(f"      eps={eps:<5} enclosed phantom mass M_ph = {M_ph:+.4e}   (Newtonian vacuum flux = {flux_newton:+.1e})")
    report("Newtonian flux through the sphere vanishes (vacuum check)",
           all(abs(v["newtonian_flux"]) < 1e-8 for v in enclosed.values()))
    sign_total = {k: ("negative" if v["M_ph"] < 0 else "positive") for k, v in enclosed.items()}
    print(f"      sign of enclosed phantom mass by radius: {sign_total}")

    # --- (c) coarse global map in the x-z plane: where is rho_ph negative? -----
    xs = np.linspace(0.0, 2.0, 201)[1:]
    zs = np.linspace(-2.0, 2.0, 401)
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    pts = np.stack([X, np.zeros_like(X), Z], axis=-1)
    # exclude small discs around the point masses
    mask = (np.hypot(X, Z - d) > 0.15) & (np.hypot(X, Z + d) > 0.15)
    rho = -divergence_central(F, pts, 1e-4) / four_pi_G
    rho = np.where(mask, rho, np.nan)
    neg_frac = float(np.nanmean(rho < 0))
    print(f"      x-z plane map: fraction of area (excluding r<0.15 around the masses) with rho_ph<0 = {neg_frac:.2f}; "
          f"min = {np.nanmin(rho):+.3e}, max = {np.nanmax(rho):+.3e}")
    RESULTS["C3"] = dict(nu=nu_name, near_null=sample, enclosed=enclosed, plane_negative_area_fraction=neg_frac)


# --------------------------------------------------------------------------
# C4: Gauss's law and negative enclosed mass
# --------------------------------------------------------------------------
G_SI = 6.67430e-11


def enclosed_effective_mass(r, mean_radial_accel_outward, G=G_SI):
    """Effective enclosed mass from the angle-averaged radial acceleration (outward positive).

    Gauss: <g_r>_outward * 4 pi r^2 = -4 pi G M_eff   =>   M_eff = -<g_r>_outward r^2 / G.
    Outward mean force  <=>  M_eff < 0.
    """
    return -mean_radial_accel_outward * r**2 / G


def check_c4():
    print("\nC4  Gauss's law: outward mean radial force <=> negative enclosed effective mass")
    r, M = 1.0e17, 2.0e30
    g_in = -G_SI * M / r**2  # Newtonian (inward) acceleration, outward-positive convention
    report("inward Newtonian acceleration gives +M", np.isclose(enclosed_effective_mass(r, g_in), M, rtol=1e-12))
    report("outward acceleration gives a NEGATIVE enclosed mass", enclosed_effective_mass(r, -g_in) < 0)
    # Arithmetic for the range quoted in the note (10-483 x the stellar Newtonian force).
    # Reading used here: outward force = k * G M*/r^2  =>  M_eff(stars+phantom) = -k M*,
    # i.e. phantom mass = -(k+1) M*.  This is arithmetic on the quoted range, NOT v7.06 data.
    table = []
    for k in (10, 483):
        a_out = k * G_SI * M / r**2
        M_eff = enclosed_effective_mass(r, a_out)
        table.append(dict(k=k, M_eff_over_Mstar=M_eff / M, M_phantom_over_Mstar=M_eff / M - 1.0))
        print(f"      k={k:<4} M_eff/M* = {M_eff/M:+.0f}   M_phantom/M* = {M_eff/M - 1.0:+.0f}")
    report("k x outward force => M_eff = -k M*", all(np.isclose(t["M_eff_over_Mstar"], -t["k"]) for t in table))
    RESULTS["C4"] = table


# --------------------------------------------------------------------------
# C5: dimensionless constant remark
# --------------------------------------------------------------------------
def check_c5():
    print("\nC5  Dimensionless-constant remark: l0 = X * (c^2/a0) requires X << 1")
    l_mond = C_LIGHT**2 / A0_MOND
    print(f"      c^2/a0 = {l_mond:.3e} m = {l_mond/AU:.3e} AU  (a0 = {A0_MOND:g} m/s^2)")
    rows = []
    for l0_AU in (1e3, 6e3, 1e5):
        X = l0_AU * AU / l_mond
        rows.append(dict(l0_AU=l0_AU, X=X))
        print(f"      l0 = {l0_AU:8.0e} AU  ->  X = l0 / (c^2/a0) = {X:.2e}")
    report("X is many orders of magnitude below unity for any l0 << Gpc", all(r["X"] < 1e-9 for r in rows))
    RESULTS["C5"] = dict(c2_over_a0_m=l_mond, rows=rows)


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nu", choices=["simple", "standard"], default="simple",
                    help="interpolating function used in C3 (default: simple)")
    ap.add_argument("--json", metavar="FILE", help="write machine-readable results to FILE")
    args = ap.parse_args(argv)

    print("UDCT v7.08 reproducibility checks (analytic statements only; no new physics claims)")
    check_c1()
    check_c2()
    check_c3(args.nu)
    check_c4()
    check_c5()

    print("\nSummary:", "all checks passed" if not FAILED else f"{len(FAILED)} check(s) FAILED: {FAILED}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(RESULTS, fh, indent=2, default=float)
        print(f"results written to {args.json}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
