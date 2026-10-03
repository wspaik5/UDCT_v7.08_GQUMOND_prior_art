# UDCT v7.08 — reproducibility checks

Companion code for **UDCT v7.08, "Prior-art correction: the action-derived Hessian source term and
screening by second derivatives of the Newtonian potential belong to Milgrom's generalized QUMOND
(GQUMOND) class"** (Won Shik Paik, 30 September 2026).

## Read this first: what this repository does and does not do

UDCT v7.08 is a **correction note with no new calculation**. There is therefore no new numerical
result to "reproduce". This script instead **independently re-derives or numerically checks the
elementary analytic statements the note relies on**, so a reader can confirm them without trusting
the prose. It is a sanity-check suite, not a new test of any theory.

**Checked here (C1–C5)**

| ID | Statement in the note | What the script does |
|----|------------------------|----------------------|
| C1 | Milgrom's example screening variable `u = (∇ψ)² / (ℓ0² ψ,ij ψ,ij)` with `f(u) = u/(1+u)` screens small systems | For a point mass, verifies by finite differences that `u = r²/(6ℓ0²)` (so `ψ,ijψ,ij = 6G²M²/r⁶`), and the limits `f→0` (u≪1), `f→1` (u≫1) |
| C2 | For `P = f(u) Q[(∇ψ)²/(f a0²)]` the Hessian term "scales with the local logarithmic derivative of f and is small where f ≈ 1" | Verifies `∂P/∂ψ,ij = f′(u)[Q − sQ′]·(−2u ψ,ij/H)` against finite differences for random symmetric Hessians and a test `Q`; shows the relative factor `u f′/f = 1/(1+u) → 0` for `u≫1` |
| C3 | "In QUMOND the phantom density is negative in regions such as the neighbourhood of the midpoint between two masses" | Computes `ρ_ph = −(1/4πG) ∇·[(ν−1) g_N]` for two equal point masses; reports local sign near the null point and the sign of the phantom mass *enclosed* by small spheres around it |
| C4 | An outward angle-averaged radial force means, by Gauss's law, a negative enclosed effective mass | Unit checks of `M_eff = −⟨g_r⟩_out r²/G`; arithmetic for the "10–483×" range quoted in the note |
| C5 | The extra dimensioned constant can be replaced by a dimensionless multiple of `c²/a0`, "necessarily very different from unity" | Computes `X = ℓ0 / (c²/a0)` for illustrative `ℓ0` values |

**Not checked here** (needs the original sources or the v7.06 released script)

- The exact form and numbering of Milgrom (2023) Eqs. 7 and 14 — read arXiv:2305.01589 directly.
- The v7.06 Table 1 values (outward force of 10–483× the stellar Newtonian force). Those come from
  the v7.06 script and data, which are not part of this repository. C4 is **arithmetic on the quoted
  range only**, under the reading "outward force = k × GM*/r² ⇒ M_eff = −k M*".
- The "≈ 6000 AU" solar balance scale quoted in the note.
- **Whether QUMOND or GQUMOND can make the enclosed effective mass negative in a strong external
  field.** The note explicitly defers this to v7.09. Nothing here answers it; see the C3 caveat below.
- Whether Milgrom's explicit `P = fQ` example shares the v7.06 obstruction (open question, v7.09).

## Requirements and usage

- Python ≥ 3.9, `numpy` only. (Tested: Python 3.12.3, numpy 2.4.4.)
- Deterministic (fixed RNG seed in C2); runs in a few seconds.

```bash
python udct_v7_08_checks.py                  # all checks, "simple" interpolating function
python udct_v7_08_checks.py --nu standard    # C3 with the "standard" interpolating function
python udct_v7_08_checks.py --json out.json  # also write machine-readable results
```

Exit code is `0` if every check passes, `1` otherwise.

## Expected output (default run)

The Newtonian vacuum-flux digits in C3 are roundoff and will differ between machines at the
`1e-17`–`1e-15` level. The check only requires `|flux| < 1e-8`. All other printed values below
matched a run on Python 3.12 / numpy 2.x.

```text
UDCT v7.08 reproducibility checks (analytic statements only; no new physics claims)

C1  Point-mass Hessian identity and screening-function limits
      r/l0=0.3   u_num=0.015000  u_ana=0.015000  f(u)=0.0148
      r/l0=1.0   u_num=0.166667  u_ana=0.166667  f(u)=0.1429
      r/l0=3.0   u_num=1.500000  u_ana=1.500000  f(u)=0.6000
  [PASS] u = r^2/(6 l0^2) for a point mass (|psi_ij psi_ij| = 6 G^2M^2/r^6) -- max rel. err = 2.03e-07
  [PASS] f(u) -> 0 for u << 1 (screening ON in small systems)
  [PASS] f(u) -> 1 for u >> 1 (standard QUMOND recovered)
  [PASS] u = 1 at r = sqrt(6) l0 for a point mass

C2  Chain rule for P = f(u) Q[(grad psi)^2/(f a0^2)] and the factor u f'/f
  [PASS] dP/dpsi_ij = f'(u)[Q - sQ'] * (-2u psi_ij/H)  (analytic vs finite difference) -- max rel. err = 3.82e-09
  [PASS] u f'(u)/f(u) = 1/(1+u)
  [PASS] Hessian-term factor -> 0 where f ~ 1 (u >> 1) -- 1/(1+u) at u=1e4 is 1.0e-04
  [PASS] Hessian-term factor -> 1 where f << 1 (u << 1)
      u        f(u)     u f'/f = 1/(1+u)
       1.0e-04  0.0001   9.9990e-01
       1.0e-03  0.0010   9.9900e-01
       1.0e-02  0.0099   9.9010e-01
       1.0e-01  0.0909   9.0909e-01
       1.0e+00  0.5000   5.0000e-01
       1.0e+01  0.9091   9.0909e-02
       1.0e+02  0.9901   9.9010e-03
       1.0e+03  0.9990   9.9900e-04
       1.0e+04  0.9999   9.9990e-05

C3  QUMOND phantom density, two equal masses at z=+/-d (nu = 'simple', G=M=d=1, a0=0.1)
      eps=0.02  rho_ph(axis)=+1.627e-01  rho_ph(equator)=-1.198e-01  solid-angle fraction with rho_ph<0: 0.33
      eps=0.05  rho_ph(axis)=+9.296e-02  rho_ph(equator)=-7.054e-02  solid-angle fraction with rho_ph<0: 0.33
      eps=0.1   rho_ph(axis)=+5.883e-02  rho_ph(equator)=-4.436e-02  solid-angle fraction with rho_ph<0: 0.33
  [PASS] phantom density is NEGATIVE in part of the neighbourhood of the midpoint
      eps=0.02  enclosed phantom mass M_ph = +1.9527e-06   (Newtonian vacuum flux = roundoff, |flux| << 1e-8)
      eps=0.05  enclosed phantom mass M_ph = +1.7912e-05   (Newtonian vacuum flux = roundoff, |flux| << 1e-8)
      eps=0.1   enclosed phantom mass M_ph = +9.3419e-05   (Newtonian vacuum flux = roundoff, |flux| << 1e-8)
  [PASS] Newtonian flux through the sphere vanishes (vacuum check)
      sign of enclosed phantom mass by radius: {0.02: 'positive', 0.05: 'positive', 0.1: 'positive'}
      x-z plane map: fraction of area (excluding r<0.15 around the masses) with rho_ph<0 = 0.01; min = -1.735e-01, max = +1.581e-01

C4  Gauss's law: outward mean radial force <=> negative enclosed effective mass
  [PASS] inward Newtonian acceleration gives +M
  [PASS] outward acceleration gives a NEGATIVE enclosed mass
      k=10   M_eff/M* = -10   M_phantom/M* = -11
      k=483  M_eff/M* = -483   M_phantom/M* = -484
  [PASS] k x outward force => M_eff = -k M*

C5  Dimensionless-constant remark: l0 = X * (c^2/a0) requires X << 1
      c^2/a0 = 7.490e+26 m = 5.007e+15 AU  (a0 = 1.2e-10 m/s^2)
      l0 =    1e+03 AU  ->  X = l0 / (c^2/a0) = 2.00e-13
      l0 =    6e+03 AU  ->  X = l0 / (c^2/a0) = 1.20e-12
      l0 =    1e+05 AU  ->  X = l0 / (c^2/a0) = 2.00e-11
  [PASS] X is many orders of magnitude below unity for any l0 << Gpc

Summary: all checks passed
```

## How to read the C3 result (important caveat)

In the two-mass QUMOND example (units `G = M = d = 1`, `a0 = 0.1`, two masses at `z = ±d`):

- The local phantom density near the midpoint is **angle dependent**: positive along the axis joining
  the masses and negative in the equatorial region (about one third of the solid angle on small
  spheres). This supports the note's statement that negative phantom density occurs in *local
  regions* near the midpoint, for both interpolating functions tested.
- The phantom mass **enclosed** by small spheres centred on the midpoint is **positive** in this
  example (≈ 2×10⁻⁶ to 9×10⁻⁵ in these units for ε = 0.02–0.1). So this configuration does **not**
  make the enclosed effective mass negative.
- This is consistent with the note's distinction between *local negative lobes* (known in MOND) and
  an *enclosed negative effective mass* inside a body (the v7.05–v7.06 result). It does **not**
  settle whether the latter can arise in QUMOND/GQUMOND in a strong external field — the v7.09
  diagnostic. Do not cite C3 as evidence for or against that question.
- The two-mass setup, `a0 = 0.1` and the two interpolating functions are choices made for this
  script; they are not taken from the note.

## Files

- `UDCT_v7.08_Prior_Art_Correction_GQUMOND_Hessian_Term_2026-09-30.pdf` — the correction note (CC BY 4.0)
- `udct_v7_08_checks.py` — the checks (self-contained, documented in its module docstring; MIT)
- `README.md` — this file (MIT)
- `LICENSE` — MIT license for the script and this README

## Provenance and acknowledgement of assistance

This script and README were written with the assistance of a large language model (Anthropic Claude)
under the author's direction, and the numerical output above was produced by running the script.
The note being checked (v7.08) was likewise drafted with LLM assistance, as stated in its
Acknowledgments. Independent verification by a physicist is welcome; please report errors via the
repository issue tracker.

## Reference

M. Milgrom, *Generalizations of quasilinear MOND*, Phys. Rev. D **108**, 084005 (2023),
arXiv:2305.01589.

## License

The note (PDF) is licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0). The check script and this README are licensed under the MIT License (see `LICENSE`).

Copyright (c) 2026 Won Shik Paik.
