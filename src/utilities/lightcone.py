"""Light-cone shape and Lorentz anisotropy of the taxispace walk.

Two limits behave completely differently:

  DIFFUSIVE  (the field, the Green's function)  -- exactly isotropic.
             Every shell's step covariance is a multiple of delta_ij, so the CLT
             erases the lattice.  This is why G(x) came out spherically symmetric.

  BALLISTIC  (the light cone, massless propagation)  -- anisotropic by sqrt(3).
             Max speed in direction n is 1/||n||_1: 1 along an axis, 1/sqrt(2)
             along a face diagonal, 1/sqrt(3) along a body diagonal.

A massless particle is the ballistic case, so it inherits the anisotropy.

CORRECTION (10 Sept 2026).  An earlier version of this docstring concluded:
"Terrestrial Lorentz tests bound c-anisotropy at the 1e-18 level; sqrt(3) is a 73%
effect.  This is a falsification-level problem, not a calibration issue."
THAT WAS WRONG and is retracted.  It compared the single-tick reachable-set
anisotropy against a bound derived from long-wavelength optics -- different
quantities.  An experiment measures a wave with lambda >> a, governed by the
DISPERSION RELATION, not by the cone's support; the observable anisotropy is
suppressed as (ka)^2 (~ k^2/18 for the axial quadratic form).  The experiment is
therefore a CONSTRAINT ON LATTICE SPACING (a < 0.34 fm at 1e-18), which any
Planck- or nuclear-scale lattice satisfies with room to spare -- a derived bound,
not a refutation.  See section 9.4 of the notes.

The numbers below are unaffected: the sqrt(3) single-tick cone anisotropy is real,
and the ballistic table shows it is exactly scale-invariant in t (it never decays),
so the (ka)^2 dispersion suppression is doing all the work.
"""
from itertools import product
import numpy as np

AX = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 1]   # 6
FD = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 2]   # 12
BD = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 3]   # 8
SETS = [("6 axial", AX), ("12 face", FD), ("8 body", BD), ("26 all", AX + FD + BD)]

print("BALLISTIC light cone: max Euclidean displacement per tick")
for label, S in SETS:
    reach = {(0, 0, 0)}
    vals = []
    for t in range(1, 7):
        reach = {tuple(a + b for a, b in zip(p, s)) for p in reach for s in S}
        vals.append(max(np.linalg.norm(p) for p in reach) / t)
    print(f"  {label:<9} " + "  ".join(f"t={t}:{v:.6f}" for t, v in enumerate(vals, 1)))

print("\nDirection-dependent max speed,  v(n) = 1 / ||n||_1")
for nm, v in (("axis      (1,0,0)", (1,0,0)),
              ("face diag (1,1,0)", (1,1,0)),
              ("body diag (1,1,1)", (1,1,1))):
    n = np.array(v, float); n /= np.linalg.norm(n)
    l1 = np.abs(n).sum()
    print(f"  {nm}   ||n||_1 = {l1:.6f}   v = {1/l1:.6f}")
print(f"  anisotropy fastest/slowest = {np.sqrt(3):.6f} = sqrt(3)   (a 73% effect)")

print("\nOne-tick reachable set (the cone shape):")
print("  6 axial -> convex hull {+-e_i}         = OCTAHEDRON  (axes fastest)")
print("  8 body  -> convex hull {(+-1,+-1,+-1)} = CUBE        (diagonals fastest)")
print("  26 all  -> hull of both                = CUBE = the L-inf ball")
print("  (i.e. the causal cone of the full stencil is the original 'cubic sphere',")
print("   and its 8 corners are the maximal-speed directions at sqrt(3) per tick)")

print("\nDIFFUSIVE limit: step covariance is exactly isotropic")
for label, S in SETS:
    A = np.array(S, float)
    C = A.T @ A / len(A)
    off = np.abs(C - np.diag(np.diag(C))).max()
    print(f"  {label:<9} cov = {np.diag(C).round(6)} * I    max off-diagonal = {off:.2e}")
print("  => the field inherits no anisotropy; only the light cone does.")
