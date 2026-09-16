"""Redo the core results on the ACTUAL geometry: alternating tetrahedra (T3_diamond).

Green  V1..V4  = (1,1,1) (1,-1,-1) (-1,1,-1) (-1,-1,1)   sign-product +1
Red   -V1..-V4                                            sign-product -1
Green is used from one sublattice, red from the other.

Site set: all-even and all-odd coordinate triples (= BCC), bonded tetrahedrally.
Coordination 4, bipartite, and this is the diamond lattice.

The key reduction proved below: the TWO-STEP process, rescaled by 2, is a LAZY WALK
ON THE 12 FACE DIAGONALS with laziness exactly 1/4.  Everything else follows.
"""
import numpy as np
from itertools import product
from collections import defaultdict
from fractions import Fraction as F

GREEN = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
RED   = [tuple(-c for c in v) for v in GREEN]
FD    = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s)) == 2]   # 12

# ---------- 1. the two-step reduction ----------
print("1. TWO-STEP REDUCTION")
pairs = defaultdict(F)
for a in GREEN:
    for b in RED:
        d = tuple(x+y for x, y in zip(a, b))
        pairs[d] += F(1, 16)
null = pairs[(0,0,0)]
moves = {d: w for d, w in pairs.items() if d != (0,0,0)}
halved = sorted({tuple(c//2 for c in d) for d in moves})
print(f"   P(null per two-step)      = {null}  = {float(null)}")
print(f"   distinct non-null moves   = {len(moves)}, each with weight {list(moves.values())[0]}")
print(f"   halved, they are the 12 face diagonals: {sorted(halved) == sorted(FD)}")
print(f"   => two-step, rescaled by 2: stay w.p. 1/4, else uniform on the 12 FCC vectors")
print()

# ---------- 2. Green's function on the rescaled sublattice ----------
print("2. GREEN'S FUNCTION AND THE FIELD")
L = 192
n = np.fft.fftfreq(L, 1.0/L); k = 2*np.pi*n/L; c = np.cos(k)
c1 = c[:, None, None]; c2 = c[None, :, None]; c3 = c[None, None, :]
phi_fcc = (c1*c2 + c2*c3 + c3*c1)/3.0
phi2 = 0.25 + 0.75*phi_fcc                 # the lazy-FCC structure function
den = 1.0 - phi2
# BUGFIX (verified 2026-09-15): the FCC structure function equals 1 at TWO points,
# k = 0 AND k = (pi,pi,pi) -- the index-2 sublattice (sec 9.11 records "12 face -> zeros 2").
# Masking only the origin left a division by zero and made all of section 2 inf/nan.
den[0,0,0] = np.inf
den[L//2, L//2, L//2] = np.inf
Ginv = 1.0/den
vol = L**3
S = Ginv.sum(axis=(1,2))                   # separable partial sum for axis direction

print(f"   G(0) on the rescaled sublattice = {Ginv.sum()/vol:.6f}")
print(f"   {'r*':>4} {'G':>12} {'r*G':>10}   (r* in rescaled units = 2 lattice units)")
Gax = {}
for r in range(0, 70):
    Gax[r] = float(np.cos(k*r) @ S)/vol
for r in (2, 4, 8, 12, 16, 24):
    print(f"   {r:>4} {Gax[r]:12.8f} {r*Gax[r]:10.6f}")
print(f"   field |grad G|, r^2*E  (should converge to the Coulomb coefficient):")
# BUGFIX (verified 2026-09-15): the walk lives on the EVEN-SUM sublattice, so G is ~0
# on odd-sum sites (measured 4e-16).  A unit-spacing gradient samples empty sites and
# returns ~0.  Use spacing 2.
for r in (8, 12, 16, 24, 32, 40):
    E = abs(Gax[r-2] - Gax[r+2])/4
    print(f"     r*={r:>3}  E={E:.8f}  r^2 E = {r*r*E:.6f}")
# CORRECTION (verified 2026-09-15): the previous two prediction lines were wrong.
#   "1/pi"        omitted the factor 2 for the index-2 sublattice (half the site density).
#   "halves to 1/(2pi)"  was wrong in DIRECTION (r = 2r* means C = 2C*, it doubles)
#                        and in base -- off by 8x against the correct value.
print(f"   naive 3/(2 pi <|s|^2>), <|s|^2> = (3/4)*2 = 1.5      : {3/(2*np.pi*1.5):.6f} = 1/pi")
print(f"   x2 for the index-2 sublattice (measured limit)       : {2/np.pi:.6f} = 2/pi")
print(f"   in ORIGINAL lattice units (r = 2 r*), C = 2 C*       : {4/np.pi:.6f} = 4/pi")
print(f"   ratio to the 6-axial 3/(2pi) = {3/(2*np.pi):.6f}        : "
      f"{(4/np.pi)/(3/(2*np.pi)):.6f} = 8/3")
print()

# ---------- 3. band structure of the diamond adjacency ----------
print("3. BAND STRUCTURE  (bipartite: E = +- |f(k)|)")
f = 2*(np.exp(1j*k)[:,None,None]*np.cos(k[None,:,None]+k[None,None,:])
       + np.exp(-1j*k)[:,None,None]*np.cos(k[None,:,None]-k[None,None,:]))
absf = np.abs(f)
print(f"   |f| range      : {absf.min():.6f} .. {absf.max():.6f}   (coordination 4)")
print(f"   zero-point sum : <|f|>/2 = {absf.mean()/2:.6f} per site")
print(f"   fraction of BZ with |f| < 0.05 : {(absf<0.05).mean():.6f}  (nodal LINE, codim 2)")
print()

# ---------- 4. moments ----------
print("4. MOMENTS of the tetrahedral step set")
V = np.array(GREEN, float); Nv = V/np.linalg.norm(V, axis=1, keepdims=True)
cov = Nv.T @ Nv / 4
cub = float((Nv**4).sum(1).mean())
print(f"   covariance = {np.diag(cov).round(6)}, max off-diagonal "
      f"{np.abs(cov-np.diag(np.diag(cov))).max():.2e}  -> isotropic")
print(f"   <n1^4+n2^4+n3^4> = {cub:.6f}   target 3/5   deviation {cub-0.6:+.6f}"
      f"   relative {(cub-0.6)/0.6:+.2%}")
print(f"   vector sum = {tuple(V.sum(0).astype(int))}  -> no drift from a full tetrahedron")
