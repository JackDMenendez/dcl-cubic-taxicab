"""How many axes are 'enough'? Direction counting vs moment isotropy.

Worry: 6/8/12 discrete directions cannot represent a continuum of directions.
Answer: directional resolution is not set by coordination number -- it grows with
how many ticks you spend -- and isotropy is not set by how MANY directions there
are but by which spherical-harmonic moments vanish.
"""
from itertools import product
from math import gcd, pi, sqrt
from fractions import Fraction as F
from collections import defaultdict
import numpy as np

AX = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s))==1]
FD = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s))==2]
BD = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s))==3]

# ---------- 1. how many distinct directions exist within radius R ----------
print("1. DISTINCT LATTICE DIRECTIONS within radius R (primitive vectors, gcd=1)")
print(f"   {'R':>4} {'directions':>12} {'ang. spacing':>14}")
for R in (1, 2, 3, 5, 10, 20, 40):
    dirs = {tuple(np.array(p)//gcd(gcd(abs(p[0]),abs(p[1])),abs(p[2])))
            for p in product(range(-R,R+1),repeat=3)
            if p != (0,0,0) and sum(c*c for c in p) <= R*R}
    n = len(dirs)
    print(f"   {R:>4} {n:>12} {sqrt(4*pi/n):>13.4f} rad  = {sqrt(4*pi/n)*180/pi:6.2f} deg")
print("   -> grows as R^3, so angular spacing falls as R^-1.5.  Coordination number")
print("      limits directions per TICK, not directions available to the process.")

# ---------- 2. moment isotropy of each shell ----------
def moments(vecs, weights=None):
    """Return <n_i n_j> deviation and the degree-4 cubic invariant <n1^4+n2^4+n3^4>."""
    V = np.array(vecs, float)
    N = V / np.linalg.norm(V, axis=1, keepdims=True)
    w = np.ones(len(N)) if weights is None else np.asarray(weights, float)
    w = w / w.sum()
    C = np.einsum('a,ai,aj->ij', w, N, N)
    off = np.abs(C - np.diag(np.diag(C))).max()
    diagspread = np.ptp(np.diag(C))
    cub = float(np.einsum('a,a->', w, (N**4).sum(1)))
    return off, diagspread, cub

ISO4 = 3/5   # <n1^4+n2^4+n3^4> for the uniform sphere
print("\n2. MOMENT ISOTROPY per shell  (unit vectors)")
print(f"   {'shell':<9}{'off-diag <nn>':>15}{'diag spread':>13}"
      f"{'<n1^4+n2^4+n3^4>':>19}{'dev from 3/5':>14}")
for nm, S in (("6 axial", AX), ("12 face", FD), ("8 body", BD)):
    off, ds, cub = moments(S)
    print(f"   {nm:<9}{off:>15.2e}{ds:>13.2e}{cub:>19.6f}{cub-ISO4:>+14.6f}")
print(f"   isotropic target {ISO4}")
print("   -> degree 2 is EXACT for all three (this is why G(x) is spherical).")
print("      degree 4 fails, and the 6-shell deviates with the OPPOSITE SIGN to 8 and 12.")

# ---------- 3. mixtures that cancel the degree-4 anisotropy ----------
print("\n3. MIXTURES that kill the degree-4 term exactly")
c6, c12, c8 = 1.0, 0.5, 1/3
w = F(2,5)
print(f"   axial + body :  w6 = {F(2,5)}  w8 = {F(3,5)}"
      f"   check = {float(F(2,5))*c6 + float(F(3,5))*c8:.6f}")
# axial + face
w6 = (ISO4 - c12)/(c6 - c12)
print(f"   axial + face :  w6 = {w6:.6f} = {F(w6).limit_denominator(50)}"
      f"   w12 = {1-w6:.6f}   check = {w6*c6+(1-w6)*c12:.6f}")
print("   -> note 2/5 : 3/5.  The SAME rationals as the hop parity split of section 5.5,")
print("      but with the shells swapped (there: 2/5 was the BODY-diagonal family).")

# ---------- 4. the actual hop distribution's degree-4 anisotropy ----------
onBCC = lambda p: len({c & 1 for c in p}) == 1
live, land = {(0,0,0): 1.0}, defaultdict(float)
for t in range(1, 300):
    nxt = defaultdict(float)
    for p, wt in live.items():
        for s in AX:
            q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
            if onBCC(q): land[q] += wt/6
            else:        nxt[q] += wt/6
    live = {k: v for k, v in nxt.items() if v > 1e-18}
    if sum(live.values()) < 1e-15: break
sites = [(p, w) for p, w in land.items() if p != (0,0,0)]
V = [p for p, _ in sites]; W = [w for _, w in sites]
off, ds, cub = moments(V, W)
print("\n4. THE ACTUAL HOP DISTRIBUTION (axial walk, first return, null hop excluded)")
print(f"   off-diagonal <n_i n_j> = {off:.2e}   diagonal spread = {ds:.2e}")
print(f"   <n1^4+n2^4+n3^4> = {cub:.6f}   target {ISO4}   deviation {cub-ISO4:+.6f}")
print(f"   relative degree-4 anisotropy = {(cub-ISO4)/ISO4:+.4%}")
print("   -> degree 2 exact; degree 4 off by a few percent, NOT the tens of percent")
print("      of any single shell.  The hop mixture is already partly self-cancelling.")
