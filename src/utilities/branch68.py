"""The six ARE half-steps: which lattice site does the walk return to first?

BCC = sites of Z^3 whose coordinates all share one parity.
  the 8 diagonals (+-1,+-1,+-1) ARE lattice vectors, |v| = sqrt(3)
  the 6 axials    (+-1,0,0)     are NOT -- they are half of (+-2,0,0), |v| = 2

So a unit-axial walk spends its time on non-lattice ("virtual") sites and only
intermittently touches BCC.  At its FIRST return to the lattice, does it land on
a 6-type axial neighbour, an 8-type diagonal neighbour, or back home?

Forward propagation of exact rational mass; absorption is geometric so a modest
horizon captures essentially all of it.
"""
from fractions import Fraction as F
from collections import defaultdict

AX = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC = lambda p: len({c & 1 for c in p}) == 1
cls   = lambda p: tuple(sorted(map(abs,p), reverse=True))
P6    = F(1,6)

HORIZON = 60
live = {(0,0,0): F(1)}          # current mass on non-lattice sites (origin at t=0)
landed = defaultdict(F)         # first-return landing site -> mass
Etime  = F(0)
mass_by_t = []

for t in range(1, HORIZON+1):
    nxt = defaultdict(F)
    for p, w in live.items():
        for s in AX:
            q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
            if onBCC(q):
                landed[q] += w*P6
                Etime     += w*P6*t
            else:
                nxt[q] += w*P6
    live = nxt
    mass_by_t.append((t, sum(landed.values())))
    if not live:
        break

resid = sum(live.values())
byc = defaultdict(F)
for p, w in landed.items():
    byc[cls(p)] += w

print(f"horizon {HORIZON} steps   absorbed mass = {float(sum(landed.values())):.12f}"
      f"   unabsorbed = {float(resid):.3e}")
print(f"E[steps to first lattice return] = {float(Etime):.6f}")
print("\nfirst-return landing site, by class:")
for k in sorted(byc, key=lambda k: -float(byc[k])):
    if float(byc[k]) < 1e-5: continue
    tag = {(0,0,0): "HOME (return)",
           (2,0,0): "6-type AXIAL  |v|=2",
           (1,1,1): "8-type DIAGONAL |v|=sqrt(3)"}.get(k, "")
    print(f"  {str(k):<12} p = {float(byc[k]):.8f}   {tag}")

home = byc[(0,0,0)]
a6, a8 = byc[(2,0,0)], byc[(1,1,1)]
print(f"\n  HOME    {float(home):.8f}   = {home}")
print(f"  6-type  {float(a6):.8f}   = {a6}")
print(f"  8-type  {float(a8):.8f}   = {a8}")
print(f"  moved-away total = {float(a6+a8):.8f}")
print(f"  6 : 8  ratio = {float(a6/a8):.8f}      8/(6+8) = {float(a8/(a6+a8)):.8f}")
print(f"  per-site: axial {float(a6/6):.8f}   diagonal {float(a8/8):.8f}"
      f"   ratio {float((a8/8)/(a6/6)):.6f}")
print("\nabsorbed mass by step:", [(t, round(float(m),6)) for t,m in mass_by_t[:8]])
