"""Confine the walk to a finite cell with a Markov chain; hop on BCC arrival.

Motivation: the unconfined walk sends 16-26% of hops beyond the 14-neighbourhood
(the 'tail'), which no truncation rule removes cleanly.  Confining the walk to a
box removes the tail BY CONSTRUCTION -- the hop can only land on a BCC site inside
the box.

Key structural fact.  Inside the 3x3x3 cell (L-inf <= 1) the BCC sites are exactly
    the centre (0,0,0)          -> null hop
    the 8 corners (+-1,+-1,+-1) -> body-diagonal hop
The 6 axial and 12 face sites are NOT BCC.  So confining to 27 nodes automatically
yields 'null + 8 axes' with no 6-type, no 12-type and no tail.

Inside the 5x5x5 cell (L-inf <= 2) the BCC sites additionally include (+-2,0,0),
(+-2,+-2,0), (+-2,+-2,+-2), restoring the shell mixture with a bounded tail.

Boundary conventions compared:
  stay   -- an illegal step leaves the walker in place and consumes a tick
  resample -- an illegal step is redrawn; only legal moves, no tick wasted
  reflect  -- an illegal step is mirrored back into the box
"""
from itertools import product
from fractions import Fraction as F
from collections import defaultdict
import numpy as np

AX = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC = lambda p: len({c & 1 for c in p}) == 1
cls   = lambda p: tuple(sorted(map(abs,p), reverse=True))
Linf  = lambda p: max(map(abs,p))


def step_targets(p, R, mode):
    """Yield (target, weight) pairs summing to 1."""
    legal, illegal = [], []
    for s in AX:
        q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
        (legal if Linf(q) <= R else illegal).append((s, q))
    if mode == "resample":
        w = F(1, len(legal))
        return [(q, w) for _, q in legal]
    out = [(q, F(1,6)) for _, q in legal]
    if mode == "stay":
        if illegal: out.append((p, F(len(illegal), 6)))
    elif mode == "reflect":
        # A unit step that would leave the box bounces straight back: take its
        # negation instead.  (For unit steps in a box this is always legal.)
        # NOTE an earlier version of this file mirrored the coordinate as c-2d,
        # which is a 2-unit point reflection, not a wall bounce.  That was a bug
        # and it made 'reflect' look like the worst convention; it is the best.
        for s, _ in illegal:
            r = (p[0]-s[0], p[1]-s[1], p[2]-s[2])
            assert Linf(r) <= R
            out.append((r, F(1,6)))
    return out


def solve(R, mode):
    trans = [p for p in product(range(-R, R+1), repeat=3)
             if Linf(p) <= R and not onBCC(p)]
    idx = {p: i for i, p in enumerate(trans)}
    n = len(trans)
    absorb = sorted({q for p in trans for q, _ in step_targets(p, R, mode) if onBCC(q)})
    lidx = {p: j for j, p in enumerate(absorb)}
    m = len(absorb)

    A = [[F(0)]*n for _ in range(n)]
    Rh = [[F(0)]*(m+1) for _ in range(n)]
    for p in trans:
        i = idx[p]; A[i][i] += F(1); Rh[i][m] = F(1)
        for q, w in step_targets(p, R, mode):
            if onBCC(q): Rh[i][lidx[q]] += w
            else:        A[i][idx[q]] -= w
    for c in range(n):
        piv = next(x for x in range(c, n) if A[x][c] != 0)
        A[c], A[piv] = A[piv], A[c]; Rh[c], Rh[piv] = Rh[piv], Rh[c]
        inv = F(1)/A[c][c]
        A[c] = [v*inv for v in A[c]]; Rh[c] = [v*inv for v in Rh[c]]
        for x in range(n):
            if x != c and A[x][c] != 0:
                f = A[x][c]
                A[x] = [u-f*v for u, v in zip(A[x], A[c])]
                Rh[x] = [u-f*v for u, v in zip(Rh[x], Rh[c])]

    out = [F(0)]*(m+1)
    for q, w in step_targets((0,0,0), R, mode):     # one step off the centre
        if onBCC(q): out[lidx[q]] += w
        else:
            row = Rh[idx[q]]
            for j in range(m+1): out[j] += w*row[j]
    out[m] += 1
    byc = defaultdict(F)
    for p, j in lidx.items(): byc[cls(p)] += out[j]
    return byc, out[m], n


ISO4 = 3/5
def cubic4(byc):
    num = den = 0.0
    for k, v in byc.items():
        if k == (0,0,0): continue
        nv = np.array(k, float); nv /= np.linalg.norm(nv)
        num += float(v) * (nv**4).sum(); den += float(v)
    return num/den if den else float('nan')


for R, name in ((1, "3x3x3 cell (27 nodes)"), (2, "5x5x5 cell (125 nodes)")):
    print(f"\n=== {name} ===")
    for mode in ("stay", "resample", "reflect"):
        byc, ET, n = solve(R, mode)
        tot = sum(byc.values())
        print(f"  boundary = {mode:<9} transient states = {n:<3} "
              f"E[steps/hop] = {str(ET):>12} = {float(ET):.6f}")
        for k in sorted(byc, key=lambda k: -float(byc[k])):
            if float(byc[k]) < 1e-6: continue
            tag = {(0,0,0):"null", (1,1,1):"body 8", (2,0,0):"axial 6",
                   (2,2,0):"face 12", (2,2,2):"body 8 doubled"}.get(k, "")
            print(f"      {str(k):<12} {float(byc[k]):.8f}  {str(byc[k]):>22}  {tag}")
        c4 = cubic4(byc)
        print(f"      total = {float(tot):.8f}   degree-4 <n^4 sum> = {c4:.6f}"
              f"   dev {c4-ISO4:+.6f}   (target {ISO4})")
