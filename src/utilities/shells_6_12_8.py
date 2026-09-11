"""The three shells of the 3x3x3 stencil = the three cubic Bravais lattices.

  6 axial       (+-1,0,0)      generate Z^3   (index 1)   -> SC,  coordination 6
 12 face diag   (+-1,+-1,0)    generate FCC   (index 2)   -> FCC, coordination 12
  8 body diag   (+-1,+-1,+-1)  generate BCC   (index 4)   -> BCC, coordination 8
                               6 + 12 + 8 = 26

Part 1  Kac's lemma: mean steps per hop = [Z^3 : BCC] = 4, exactly, for ANY step
        set that generates the order-4 quotient.  (For the 8 alone it is 1:
        every step is already a lattice step, so step = hop and there is no
        trembling.)
Part 2  Where the first return lands, by step set.  Note the hard zero: a
        face-diagonal walk preserves coordinate-sum parity, stays inside FCC
        forever, and can NEVER reach an all-odd body-diagonal site.
Part 3  Voronoi (Wigner-Seitz) cells: cube / rhombic dodecahedron / truncated
        octahedron, and how 6 and 8 swap between faces and vertices.
"""
from itertools import product
from collections import defaultdict
from fractions import Fraction as F

AX = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 1]   # 6
FD = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 2]   # 12
BD = [s for s in product((-1, 0, 1), repeat=3) if sum(map(abs, s)) == 3]   # 8

onBCC = lambda p: len({c & 1 for c in p}) == 1          # all coords same parity
onFCC = lambda p: sum(p) % 2 == 0                        # even coordinate sum
cls   = lambda p: tuple(sorted(map(abs, p), reverse=True))
qcos  = lambda p: ((p[0] - p[1]) % 2, (p[1] - p[2]) % 2)  # coset of BCC in Z^3

SETS = [("6  axial only ", AX), ("12 face only  ", FD), ("8  body only  ", BD),
        ("18 axial+face ", AX + FD), ("26 all shells ", AX + FD + BD)]

# ---------------- Part 1: Kac's lemma on the order-4 quotient ----------------
print("Part 1  mean steps per hop = mean return time to the identity coset")
print("        (quotient Z^3 / BCC has order 4)")
states = [(0, 0), (1, 0), (0, 1), (1, 1)]
for name, S in SETS:
    P = {a: defaultdict(F) for a in states}
    for a in states:
        for s in S:
            b = ((a[0] + qcos(s)[0]) % 2, (a[1] + qcos(s)[1]) % 2)
            P[a][b] += F(1, len(S))
    idx = {a: i for i, a in enumerate(states)}
    n, A, r = 4, [[F(0)] * 4 for _ in range(4)], [F(1)] * 4
    for a in states:
        i = idx[a]; A[i][i] += F(1)
        for b, p in P[a].items():
            if b != (0, 0): A[i][idx[b]] -= p
    for c in range(n):
        piv = next((x for x in range(c, n) if A[x][c] != 0), None)
        if piv is None: continue
        A[c], A[piv] = A[piv], A[c]; r[c], r[piv] = r[piv], r[c]
        inv = F(1) / A[c][c]; A[c] = [v * inv for v in A[c]]; r[c] *= inv
        for x in range(n):
            if x != c and A[x][c] != 0:
                f = A[x][c]
                A[x] = [u - f * v for u, v in zip(A[x], A[c])]; r[x] -= f * r[c]
    print(f"  {name}  E[steps per hop] = {r[idx[(0,0)]]}")

# ---------------- Part 2: first-return landing distribution ----------------
def first_return(steps, H=22, prune=1e-11):
    w0 = 1.0 / len(steps)
    live, land, ET = {(0, 0, 0): 1.0}, defaultdict(float), 0.0
    for t in range(1, H):
        nxt = defaultdict(float)
        for p, w in live.items():
            for s in steps:
                q = (p[0] + s[0], p[1] + s[1], p[2] + s[2])
                if onBCC(q): land[q] += w * w0; ET += w * w0 * t
                else:        nxt[q] += w * w0
        live = {k: v for k, v in nxt.items() if v > prune}
        if sum(live.values()) < 1e-12: break
    byc = defaultdict(float)
    for p, v in land.items(): byc[cls(p)] += v
    return byc, sum(byc.values())

print("\nPart 2  first-return landing distribution")
print(f"  {'step set':<15} {'null':>8} {'6-type':>8} {'8-type':>8} "
      f"{'12-type':>8} {'other':>8} {'8/(6+8)':>9}")
for name, S in SETS:
    byc, tot = first_return(S)
    h, a6, a8, a12 = byc[(0,0,0)], byc[(2,0,0)], byc[(1,1,1)], byc[(2,2,0)]
    other = tot - h - a6 - a8 - a12
    sh = f"{a8/(a6+a8):9.3f}" if a6 + a8 else f"{'--':>9}"
    print(f"  {name:<15} {h:8.4f} {a6:8.4f} {a8:8.4f} {a12:8.4f} {other:8.4f} {sh}")
print("  (the 12-only zero on 8-type is structural: even-sum steps can never"
      " reach an all-odd site)")

# ---------------- Part 3: Voronoi cells ----------------
try:
    import numpy as np
    from scipy.optimize import linprog
    def cell(name, member, cap):
        pts = [p for p in product(range(-2, 3), repeat=3)
               if p != (0, 0, 0) and member(p) and sum(c*c for c in p) <= cap]
        V = np.array(pts, float); b = (V ** 2).sum(1) / 2; keep = []
        for i in range(len(pts)):
            r = linprog(-V[i], A_ub=np.delete(V, i, 0), b_ub=np.delete(b, i),
                        bounds=[(-12, 12)] * 3)
            if r.success and -r.fun > b[i] + 1e-9: keep.append(i)
        H, hb, n, vs = V[keep], b[keep], len(keep), []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    M = H[[i, j, k]]
                    if abs(np.linalg.det(M)) < 1e-9: continue
                    x = np.linalg.solve(M, hb[[i, j, k]])
                    if np.all(H @ x <= hb + 1e-7): vs.append(x)
        vs = np.unique(np.round(np.array(vs), 6), axis=0)
        fs, vd = defaultdict(int), defaultdict(int)
        for i in range(n): fs[int(np.sum(np.abs(vs @ H[i] - hb[i]) < 1e-6))] += 1
        for v in vs:       vd[int(np.sum(np.abs(H @ v - hb) < 1e-6))] += 1
        print(f"  {name:<5} faces={n:<3} vertices={len(vs):<3} "
              f"face sizes {dict(sorted(fs.items()))}  vertex degrees {dict(sorted(vd.items()))}")
    print("\nPart 3  Voronoi (Wigner-Seitz) cells")
    cell("SC",  lambda p: True, 4)
    cell("FCC", onFCC, 8)
    cell("BCC", onBCC, 8)
except ImportError:
    print("\nPart 3 skipped (needs numpy + scipy)")
