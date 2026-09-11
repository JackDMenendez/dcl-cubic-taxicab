"""First-passage distribution of a taxispace walk onto the cubic-radius-3 shell.

Steps  = fine moves inside the sphere (the 'intermediate nodes').
Hop    = first time the walk attains L1 >= 3; the sphere then re-centers there.
Question: does the landing distribution favour the 8 body-diagonal corners
          over the 24 mixed (2,1,0) and 6 axial (3,0,0) sites?
Exact rational arithmetic, absorbing-Markov-chain solve.
"""
from fractions import Fraction as F
from itertools import product

AX   = [s for s in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))]
FD   = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s)) == 2]   # 12
BD   = [s for s in product((-1,0,1),repeat=3) if sum(map(abs,s)) == 3]   # 8

L1 = lambda p: sum(map(abs, p))
cls = lambda p: tuple(sorted(map(abs, p), reverse=True))


def solve(steps, weights=None):
    """steps: list of step vectors. weights: matching list (default uniform)."""
    if weights is None:
        weights = [F(1)] * len(steps)
    tot = sum(weights)
    prob = [F(w, 1) / tot for w in weights]

    interior = [p for p in product(range(-3, 4), repeat=3) if L1(p) <= 2]
    idx = {p: i for i, p in enumerate(interior)}
    n = len(interior)

    # absorbing landing sites encountered
    land = sorted({tuple(a + b for a, b in zip(p, s))
                   for p in interior for s in steps
                   if L1(tuple(a + b for a, b in zip(p, s))) >= 3})
    lidx = {p: j for j, p in enumerate(land)}
    m = len(land)

    # augmented system: rows = interior states, unknowns = m landing probs + 1 exp. time
    # (I - Q) X = R  , (I - Q) t = 1
    A = [[F(0)] * n for _ in range(n)]
    R = [[F(0)] * (m + 1) for _ in range(n)]
    for p in interior:
        i = idx[p]
        A[i][i] = F(1)
        R[i][m] = F(1)                                    # +1 step taken
        for s, pr in zip(steps, prob):
            q = tuple(a + b for a, b in zip(p, s))
            if L1(q) >= 3:
                R[i][lidx[q]] += pr
            else:
                A[i][idx[q]] -= pr

    # Gaussian elimination over Fractions
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        R[c], R[piv] = R[piv], R[c]
        inv = F(1) / A[c][c]
        A[c] = [v * inv for v in A[c]]
        R[c] = [v * inv for v in R[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
                R[r] = [a - f * b for a, b in zip(R[r], R[c])]

    row = R[idx[(0, 0, 0)]]
    by_class, overshoot = {}, F(0)
    for p, j in lidx.items():
        by_class[cls(p)] = by_class.get(cls(p), F(0)) + row[j]
        if L1(p) > 3:
            overshoot += row[j]
    return by_class, row[m], overshoot


def report(name, steps, weights=None):
    bc, t, over = solve(steps, weights)
    print(f"\n{name}")
    print(f"  steps per hop (exact E[T]) : {t}  = {float(t):.4f}")
    for k in sorted(bc, reverse=True):
        print(f"  land on |coords|={k!s:<10} p = {float(bc[k]):.6f}   ({bc[k]})")
    corner = bc.get((1, 1, 1), F(0))
    print(f"  CORNER (body-diagonal) share : {float(corner):.6f}")
    print(f"  uniform-over-38-sphere would give 8/38 = {8/38:.6f}")
    if over:
        print(f"  overshoot past L1=3          : {float(over):.6f}")


report("A. axial steps only (6-neighbourhood)", AX)
report("B. axial + face diagonals (18-nb)", AX + FD)
report("C. full 26-neighbourhood", AX + FD + BD)
# weighted: down-weight axial relative to face diagonals
report("D. 18-nb, face diagonals weighted 3x axial", AX + FD, [F(1)]*6 + [F(3)]*12)
report("E. 18-nb, axial weighted 3x face diagonals", AX + FD, [F(3)]*6 + [F(1)]*12)
