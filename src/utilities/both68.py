"""Is the '6-axis' and '8-axis' picture one object? Check the BCC Voronoi cell.

BCC on the integer scaling: points with all-even or all-odd coordinates.
  8 nearest neighbours  (+-1,+-1,+-1)  dist sqrt(3)   <- the diagonal 'hops'
  6 next-nearest        (+-2, 0, 0)    dist 2         <- the axial 'steps'
Claim: BOTH contribute facets to the Voronoi cell, giving 14 faces
       (8 hexagons + 6 squares) = the truncated octahedron = 3D permutohedron.
"""
import numpy as np
from itertools import product
from scipy.spatial import ConvexHull
from scipy.optimize import linprog

BCC = [p for p in product(range(-4, 5), repeat=3)
       if p != (0, 0, 0) and (all(c % 2 == 0 for c in p) or all(c % 2 for c in p))]
BCC.sort(key=lambda p: sum(c * c for c in p))

# halfspaces  x . v <= |v|^2 / 2  for every lattice neighbour v
V = np.array(BCC, float)
b = (V ** 2).sum(1) / 2


def is_facet(i):
    """Non-redundant iff we can strictly beat constraint i while satisfying the rest."""
    A = np.delete(V, i, 0)
    rhs = np.delete(b, i)
    # maximise x.v_i  subject to the other halfspaces  ->  does it exceed b_i ?
    r = linprog(-V[i], A_ub=A, b_ub=rhs, bounds=[(-10, 10)] * 3)
    return r.success and (-r.fun) > b[i] + 1e-9


facets = [BCC[i] for i in range(len(BCC)) if is_facet(i)]
shell = {}
for f in facets:
    shell.setdefault(round(sum(c * c for c in f), 6), []).append(f)

print("Voronoi facets of the BCC cell, by neighbour shell:")
for d2 in sorted(shell):
    ex = shell[d2][0]
    print(f"  |v|^2 = {d2:>2}  ({len(shell[d2]):>2} facets)  e.g. {ex}")
print(f"  TOTAL FACES = {len(facets)}")

# vertex count / face polygon sizes
Hn = np.array(facets, float)
Hb = (Hn ** 2).sum(1) / 2
# enumerate vertices as intersections of triples of facet planes inside the cell
verts = []
n = len(facets)
for i in range(n):
    for j in range(i + 1, n):
        for k in range(j + 1, n):
            M = Hn[[i, j, k]]
            if abs(np.linalg.det(M)) < 1e-9:
                continue
            x = np.linalg.solve(M, Hb[[i, j, k]])
            if np.all(Hn @ x <= Hb + 1e-7):
                verts.append(x)
verts = np.unique(np.round(np.array(verts), 6), axis=0)
print(f"  vertices = {len(verts)}")

sizes = {}
for i in range(n):
    on = np.sum(np.abs(verts @ Hn[i] - Hb[i]) < 1e-6)
    sizes[on] = sizes.get(on, 0) + 1
print("  face polygon sizes (verts per face -> how many faces):", dict(sorted(sizes.items())))

# permutohedron: hull of permutations of (1,2,3,4)
from itertools import permutations
P = np.array(list(permutations([1, 2, 3, 4])), float)
Pc = P - P.mean(0)
# drop to 3D via SVD (it lies in a hyperplane)
U, S, Wt = np.linalg.svd(Pc, full_matrices=False)
P3 = Pc @ Wt[:3].T
h = ConvexHull(P3)
print(f"\n3D permutohedron (hull of perms of 1..4): vertices = {len(np.unique(h.vertices))}, "
      f"facets(merged planes) = {len(np.unique(np.round(h.equations, 6), axis=0))}")
print("  proper non-empty subsets of {1,2,3,4} = 2^4-2 =", 2 ** 4 - 2,
      "  split C(4,2)=6 squares, C(4,1)+C(4,3)=8 hexagons")
