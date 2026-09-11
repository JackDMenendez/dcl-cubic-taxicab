import numpy as np, itertools
from scipy.spatial import HalfspaceIntersection, ConvexHull

def poly(normals):
    V = np.array(normals, float); b = (V**2).sum(1)/2
    hs = np.hstack([V, -b[:,None]])
    try:
        hi = HalfspaceIntersection(hs, np.zeros(3))
    except Exception as e:
        return None, f"unbounded/failed: {type(e).__name__}"
    h = ConvexHull(hi.intersections)
    return h.volume, f"{len(np.unique(np.round(h.equations,6),axis=0))} facets, {len(h.vertices)} verts"

BD = [v for v in itertools.product((-1,1),repeat=3)]
AX = [v for v in itertools.permutations([2,0,0])] + [v for v in itertools.permutations([-2,0,0])]
AX = sorted(set(AX))

for name, nrm in [("8 body diagonals ONLY", BD),
                  ("6 axial ONLY", AX),
                  ("both (true BCC Voronoi cell)", BD+AX)]:
    vol, desc = poly(nrm)
    print(f"{name:<32} volume = {vol!s:<8.8}  {desc}")

print("\nBCC lattice index in Z^3 = 4  ->  fundamental-domain volume must be 4")
