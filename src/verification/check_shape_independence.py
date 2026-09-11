from itertools import product
from fractions import Fraction as F
import numpy as np, sys
sys.path.insert(0,'.')
AX=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC=lambda p: len({c&1 for c in p})==1

def solve(inside, pts, mode):
    def targets(p):
        legal,illegal=[],[]
        for s in AX:
            q=(p[0]+s[0],p[1]+s[1],p[2]+s[2])
            (legal if inside(q) else illegal).append((s,q))
        if mode=="resample":
            w=F(1,len(legal)); return [(q,w) for _,q in legal]
        out=[(q,F(1,6)) for _,q in legal]
        if mode=="stay":
            if illegal: out.append((p,F(len(illegal),6)))
        else:
            for s,_ in illegal:
                r=(p[0]-s[0],p[1]-s[1],p[2]-s[2])
                assert inside(r), f"reflect illegal from {p} step {s}"
                out.append((r,F(1,6)))
        return out
    trans=[p for p in pts if not onBCC(p)]
    idx={p:i for i,p in enumerate(trans)}; n=len(trans)
    A=[[F(0)]*n for _ in range(n)]; b=[F(1)]*n
    for p in trans:
        i=idx[p]; A[i][i]+=F(1)
        for q,w in targets(p):
            if not onBCC(q): A[i][idx[q]]-=w
    for c in range(n):
        piv=next(x for x in range(c,n) if A[x][c]!=0)
        A[c],A[piv]=A[piv],A[c]; b[c],b[piv]=b[piv],b[c]
        inv=F(1)/A[c][c]; A[c]=[v*inv for v in A[c]]; b[c]*=inv
        for x in range(n):
            if x!=c and A[x][c]!=0:
                f=A[x][c]; A[x]=[u-f*v for u,v in zip(A[x],A[c])]; b[x]-=f*b[c]
    ET=F(0)
    for q,w in targets((0,0,0)):
        ET += w*(F(0) if onBCC(q) else b[idx[q]])
    return ET+1

SHAPES=[]
for R in (1,2,3):
    ins=(lambda R: (lambda p: max(map(abs,p))<=R))(R)
    SHAPES.append((f"cube R={R}", ins, [p for p in product(range(-R,R+1),repeat=3) if ins(p)]))
ins=lambda p: abs(p[0])<=1 and abs(p[1])<=3 and abs(p[2])<=3
SHAPES.append(("slab |x|<=1,|y|,|z|<=3", ins,
               [p for p in product(range(-1,2),range(-3,4),range(-3,4)) if ins(p)]))
ins=lambda p: sum(map(abs,p))<=3
SHAPES.append(("L1 ball r<=3", ins,
               [p for p in product(range(-3,4),repeat=3) if ins(p)]))

print(f"{'box':<26}{'reflect':>12}{'stay':>12}{'resample':>12}")
for nm, ins, pts in SHAPES:
    row=[]
    for mode in ("reflect","stay","resample"):
        try: row.append(f"{float(solve(ins,pts,mode)):.6f}")
        except AssertionError as e: row.append("REFL-ILLEGAL")
    print(f"{nm:<26}" + "".join(f"{v:>12}" for v in row))
print("\ndocument sec 5.8:")
print("  cube R=1 4.000000 5.000000 4.000000 / R=2 4.000000 4.120000 3.960000")
print("  cube R=3 4.000000 4.029224 3.993112 / slab 4.000000 4.433656 4.078532")
print("  L1 ball  4.000000 4.545455 3.818182")
