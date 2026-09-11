from itertools import product
from fractions import Fraction as F
AX=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC=lambda p: len({c&1 for c in p})==1
ins=lambda p: sum(map(abs,p))<=3
pts=[p for p in product(range(-3,4),repeat=3) if ins(p)]

def targets(p, variant):
    legal,illegal=[],[]
    for s in AX:
        q=(p[0]+s[0],p[1]+s[1],p[2]+s[2])
        (legal if ins(q) else illegal).append((s,q))
    out=[(q,F(1,6)) for _,q in legal]; stuck=[]
    for s,_ in illegal:
        r=(p[0]-s[0],p[1]-s[1],p[2]-s[2])
        if ins(r): out.append((r,F(1,6)))
        else:      stuck.append(s)
    if stuck:
        if variant=="refl-else-stay":
            out.append((p,F(len(stuck),6)))
        elif variant=="refl-else-resample":
            tot=F(len(stuck),6); n=len(out)
            out=[(q,w+tot/n) for q,w in out]
    return out

def ET(variant):
    trans=[p for p in pts if not onBCC(p)]
    idx={p:i for i,p in enumerate(trans)}; n=len(trans)
    A=[[F(0)]*n for _ in range(n)]; b=[F(1)]*n
    for p in trans:
        i=idx[p]; A[i][i]+=F(1)
        for q,w in targets(p,variant):
            if not onBCC(q): A[i][idx[q]]-=w
    for c in range(n):
        piv=next(x for x in range(c,n) if A[x][c]!=0)
        A[c],A[piv]=A[piv],A[c]; b[c],b[piv]=b[piv],b[c]
        inv=F(1)/A[c][c]; A[c]=[v*inv for v in A[c]]; b[c]*=inv
        for x in range(n):
            if x!=c and A[x][c]!=0:
                f=A[x][c]; A[x]=[u-f*v for u,v in zip(A[x],A[c])]; b[x]-=f*b[c]
    tot=F(0)
    for q,w in targets((0,0,0),variant):
        tot += w*(F(0) if onBCC(q) else b[idx[q]])
    return tot+1

print("L1 ball r<=3, reflect variants (the pure s -> -s rule is UNDEFINED here):")
for v in ("refl-else-stay","refl-else-resample"):
    e=ET(v); print(f"  {v:<22} E[T] = {str(e):>14} = {float(e):.6f}")
print(f"  document sec 5.8 claims: reflect  4.000000")
print("\nWhen IS s -> -s always legal?  Test regions:")
def always_legal(inside, pts):
    for p in pts:
        if onBCC(p): continue
        for s in AX:
            q=(p[0]+s[0],p[1]+s[1],p[2]+s[2])
            if not inside(q):
                r=(p[0]-s[0],p[1]-s[1],p[2]-s[2])
                if not inside(r): return False
    return True
tests=[("cube R=2", lambda p: max(map(abs,p))<=2,
        [p for p in product(range(-2,3),repeat=3)]),
       ("slab |x|<=1,|y|,|z|<=3", lambda p: abs(p[0])<=1 and abs(p[1])<=3 and abs(p[2])<=3,
        [p for p in product(range(-1,2),range(-3,4),range(-3,4))]),
       ("box 3x5x7", lambda p: abs(p[0])<=1 and abs(p[1])<=2 and abs(p[2])<=3,
        [p for p in product(range(-1,2),range(-2,3),range(-3,4))]),
       ("L1 ball r<=3", ins, pts),
       ("L2 ball r<=3", lambda p: sum(c*c for c in p)<=9,
        [p for p in product(range(-3,4),repeat=3) if sum(c*c for c in p)<=9])]
for nm, f, P in tests:
    print(f"  {nm:<24} s -> -s always legal? {'YES' if always_legal(f,[p for p in P if f(p)]) else 'NO'}")
