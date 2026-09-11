import numpy as np, itertools
from fractions import Fraction as F
from itertools import product

def expand(cls_probs, keep_null):
    out=[]
    for k,p in cls_probs.items():
        if k==(0,0,0):
            if keep_null: out.append((np.zeros(3), float(p)))
            continue
        S=set()
        for perm in set(itertools.permutations(k)):
            for sg in product(*[(1,-1) if c else (1,) for c in perm]):
                S.add(tuple(c*s for c,s in zip(perm,sg)))
        for v in S: out.append((np.array(v,float), float(p)/len(S)))
    return out

REFLECT5 = {(1,1,1):F(2,5), (2,0,0):F(61,210), (0,0,0):F(47,210),
            (2,2,0):F(8,105), (2,2,2):F(1,105)}
STAY5    = {(1,1,1):F(9,25), (2,0,0):F(3,10), (0,0,0):F(11,50),
            (2,2,0):F(1,10), (2,2,2):F(1,50)}
BODY3    = {(1,1,1):F(2,5), (0,0,0):F(3,5)}

def moments(steps):
    W=np.array([p for _,p in steps]); V=np.array([v for v,_ in steps])
    n=V/np.linalg.norm(V,axis=1,keepdims=True); s4=(V**2).sum(1)**2
    return (W*(n**4).sum(1)).sum()/W.sum(), (W*s4*(n**4).sum(1)).sum()/(W*s4).sum()

def hop_dispersion(steps):
    res=[]
    tot=sum(p for _,p in steps)
    for km in (0.02, 0.005):
        E=lambda kv: tot-sum(p*np.cos(np.array(kv)@v) for v,p in steps)
        ax=E((km,0,0))/km**2; d=km/np.sqrt(3); bd=E((d,d,d))/km**2
        res.append(((ax-bd)/ax)/km**2)
    return res

print(f"{'construction':<24}{'unit dev':>11}{'|v|^4 dev':>11}"
      f"{'unit rel':>10}{'|v|^4 rel':>11}  hop-dispersion aniso/k^2")
for nm, cp in (("5x5x5 reflect", REFLECT5), ("5x5x5 stay/resample", STAY5),
               ("3x3x3 (any conv.)", BODY3)):
    u,w = moments(expand(cp, False))
    disp = hop_dispersion(expand(cp, True))
    print(f"{nm:<24}{u-0.6:>+11.6f}{w-0.6:>+11.6f}"
          f"{(u-0.6)/0.6:>+9.2%}{(w-0.6)/0.6:>+10.2%}  "
          + "  ".join(f"{x:+.6f}" for x in disp))
print("\nsanity: E(0) must be 0 ->",
      f"{sum(p for _,p in expand(REFLECT5,True)) - sum(p*1 for _,p in expand(REFLECT5,True)):.1e}")
