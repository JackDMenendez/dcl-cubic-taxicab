import numpy as np, itertools
from itertools import product, combinations

stencil = [s for s in product((-1,0,1),repeat=3) if any(s)]
q = lambda s: int(sum(np.sign(s)))

# CLAIM A: on this stencil q(V) is exactly the linear functional (1,1,1).V
bad = [s for s in stencil if q(s) != s[0]+s[1]+s[2]]
print(f"q(V) == (1,1,1).V for all 26 stencil vectors?  {'YES' if not bad else 'NO '+str(bad)}")

# CLAIM B (document): "every inversion-breaking subset carries a net drift"
shell8 = [s for s in stencil if sum(map(abs,s)) == 3]
n = len(shell8)
counterex = []
stats = {"inv-breaking":0, "inv-breaking & zero drift":0,
         "nonzero charge":0, "nonzero charge & zero drift":0}
for mask in range(1, 1<<n):
    S = [shell8[i] for i in range(n) if mask>>i & 1]
    Sset = set(S); negset = {tuple(-c for c in v) for v in S}
    inv_breaking = Sset != negset
    drift = np.sum(np.array(S), axis=0)
    zero_drift = not drift.any()
    charge = sum(q(v) for v in S)
    if inv_breaking:
        stats["inv-breaking"] += 1
        if zero_drift:
            stats["inv-breaking & zero drift"] += 1
            if len(counterex) < 3: counterex.append((S, int(charge)))
    if charge: 
        stats["nonzero charge"] += 1
        if zero_drift: stats["nonzero charge & zero drift"] += 1

print(f"\nover all {(1<<n)-1} non-empty subsets of the 8-body shell:")
for k,v in stats.items(): print(f"  {k:<32} {v}")

print("\nCOUNTEREXAMPLES to 'every inversion-breaking subset carries a net drift':")
for S,c in counterex:
    print(f"  {S}   drift = (0,0,0)   net charge = {c}")

print("\nCLAIM C: net charge == (1,1,1) . (vector sum), for every subset?")
ok = all(sum(q(v) for v in [shell8[i] for i in range(n) if m>>i&1])
         == int(np.sum(np.array([shell8[i] for i in range(n) if m>>i&1]))) 
         for m in range(1, 1<<n))
print(f"  {'YES' if ok else 'NO'}  -> charge != 0  <=>  the drift has a nonzero (1,1,1) component")
