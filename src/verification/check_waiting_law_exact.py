from fractions import Fraction as F
from collections import defaultdict
AX=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC=lambda p: len({c&1 for c in p})==1
P6=F(1,6); live={(0,0,0):F(1)}; Et=F(0)
print(" t | P(T=t) exact      | (1/3)(2/3)^(t-2)  | match")
for t in range(1,15):
    nxt=defaultdict(F); inc=F(0)
    for p,w in live.items():
        for s in AX:
            q=(p[0]+s[0],p[1]+s[1],p[2]+s[2])
            if onBCC(q): inc+=w*P6; Et+=w*P6*t
            else: nxt[q]+=w*P6
    live=nxt
    pred=F(1,3)*F(2,3)**(t-2) if t>=2 else F(0)
    print(f"{t:>2} | {str(inc):<17} | {str(pred):<17} | {'YES' if inc==pred else 'NO'}")
print("\nE[T] partial through t=14:", float(Et))
print("Kac's lemma: E[T] = [Z^3 : Lambda_BCC] = 4 exactly (independent of this sum)")
