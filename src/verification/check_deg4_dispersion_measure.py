import numpy as np
from itertools import product
AX=[s for s in product((-1,0,1),repeat=3) if sum(map(abs,s))==1]
BD=[s for s in product((-1,0,1),repeat=3) if sum(map(abs,s))==3]

def disp_aniso(a):
    """a = total probability on the axial shell, 1-a on the body shell.
       Return relative anisotropy of E(k)/k^2 between axis and body diagonal, /k^2."""
    steps = [(s, a/len(AX)) for s in AX] + [(s, (1-a)/len(BD)) for s in BD]
    def E(kvec):
        kv = np.array(kvec, float)
        return 1.0 - sum(p*np.cos(kv @ np.array(s)) for s,p in steps)
    out=[]
    for kmag in (0.05, 0.01, 0.002):
        ax = E((kmag,0,0))/kmag**2
        d  = kmag/np.sqrt(3)
        bd = E((d,d,d))/kmag**2
        out.append(((ax-bd)/ax)/kmag**2)
    return out

# unit-vector moment nulling mixture vs |s|^4-weighted nulling mixture
print("mixture axial:body   -> relative degree-4 dispersion anisotropy / k^2")
print(f"{'a (axial share)':>16} {'k=0.05':>12} {'k=0.01':>12} {'k=0.002':>12}   note")
for a, note in ((1.0, "pure axial"),
                (2/5, "2/5 : 3/5  <- sec 9.5 unit-vector nulling"),
                (6/7, "6/7 : 1/7  <- |s|^4-weighted nulling"),
                (0.0, "pure body")):
    v = disp_aniso(a)
    print(f"{a:>16.6f} " + " ".join(f"{x:>12.6f}" for x in v) + f"   {note}")

print("\ncubic invariant <n1^4+n2^4+n3^4> under the two weightings:")
for a in (2/5, 6/7):
    N_ax=np.array(AX,float); N_bd=np.array(BD,float)
    n_ax=N_ax/np.linalg.norm(N_ax,axis=1,keepdims=True)
    n_bd=N_bd/np.linalg.norm(N_bd,axis=1,keepdims=True)
    for tag, s4ax, s4bd in (("unit-vector (sec 9.5)",1.0,1.0), ("|s|^4-weighted",1.0,9.0)):
        wa, wb = a*s4ax, (1-a)*s4bd
        cub = (wa*(n_ax**4).sum(1).mean() + wb*(n_bd**4).sum(1).mean())/(wa+wb)
        print(f"  a={a:.4f}  {tag:<24} <n^4> = {cub:.6f}  dev {cub-0.6:+.6f}")
