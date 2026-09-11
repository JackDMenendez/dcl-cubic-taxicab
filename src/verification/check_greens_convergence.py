import numpy as np
TARGET = 3/(2*np.pi)
def axis_S(L):
    """S(k1) = sum over k2,k3 of 1/(1-phi), computed slice-wise (low memory)."""
    n = np.fft.fftfreq(L, 1.0/L); k = 2*np.pi*n/L; c = np.cos(k)
    S = np.empty(L); tot = 0.0
    for i in range(L):
        d = 1.0 - (c[i] + c[:,None] + c[None,:])/3.0
        if i == 0: d[0,0] = np.inf
        inv = 1.0/d
        S[i] = inv.sum(); tot += S[i]
    return k, S, tot/L**3

print(f"target 3/(2pi) = {TARGET:.8f}\n")
print(f"{'L':>5} {'G(0)':>12} {'relerr':>9} | " + " ".join(f"r={r:<3}" for r in (8,16,24,32,48,64)))
for L in (192, 384, 768):
    k, S, G0 = axis_S(L)
    vol = L**3
    G = lambda r: float(np.cos(k*r) @ S)/vol
    cache = {r: G(r) for r in range(1, 70)}
    row = []
    for r in (8,16,24,32,48,64):
        E = abs(cache[r-1]-cache[r+1])/2
        row.append(f"{r*r*E:.5f}")
    print(f"{L:>5} {G0:12.8f} {abs(G0-1.5163860591)/1.5163860591:9.2e} | " + " ".join(f"{v:>6}" for v in row))

print("\nbest-plateau value and its distance from 3/(2pi), per L:")
for L in (192, 384, 768):
    k, S, G0 = axis_S(L)
    vol = L**3
    cache = {r: float(np.cos(k*r) @ S)/vol for r in range(1, 90)}
    best = min(((abs(r*r*abs(cache[r-1]-cache[r+1])/2 - TARGET), r,
                 r*r*abs(cache[r-1]-cache[r+1])/2) for r in range(4, 88)))
    print(f"  L={L:<4} closest at r={best[1]:<3} value {best[2]:.6f}  |diff| {best[0]:.2e}")
