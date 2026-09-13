"""Noise contributed by the random walk inside the confined cell.

Everything needed is already fixed: the reflecting transition matrix, its stationary
measure pi ~ 2^-(walls touched), and the hop statistics.  So the noise is fully
characterised -- no new assumptions.

Two components, on different timescales:

  INTRA-HOP  the walker jitters inside a BOUNDED cell.  A finite reversible chain has
             an autocorrelation that is a sum of exponentials, so the power spectrum is
             a sum of LORENTZIANS -- flat below the corner frequency, rolling off above.
             Bounded means NO low-frequency divergence and NO 1/f.

  INTER-HOP  the cell centre itself performs a random walk.  Position of a diffusing
             object has a 1/f^2 spectrum at low frequency.

Total: 1/f^2 at low frequency from cell-centre diffusion, with a Lorentzian trembling
shoulder at the intra-cell relaxation rate, flat between, Nyquist cutoff at the tick.
"""
import numpy as np
from itertools import product

AX = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
Linf = lambda p: max(map(abs, p))


def cell(R):
    sites = [p for p in product(range(-R, R+1), repeat=3)]
    idx = {p: i for i, p in enumerate(sites)}
    n = len(sites)
    P = np.zeros((n, n))
    for p in sites:
        i = idx[p]
        for s in AX:
            q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
            if Linf(q) > R:                      # reflect: take the negated step
                q = (p[0]-s[0], p[1]-s[1], p[2]-s[2])
            P[i, idx[q]] += 1/6
    return np.array(sites, float), P


for R in (1, 2, 3):
    S, P = cell(R)
    n = len(S)
    w, V = np.linalg.eig(P.T)
    k = np.argmin(abs(w - 1))
    pi = np.real(V[:, k]); pi /= pi.sum()

    walls = np.sum(np.abs(S) == R, axis=1)
    pred = 2.0 ** (-walls); pred /= pred.sum()

    varx = float((pi * S[:, 0]**2).sum())
    ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
    lam2 = ev[1]                                  # second-largest modulus
    tau_relax = -1/np.log(lam2)

    # position autocorrelation of x under stationarity
    x = S[:, 0]
    C = []
    Pt = np.eye(n)
    for t in range(40):
        C.append(float((pi * x * (Pt @ x)).sum()))
        Pt = Pt @ P
    C = np.array(C)
    tau_int = C.sum()/C[0]                        # integrated correlation time

    print(f"=== reflecting cell, L-inf <= {R}  ({n} sites) ===")
    print(f"  stationary measure matches 2^-walls : {np.allclose(pi, pred)}")
    print(f"  Var(x) = {varx:.6f}   RMS jitter = {np.sqrt(varx):.6f} lattice units")
    print(f"  |lambda_2| = {lam2:.6f}   relaxation time = {tau_relax:.4f} ticks")
    print(f"  integrated correlation time = {tau_int:.4f} ticks")
    print(f"  corner frequency f_c = 1/(2 pi tau) = {1/(2*np.pi*tau_int):.6f} per tick")
    print(f"  C(t)/C(0) for t=0..6: " + " ".join(f"{v:.4f}" for v in (C/C[0])[:7]))
    # white-noise level of the bounded jitter: S(0) = 2 * Var * tau_int
    print(f"  low-frequency plateau  S(0) = 2*Var*tau = {2*varx*tau_int:.6f}"
          f"   (units: lattice^2 per unit frequency)")
    print()

print("Interpretation")
print("  * the intra-cell noise is BOUNDED: S(f) is flat below f_c and rolls off as")
print("    1/f^2 above it -- a Lorentzian.  A finite cell cannot produce 1/f noise.")
print("  * the cell-centre random walk contributes 1/f^2 in POSITION at low frequency,")
print("    which dominates below f_c.  Its velocity spectrum is white.")
print("  * so the predicted shape is: 1/f^2 (centre diffusion) -> plateau -> 1/f^2")
print("    (trembling roll-off) -> Nyquist cutoff at the tick rate.  Three regimes,")
print("    two corners, no free parameters beyond the cell size.")
