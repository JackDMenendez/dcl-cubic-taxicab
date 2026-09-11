"""Does the virtual excursion measure already contain the 1/r field?

G(x) = expected number of visits to x by a walk started at 0
     = (2pi)^-3 * integral cos(k.x) / (1 - phi(k)) d^3k ,   phi = (cos k1+cos k2+cos k3)/3

Classical result: for simple random walk on Z^3,  G(x) ~ 3 / (2 pi |x|).
So the visit density of the walk IS a Coulomb potential -- no 1/r^2 postulate needed;
the inverse-square law is the Green's function of a 3D walk.

Validation: G(0) should reproduce 1/(1 - p_return) = 1.5163860591...
"""
import numpy as np

L = 192
n = np.fft.fftfreq(L, 1.0 / L)          # integer modes
k = 2 * np.pi * n / L
c = np.cos(k)

# 1 / (1 - phi) on the full 3D mode grid, k = 0 excluded
denom = 1.0 - (c[:, None, None] + c[None, :, None] + c[None, None, :]) / 3.0
denom[0, 0, 0] = np.inf
Ginv = 1.0 / denom
vol = L ** 3

G0 = Ginv.sum() / vol
print(f"L = {L}")
print(f"G(0) numeric = {G0:.8f}      exact 1/(1-p_ret) = 1.51638606")
print(f"   relative error {abs(G0-1.5163860591)/1.5163860591:.2e}   (expected O(1/L))\n")

# --- axis direction, separable: sum out k2,k3 first
S = Ginv.sum(axis=(1, 2))                                  # S(k1)
print("axis direction  x = (r,0,0):")
print(f"  {'r':>3} {'G(r)':>12} {'r*G(r)':>10}   target 3/(2pi) = {3/(2*np.pi):.6f}")
for r in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48):
    G = float(np.cos(k * r) @ S) / vol
    print(f"  {r:>3} {G:12.8f} {r*G:10.6f}")

# --- body-diagonal direction, x = (r,r,r): |x| = r*sqrt(3)
print("\nbody-diagonal direction  x = (r,r,r):")
print(f"  {'r':>3} {'|x|':>8} {'G(x)':>12} {'|x|*G':>10}")
for r in (1, 2, 3, 4, 6, 8, 12, 16, 24):
    e = np.cos(k * r)
    s = np.sin(k * r)
    # Re[ exp(i r(k1+k2+k3)) ] expanded over the separable product
    E = (e[:, None, None] * e[None, :, None] * e[None, None, :]
         - e[:, None, None] * s[None, :, None] * s[None, None, :]
         - s[:, None, None] * e[None, :, None] * s[None, None, :]
         - s[:, None, None] * s[None, :, None] * e[None, None, :])
    G = float((E * Ginv).sum()) / vol
    d = r * np.sqrt(3)
    print(f"  {r:>3} {d:8.4f} {G:12.8f} {d*G:10.6f}")

# --- the field: discrete gradient magnitude along the axis
print("\nfield strength  |G(r-1)-G(r+1)|/2  along the axis  (expect ~ 1/r^2):")
Gax = {r: float(np.cos(k * r) @ S) / vol for r in range(1, 52)}
print(f"  {'r':>3} {'E(r)':>12} {'r^2*E(r)':>10}   target 3/(2pi) = {3/(2*np.pi):.6f}")
for r in (2, 3, 4, 6, 8, 12, 16, 24, 32, 40):
    E = abs(Gax[r - 1] - Gax[r + 1]) / 2
    print(f"  {r:>3} {E:12.8f} {r*r*E:10.6f}")
