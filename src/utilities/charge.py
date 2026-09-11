"""Charge on the taxispace stencil: why a step-linear charge functional cannot work,
and what the coupling constant turns out to be.

Part 1  q(V) = sum of component signs ('plus to plus, minus to minus, none to zero')
        is ODD under inversion V -> -V.  Every shell is closed under inversion, so
        every shell has exactly zero net charge.  No accident, no tuning available.

Part 2  Restoring nonzero charge requires an inversion-breaking subset -- and every
        such subset carries a net drift.  A 'charged' particle would self-propel and
        its field would be anisotropic.  This is a no-go, not a difficulty.

Part 3  The Coulomb coefficient of the walk's Green's function is fixed by the mean
        square step length:   G(x) -> 3 / (2 pi <|s|^2> |x|).
        The measured 0.477465 for the axial walk is exactly 6/(4pi).

Part 4  Charge as 'which BCC sublattice' fails: 8-type hops flip parity and they are
        33.6% of hops, so charge would flip on a third of them.
"""
from itertools import product
from collections import Counter
import numpy as np

shell = {1: [], 2: [], 3: []}
for s in product((-1, 0, 1), repeat=3):
    w = sum(map(abs, s))
    if w: shell[w].append(s)
NAMES = {1: "6 axial", 2: "12 face", 3: "8 body"}

q = lambda s: int(sum(np.sign(s)))          # the proposed charge functional

print("Part 1  q(V) = sum of component signs, per shell")
for w in (1, 2, 3):
    S = shell[w]; qs = [q(s) for s in S]
    vs = tuple(int(v) for v in np.sum(np.array(S), axis=0))
    print(f"  {NAMES[w]:<8} n={len(S):>2}  q values {dict(sorted(Counter(qs).items()))}"
          f"   NET = {sum(qs)}   shell vector sum = {vs}")
print("  q is odd under V -> -V and every shell is inversion-closed => NET = 0 always.")

print("\nPart 2  the inversion-breaking subset needed for nonzero charge")
for w in (2, 3):
    P = [s for s in shell[w] if q(s) > 0]
    v = np.sum(np.array(P), axis=0)
    print(f"  {NAMES[w]:<8} q>0 subset {P}")
    note = "= 2/sqrt(3)" if w == 2 else "= sqrt(3)/2"
    print(f"           vector sum = {tuple(int(x) for x in v)}   "
          f"drift per step = {np.linalg.norm(v)/len(P):.6f}   ({note})")
print("  => net drift along the body diagonal: the particle self-propels, field anisotropic.")

print("\nPart 3  Coulomb coefficient from the mean square step length")
print(f"  {'step set':<9} {'<|s|^2>':>8} {'coefficient':>13}   as flux/(4pi)")
for w in (1, 2, 3):
    S = np.array(shell[w]); m2 = (S ** 2).sum(1).mean()
    coef = 3 / (2 * np.pi * m2)
    print(f"  {NAMES[w]:<9} {m2:8.1f} {coef:13.6f}   {2*3/m2:.0f}/(4pi)")
print(f"  6/(4pi) = {6/(4*np.pi):.6f}   <- matches the measured r^2|grad G| = 0.4775")
print("  derivation: (I-P)G = delta  =>  -(<|s|^2>/6) Lap G = delta"
      "  =>  Lap G = -(6/<|s|^2>) delta")

print("\nPart 4  does a hop preserve the BCC sublattice (all-even vs all-odd)?")
for nm, v, p in (("6-type  (2,0,0)", (2,0,0), 0.2828),
                 ("8-type  (1,1,1)", (1,1,1), 0.3358),
                 ("null    (0,0,0)", (0,0,0), 0.2187),
                 ("        (2,2,0)", (2,2,0), 0.0735)):
    keep = all(c % 2 == 0 for c in v)
    print(f"  {nm:<16} {'preserves' if keep else 'FLIPS    '}   hop share {p:.4f}")
print("  => sublattice-as-charge flips on 33.6% of hops. Not viable as charge.")
print("     (and 0.336 != 0.5, so it is not the f_beat = 0.5 - f_zitt alternation either)")
