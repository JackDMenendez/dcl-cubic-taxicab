"""Momentum = bias.  What survives, and what breaks.

An unbiased walk has zero mean displacement: a particle AT REST.  Give the step
distribution a drift and you have a particle with momentum.  The natural way to add
it is an exponential tilt (maximum entropy at fixed mean displacement):

    p(s) = exp(beta . s) / Z        Z = 2 cosh(beta) + 4   for a tilt along x

Questions:
  1. does the hop clock E[T] = 4 survive a bias?
  2. what is the drift velocity, and is it capped by the light cone?
  3. does the trembling rate slow down with speed, i.e. is there time dilation?
     Relativity wants internal cycles per coordinate tick to fall as 1/gamma.
"""
import numpy as np
from collections import defaultdict

AX = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC = lambda p: len({c & 1 for c in p}) == 1


def run(beta, H=4000, prune=1e-16):
    w = np.array([np.exp(beta), np.exp(-beta), 1, 1, 1, 1], float)
    w /= w.sum()
    mean_step = sum(wi*np.array(s) for wi, s in zip(w, AX))

    live = {(0,0,0): 1.0}
    ET = 0.0; disp = np.zeros(3); pnull = 0.0; mass = 0.0
    for t in range(1, H):
        nxt = defaultdict(float)
        for p, m in live.items():
            for wi, s in zip(w, AX):
                q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
                if onBCC(q):
                    mass += m*wi; ET += m*wi*t; disp += m*wi*np.array(q)
                    if q == (0,0,0): pnull += m*wi
                else:
                    nxt[q] += m*wi
        live = {k: v for k, v in nxt.items() if v > prune}
        if sum(live.values()) < 1e-14: break
    return dict(mass=mass, ET=ET/mass, disp=disp/mass, pnull=pnull/mass,
                mean_step=mean_step)


print("unconfined walk, exponential tilt along x")
print(f"  {'beta':>6}{'mass':>10}{'E[T]':>10}{'E[disp_x]':>12}{'v=disp/T':>11}"
      f"{'P(null)':>10}{'P0/P(0)':>10}{'1/gamma':>10}")
base = None
for beta in (0.0, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0):
    r = run(beta)
    v = r['disp'][0]/r['ET']
    if base is None: base = r['pnull']
    g = np.sqrt(max(0.0, 1 - v*v))
    print(f"  {beta:>6.2f}{r['mass']:>10.6f}{r['ET']:>10.6f}{r['disp'][0]:>12.6f}"
          f"{v:>11.6f}{r['pnull']:>10.6f}{r['pnull']/base:>10.6f}{g:>10.6f}", flush=True)

print("""
Reading:
  E[T] stays 4 for every finite bias.  Reason: on a finite group the uniform
  distribution is stationary for ANY step distribution (the transition matrix is
  doubly stochastic by translation invariance), so Kac's lemma needs only
  irreducibility -- and every strictly positive bias still reaches all three
  non-identity cosets of Z^3/BCC.  Only the singular limit beta -> infinity, where
  the walk moves in one direction alone, drops to E[T] = 2.

  So the hop clock is INDEPENDENT OF SPEED.  Internal cycles per coordinate tick
  stay at 1/4 no matter how fast the particle drifts.  Relativity wants 1/(4 gamma).
  Compare the last three columns: P(null) does fall with speed, but check whether it
  falls like 1/gamma before reading it as time dilation.
""")
