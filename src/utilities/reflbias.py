"""Bias carried by the WALLS instead of by the bulk step distribution.

Two ways to give the walker momentum:

  TILT      -- reweight the six bulk steps, p(s) ~ exp(beta.s).
               Momentum lives in the bulk.  Side effect: the step covariance becomes
               anisotropic, and (see below) with the WRONG SIGN for a moving charge.

  WALL BIAS -- leave the bulk perfectly uniform and isotropic; make the reflection
               probability direction-dependent at the cell boundary.  Momentum lives
               entirely in the boundary condition.

Two wall models:
  R  reflect (s -> -s) with prob r_w, else STAY.  Hop only on BCC arrival.
  T  reflect (s -> -s) with prob r_w, else TRANSMIT: the walker leaves and that exit
     IS the hop.

CORRECTION (verified 2026-09-10).  An earlier version of this docstring claimed that
because -s == s mod 2, "the induced quotient walk is untouched and E[T] = 4 exactly"
for both wall models.  THE OUTPUT BELOW REFUTES THAT: model T gives E[T] = 3.959316
and model R gives 4.022767.  The section 5.8 theorem covers PURE reflection only --
transmission adds an absorption channel and 'stay' adds a quotient self-loop, and
either one perturbs the hop clock in a box-dependent way.  Note the r = 1.00/1.00 row
does return exactly 4.000000: that IS pure reflection, where the theorem applies.
So the theorem's scope is right; only the old docstring overreached.
"""
import numpy as np
from collections import defaultdict

AX = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
onBCC = lambda p: len({c & 1 for c in p}) == 1
Linf  = lambda p: max(map(abs,p))


# ---------- the tilt, for comparison: covariance anisotropy and its sign ----------
print("TILT: does the bulk step covariance stay isotropic?  (it does not)")
print(f"  {'beta':>6}{'<s_x^2>':>10}{'<s_y^2>':>10}{'ratio':>9}   longitudinal spread")
for beta in (0.0, 0.5, 1.0, 2.0):
    w = np.array([np.exp(beta), np.exp(-beta), 1,1,1,1], float); w /= w.sum()
    sx = sum(wi*s[0]**2 for wi, s in zip(w, AX))
    sy = sum(wi*s[1]**2 for wi, s in zip(w, AX))
    tag = "same" if abs(sx-sy) < 1e-12 else ("STRETCHED along motion" if sx > sy else "compressed")
    print(f"  {beta:>6.2f}{sx:>10.6f}{sy:>10.6f}{sx/sy:>9.4f}   {tag}")
print("  A moving charge's field is COMPRESSED longitudinally (E_par ~ 1/gamma^2,")
print("  E_perp ~ gamma).  The tilt stretches it instead -- wrong sign.\n")


# ---------- wall-biased walks ----------
def run(R, rplus, rminus, model, H=3000, prune=1e-15):
    """r_w for the +x and -x walls; all transverse walls perfectly reflecting."""
    def moves(p):
        out = []
        for s in AX:
            q = (p[0]+s[0], p[1]+s[1], p[2]+s[2])
            if Linf(q) <= R:
                out.append((q, 1/6, False)); continue
            r = rplus if s == (1,0,0) else rminus if s == (-1,0,0) else 1.0
            back = (p[0]-s[0], p[1]-s[1], p[2]-s[2])
            out.append((back, r/6, False))
            if r < 1.0:
                if model == "T": out.append((q, (1-r)/6, True))    # transmit = hop
                else:            out.append((p, (1-r)/6, False))   # stay
        return out

    live = {(0,0,0): 1.0}
    ET = 0.0; disp = np.zeros(3); mass = 0.0
    for t in range(1, H):
        nxt = defaultdict(float)
        for p, m in live.items():
            for q, w, transmit in moves(p):
                if onBCC(q) or transmit:
                    mass += m*w; ET += m*w*t; disp += m*w*np.array(q)
                else:
                    nxt[q] += m*w
        live = {k: v for k, v in nxt.items() if v > prune}
        if sum(live.values()) < 1e-13: break
    return ET/mass, disp/mass, mass


for model, label in (("T", "reflect-or-TRANSMIT (exit = hop)"),
                     ("R", "reflect-or-STAY (hop on BCC only)")):
    print(f"5x5x5 cell, model {model}: {label}")
    print(f"  {'r(+x)':>7}{'r(-x)':>7}{'E[T]':>11}{'drift_x':>11}{'v=d/T':>10}{'mass':>10}")
    for rp, rm in ((1.0,1.0), (0.9,0.9), (0.8,1.0), (0.6,1.0), (0.3,1.0), (0.0,1.0)):
        try:
            ET, d, mass = run(2, rp, rm, model)
            print(f"  {rp:>7.2f}{rm:>7.2f}{ET:>11.6f}{d[0]:>11.6f}{d[0]/ET:>10.6f}{mass:>10.6f}")
        except Exception as e:
            print(f"  {rp:>7.2f}{rm:>7.2f}   failed: {e}")
    print()
