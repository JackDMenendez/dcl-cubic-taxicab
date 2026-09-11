import numpy as np
from collections import defaultdict
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'..','utilities','biased.py')).read().split('print("unconfined')[0])
base = run(0.0)['pnull']
print(f"{'beta':>6}{'v':>11}{'P(null)':>11}{'P0/P0(rest)':>13}{'1/gamma':>10}"
      f"{'(1-ratio)/v^2':>15}{'(1-1/g)/v^2':>13}")
for b in (0.05, 6.0):
    r = run(b); v = r['disp'][0]/r['ET']; g = np.sqrt(max(0.,1-v*v)); rat = r['pnull']/base
    print(f"{b:>6.2f}{v:>11.6f}{r['pnull']:>11.6f}{rat:>13.6f}{g:>10.6f}"
          f"{(1-rat)/v**2:>15.4f}{(1-g)/v**2:>13.4f}   E[T]={r['ET']:.6f}")
print("document: v=0.016667 -> P0/P0=0.998885, 1/g=0.999861, 4.0153, 0.5000")
print("document: beta=6.0   -> v=0.990170, P(null)=0.000036, E[T]=4.000000")
