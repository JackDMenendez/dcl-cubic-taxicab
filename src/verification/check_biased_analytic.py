import numpy as np
# tilt along x: p(+x) = e^b/Z, p(-x) = e^-b/Z, four transverse at 1/Z, Z = 2cosh b + 4
v = lambda b: 2*np.sinh(b)/(2*np.cosh(b)+4)
print("Wald:  E[disp_x] = E[T] * mean_step_x = 4 * v   (since E[T] = 4 for any finite bias)")
print(f"{'beta':>6}{'v = mean_step_x':>18}{'4v = E[disp_x]':>17}   document")
doc = {0.0:(0.000000,0.000000), 0.4:(0.533259,0.133315), 1.0:(1.326756,0.331689),
       2.0:(2.517693,0.629423), 3.0:(3.320569,0.830142), 6.0:(None,0.990170)}
for b in (0.0,0.4,1.0,2.0,3.0,6.0):
    vv=v(b); d=doc[b]
    dd = f"disp {d[0]}, v {d[1]}" if d[0] is not None else f"v {d[1]}"
    print(f"{b:>6.1f}{vv:>18.6f}{4*vv:>17.6f}   {dd}")
print(f"\nbeta -> infinity:  v -> e^b/(e^b+4) -> 1  (the axial light-cone speed)")
for b in (10,20,40): print(f"   beta={b:<3} v = {v(b):.10f}")

print("\nWhich beta produced the document's small-v dilation table?  (v = mean_step_x)")
for b in (0.05,0.1,0.2):
    print(f"   beta={b:<5} v = {v(b):.6f}   <- document rows 0.016667 / 0.033333 / 0.066666")

print("\nSmall-v expansion of the two rates:")
print("   1/gamma = sqrt(1-v^2)  ->  1 - v^2/2      so (1 - 1/gamma)/v^2 -> 0.5")
print("   document measures      ->  1 - 4 v^2      so (1 - ratio)/v^2   -> 4.0")
print("   ratio of the two coefficients = 4 / 0.5 = 8   (the flat factor of 8)")
