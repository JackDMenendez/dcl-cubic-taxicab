from fractions import Fraction as F
# x-marginal of the reflecting 3-D walk is itself a Markov chain on {-R..R}:
# +-x step w.p. 1/6 each (reflected at a wall -> inward); the 4 transverse steps
# leave x unchanged, so w.p. 4/6 x stays.  Aperiodic (unlike the full bipartite chain),
# which is why the -1 eigenvalue is invisible to the observable x.
def tau_exact(R):
    xs = list(range(-R, R+1)); n = len(xs); idx = {x:i for i,x in enumerate(xs)}
    P = [[F(0)]*n for _ in range(n)]
    for x in xs:
        i = idx[x]
        for d in (1,-1):
            y = x+d
            if abs(y) > R: y = x-d          # reflect: negate the step
            P[i][idx[y]] += F(1,6)
        P[i][i] += F(4,6)
    # stationary measure: w(+-R)=1/2, w=1 elsewhere, normalised
    w = [F(1,2) if abs(x)==R else F(1) for x in xs]; Z=sum(w); pi=[a/Z for a in w]
    # check stationarity
    for j in range(n):
        assert sum(pi[i]*P[i][j] for i in range(n)) == pi[j], "not stationary"
    var = sum(pi[i]*F(xs[i])**2 for i in range(n))
    # fundamental matrix: solve (I - P + Pi) y = x   (Pi_ij = pi_j)
    A=[[ (F(1) if i==j else F(0)) - P[i][j] + pi[j] for j in range(n)] for i in range(n)]
    b=[F(xs[i]) for i in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0)
        A[c],A[p]=A[p],A[c]; b[c],b[p]=b[p],b[c]
        inv=F(1)/A[c][c]; A[c]=[v*inv for v in A[c]]; b[c]*=inv
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[u-f*v for u,v in zip(A[r],A[c])]; b[r]-=f*b[c]
    tau = sum(pi[i]*F(xs[i])*b[i] for i in range(n))/var
    return var, tau, 2*var*tau

print(f"{'cell':<6}{'Var(x)':>10}{'tau (exact)':>16}{'S(0)=2*Var*tau':>18}   document")
doc = {1:("1/2","3","3"), 2:("3/2","10","30"), 3:("19/6","417/19","139")}
for R in (1,2,3):
    v,t,s = tau_exact(R)
    d = doc[R]
    print(f"{2*R+1}^3{'':<2}{str(v):>10}{str(t):>16}{str(s):>18}   "
          f"Var {d[0]}, tau {d[1]}, S(0) {d[2]}")
print()
print("tau / L^2 :", "  ".join(f"{float(tau_exact(R)[1])/(2*R+1)**2:.4f}" for R in (1,2,3)))
print("RMS / L   :", "  ".join(f"{float(tau_exact(R)[0])**0.5/(2*R+1):.4f}" for R in (1,2,3)))
print("f_c       :", "  ".join(f"{1/(2*3.141592653589793*float(tau_exact(R)[1])):.6f}" for R in (1,2,3)))
