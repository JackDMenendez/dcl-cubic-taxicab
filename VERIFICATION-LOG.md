# Appendix A script verification log — V1.4 cubic taxicab architecture

Board #32 (031). Run with `C:\Users\jackd\.venv-win\Scripts\python.exe` (Py 3.14.7,
numpy 2.5.0, scipy 1.18.0, sympy 1.14.0; networkx MISSING).
Scripts staged in `scratchpad/scripts/`. NOT yet committed to dcl-mathematics.

## Status: 11 of 11 recovered — PASS COMPLETE

| script | covers | verdict |
|---|---|---|
| `hop_fpd.py` | §5.1 rule A first-passage | ✅ exact reproduction |
| `both68.py` | §6 BCC Voronoi + permutohedron | ✅ numbers reproduce; ⚠️ one prose claim refuted |
| `branch68.py` | §5.3/§5.4 first-return | ✅ exact; waiting law verified past documented range |
| `shells_6_12_8.py` | §7.1/7.2/7.3 Kac + all three cells | ✅ exact; no findings |
| `greens.py` | §8.2 Green's fn, 3/(2π) | ✅ exact; ⭐ document UNDERclaims (F7) |
| `charge.py` | §8.5 no-go, §8.3 coupling | ✅ numbers exact; ⚠️ universal refuted (F8), ⭐ F9 |
| `lightcone.py` | §9.2 cone, §9.3 covariance | ✅ numbers exact; 🚨 PRE-CORRECTION docstring (F10) |
| `isotropy.py` | §9.5 moments + directions | ✅ numbers exact; 🚨 §9.5 route is WRONG (F11) |
| `confined.py` | §5.6 confined cells | ✅ POST-correction, exact incl. -0.14%; F12, F13 |
| `biased.py` | §9.6 tilt, time-dilation failure | ✅ exact, all rows; no findings |
| `reflbias.py` | §9.7 tilt covariance, wall bias | ✅ exact; 🚨 stale docstring (F14), F15 |

## Verified exactly

**`hop_fpd.py` → §5.1.** All rationals match: (3,0,0) = 1/21, (2,1,0) = 4/7,
(1,1,1) = 8/21, E[T] = 31/7, overshoot = 0 exactly. Weighted 3×-face corner share
= 3017/18471 = 0.163337, matching "0.163, below the uniform 0.211".
Multinomial rule confirmed by hand from output: p×126 = 1 / 3 / 6 for
axial / mixed / corner. Enhancement 1.81×, suppression 3.3× both check.
Runtime 0.54 s.

*Unrecorded extras:* the script computes 5 configurations (A–E); §5.1 reports only
A and D. B (18-nb) gives E[T] = 13683/4414 ≈ 3.10, overshoot 0.402; C (26-nb) gives
E[T] = 163189/71698 ≈ 2.28, overshoot 0.297. That is the §7.4 dial under rule A,
currently undocumented.

**`both68.py` → §6.** Facets 8 (|v|²=3) + 6 (|v|²=4) = 14, vertices 24, polygon
sizes {4: 6 faces, 6: 8 faces}, permutohedron 24/14. Example vectors match too.
Runtime 3.5 s.

**`branch68.py` → §5.3, §5.4, §5.5 P/site column.** All eight-decimal figures match:
body 0.33580514, axial 0.28282414, HOME 0.21871420, farther 0.16265652,
8/(6+8) = 0.54282128, per-site body/axial = 0.890496, face/site 0.00613112,
(4,0,0)/(3,1,1)/(2,2,2) per site 7.909e-4 / 2.259e-3 / 1.173e-3. Runtime 38 s
(exact Fractions, 60-step horizon, unabsorbed 4.08e-11).

Waiting law `P(T=t) = (1/3)(2/3)^(t-2)` confirmed as **exact rationals to t = 14**
(`_check_geom.py`) — document claimed verification only to t = 8. Two-step dead time
`P(T=1) = 0` exact. Memorylessness / renewal-process argument is solid.

**`shells_6_12_8.py` → §7.1, §7.2, §7.3.** Clean reproduction, no findings.
Part 1 (exact Fractions): 6→4, 12→4, 8→1, 18→4, 26→4. Part 2: all five branching
rows match to 4 dp incl. the structural zero (12-face-only reaches 8-type at exactly
0.0000) and the §7.4 dial (diagonal share 0 → 0.543 → 0.616 → 0.813; null 0.2187 → 0.0481).
Part 3: SC 6/8 all deg-3; FCC 12 rhombi/14 verts (8 deg-3 + 6 deg-4); BCC 14 {4:6, 6:8}/24
all deg-3. Runtime 3.8 s.

Two free cross-checks fell out: the BCC cell is re-derived through a different code path
than `both68.py` and agrees (14 faces, 24 verts); and the duality claim is confirmed
directly — the 8 deg-3 and 6 deg-4 vertices of the FCC cell are the same 8 and 6 as the
BCC cell's hexagons and squares.

**`greens.py` → §8.2.** Exact reproduction: G(0) = 1.50933032, relerr 4.65e-03, and all
four field rows (r=8/12/16/24 → 0.491382 / 0.482949 / 0.479601 / 0.474837). Wrap-around
caveat confirmed: r*G(r) sags 0.509 → 0.156 by r=48 while the gradient estimator holds.
Runtime 0.9 s.

*Unrecorded extra:* the script computes a body-diagonal column the document never shows.
It converges to 3/(2pi) from BELOW (0.4407 at r=1) while the axis converges from ABOVE
(0.5093) — i.e. the sub-leading corrections to G are direction-dependent even though the
leading term is isotropic. Consistent with §9.3 (whose "none" is a leading-order claim),
but this is the data that would quantify sub-leading anisotropy, and it is unused.

**`charge.py` → §8.5, §8.3.** All numbers exact. Part 1 q-multisets
{-1:3,+1:3} / {-2:3,0:6,+2:3} / {-3:1,-1:3,+1:3,+3:1}, all NET = 0. Part 2 both drift
subsets and their values 2/sqrt(3) = 1.154701 and sqrt(3)/2 = 0.866025, both summing to
(2,2,2). Part 3 coupling 0.477465 / 0.238732 / 0.159155 = 6:3:2 as flux/(4pi).
Part 4 sublattice dead end: 8-type FLIPS at 0.3358, and 0.336 != 0.5. Runtime 0.14 s.

**`lightcone.py` → §9.2, §9.3.** Numbers exact: v(n) = 1.000000 / 0.707107 / 0.577350,
anisotropy sqrt(3) = 1.732051, cone shapes (octahedron / cube / cube = L-inf ball), and
covariances 0.333333 / 0.666667 / 1.0 / 0.692308 with off-diagonal **exactly** 0.00e+00
for all four sets. Runtime 0.2 s. See F10 — the docstring is stale.

*Unrecorded extra that strengthens §9.4:* the ballistic table shows max Euclidean
displacement per tick is CONSTANT in t (1.000000, 1.414214, 1.732051, 1.732051 at every
t = 1..6). The cone anisotropy is exactly scale-invariant from the first tick and never
decays, so one cannot argue the sqrt(3) washes out over many ticks. The (ka)^2 dispersion
suppression is therefore doing ALL the work in §9.4 — which makes E4 (board #41) more
load-bearing than the document's phrasing suggests.

**`isotropy.py` → §9.5.** Full reproduction. Direction table 6 / 26 / 98 / 410 / 3458 /
27818 / 222962 with matching angular spacings (82.92 / 39.83 / 20.52 / 10.03 / 3.45 / 1.22
/ 0.43 deg). Moment table +0.400000 / -0.100000 / -0.266667, all off-diagonal exactly 0.
Both mixtures (2/5:3/5 axial:body, 1/5:4/5 axial:face). Hop distribution 0.626046,
+4.3410%, off-diagonal 5.35e-18. Runtime 4.9 s. But see F11.

**`confined.py` → §5.6. THIS IS THE CORRECTED VERSION** — its docstring carries the fix
note ("an earlier version mirrored the coordinate as c-2d ... that was a bug"). The F10
inference did NOT hold for this script. Full reproduction: 3x3x3 (18 transient) E[T] =
5 / 4 / 4, null 1/3 / 1/3 / 3/5, body 2/3 / 2/3 / 2/5, all -44.44% = -4/9. 5x5x5 (90
transient) reflect = 2/5, 61/210, 47/210, 8/105, 1/105 with E[T] = 4 exactly (denominator
210 confirmed); stay 103/25 = 4.12, resample 99/25 = 3.96, both +1.85%; reflect
**-0.1364% = the documented -0.14%**. Runtime 3.2 s.

**`biased.py` → §9.6.** Exact reproduction of every row and column: E[T] = 4.000000 at
beta = 0.0 / 0.4 / 1.0 / 2.0 / 3.0; E[disp_x] 0.000000 / 0.533259 / 1.326756 / 2.517693 /
3.320569; v 0.000000 / 0.133315 / 0.331689 / 0.629423 / 0.830142; P(null) 0.218714 /
0.203783 / 0.143407 / 0.048122 / 0.010444; P0/P0(rest) and 1/gamma columns likewise.
Docstring is clean (stops short of a conclusion — no F10-style staleness).

**Analytic closure (`_check_biased_analytic.py`), no propagation needed.** By Wald,
v = E[disp]/E[T] is exactly the mean step, so v(beta) = 2 sinh(beta)/(2 cosh(beta) + 4).
That closed form reproduces every documented v including **beta = 6 -> 0.990170**, and
saturates at 1 (beta = 40 gives 1.0000000000). It also identifies the beta values behind
§9.6's small-v dilation table: **0.05, 0.1, 0.2**. Rows for 0.1 and 0.2 appear directly in
the main run and match exactly (0.995548 / 0.999444 and 0.982346 / 0.997775); recomputing
the derived columns gives 4.0068 and 3.9722 vs the document's 4.0066 and 3.9723.

**The time-dilation failure is confirmed and the factor of 8 is structural:**
1/gamma = sqrt(1-v^2) -> 1 - v^2/2, so (1 - 1/gamma)/v^2 -> 0.5; measured -> 1 - 4v^2, so
(1 - ratio)/v^2 -> 4.0; ratio 4/0.5 = 8. The hop clock is coordinate time and does not
dilate. This is the verified input for E0.3 (board #36).

**Both remaining rows closed** (`_biased_missing_rows.py`): beta = 0.05 gives
v = 0.016667, P0/P0(rest) = 0.998885, 1/gamma = 0.999861, (1-ratio)/v^2 = **4.0153**,
(1-1/gamma)/v^2 = **0.5000** — every digit of §9.6's small-v table row 1. beta = 6.0 gives
v = 0.990170, P(null) = 0.000036, **E[T] = 4.000000 at 99% of the axial light speed** —
the strongest form of the speed-independence claim. §9.6 is COMPLETELY verified, both
tables, every row.

**`reflbias.py` → §9.7.** Exact: tilt covariance 0.360537/0.319731 = 1.1276,
0.435519/0.282240 = 1.5431, 0.652910/0.173545 = 3.7622, all STRETCHED (wrong sign for a
moving charge). 5^3 wall bias: transmit E[T] = 3.959316, v_max = 0.009596; stay
E[T] = 4.022767, v_max = 0.002851. Both match. Runtime 0.6 s.

**F11/F12 do NOT recur here.** §9.7 works within the single axial shell where all |s| = 1,
so unit-vector and |v|^4 weighting coincide trivially, and the quantity at issue is the
degree-2 covariance anyway. The measure error is specific to MULTI-SHELL mixtures
(§5.6 landing distribution, §9.5 shell mixtures) — two instances of one error, not a habit.

## FINDINGS

### F1 — §6 prose claim refuted (wording, not conclusion)
Document: *"Neither shell is redundant — drop the six and the cell does not close."*
**The cell does close.** Measured (`_check_close.py`, HalfspaceIntersection + ConvexHull):

```
8 body diagonals ONLY     volume = 4.5   8 facets,  6 verts   (an octahedron)
6 axial ONLY              volume = 8.0   6 facets,  8 verts   (a cube)
both (true Voronoi cell)  volume = 4.0  14 facets, 24 verts
```

Both single-shell regions are bounded. What fails is that the diagonal-only octahedron
has volume 4.5 against the required fundamental-domain volume 4, so it overshoots and
cannot tile. **Correct wording: neither shell is redundant because the cell must have
volume 4 to tile, not because it fails to close.** Conclusion survives; sentence must change.

### F2 — new cross-check worth adding to the paper
The true cell's volume is **exactly 4** = the lattice index [Z³ : Λ_BCC] = the same 4 that
Kac's lemma gives as E[T] in §7.1. The Voronoi volume and the mean hop time are the same
integer for the same reason. Independent geometric confirmation of §7.1's headline number,
not currently in the document.

### F3 — provenance tag too strong (minor)
§6's truncated-octahedron ↔ permutohedron identification is tagged `[computed]`, but the
script only shows matching face/vertex counts (14/24), which is necessary not sufficient.
The identification itself is a standard result. Retag `[standard]`.

### F4 — method note (not a defect)
`both68.py::is_facet` uses `linprog` with a 1e-9 tolerance: integer answers via floating
point. Robust at this degeneracy distance, but note it if the script is ported.

### F5 — §5.6 mis-cites §5.4 for the tail figure
§5.6 opens *"The tail of §5.4 (16-26% of hops beyond the 14-neighbourhood)"*. But §5.4's
tail row is **0.08908304**, not 16%. The 16% figure is face + tail
= 0.07357348 + 0.08908304 = **0.16265652**, which is §5.3's "farther" row. The phrase
"beyond the 14-neighbourhood" is itself correct (the 14-neighbourhood is 6 axial + 8 body,
so the face shell is beyond it). The 26% upper end is §7.2's 26-shell row — a different
step set. Fix: cite §5.3, or say "the face shell plus §5.4's tail".

### F6 [RESOLVED] — E[T] = 4 exactness is Kac's lemma, not `branch68.py` (attribution)
**Resolved by `shells_6_12_8.py` Part 1**, which solves the 4-state quotient chain in exact
rationals and returns exactly 4 (and exactly 1 for the body-only set). The document's
"exact" claim is properly earned; the attribution belongs to that script, not to the
truncated forward propagation. Original note follows.

`branch68.py` prints `E[T] = 4.000000`, but that is horizon-truncated and biased low by
~2e-9 (unabsorbed 4e-11 at t=60). The exactness claim rests on **Kac's lemma** in §7.1
(E[T] = [Z³ : Λ_BCC] = 4), which is far stronger — it is what makes the 4 survive bias,
weighting and reflecting confinement. Keep the attribution straight in the write-up; do
not let a reader think the numeric run establishes exactness. Cross-links to F2.

### F7 — §8.2 UNDERclaims; the coefficient is good to six figures, not three
Document: *"Inverse-square with the coefficient 3/(2pi), to three figures."* That is an
artifact of the L=192 box, not the physics. Convergence test (`_check_greens_L.py`,
slice-wise low-memory S(k1)):

```
   L        G(0)     relerr  | best plateau      |diff from 3/(2pi)|
 192   1.50933032  4.65e-03  | 0.477242 (r=20)      2.23e-04
 384   1.51285817  2.33e-03  | 0.477426 (r=30)      3.84e-05
 768   1.51462211  1.16e-03  | 0.477472 (r=45)      7.02e-06
```

Target 3/(2pi) = 0.47746483. At L=768 agreement is **7e-06 — six figures**. G(0) error
halves exactly as L doubles (confirms the O(1/L) parenthetical); the plateau's best radius
tracks r ~ L/17, so the large-r degradation is pure torus wrap-around pushed outward.

**Why it matters:** §8.3 rests the entire coupling-is-1/D argument on this coefficient being
6/(4pi) exactly. At three figures that is suggestive; at six it is settled. Strengthen the
claim in the write-up — this is the one result that bites on the calibration gap rather
than relocating it. An UNDERclaim, not an overclaim.

### F8 — §8.5's universal is FALSE as written (conclusion survives, reasoning must change)
Document: *"Restoring nonzero charge therefore requires an inversion-breaking subset, and
**every such subset carries a net drift**."* Exhaustive over all 255 non-empty subsets of
the 8-body shell (`_check_charge_nogo.py`):

```
inv-breaking                     240
inv-breaking & zero drift          2   <- counterexamples
nonzero charge                   214
nonzero charge & zero drift        0
```

Two inversion-breaking subsets have drift exactly (0,0,0), so the universal is refuted.
`charge.py` only ever tested the single subset {q > 0} per shell — one instance, not the
universal it is quoted to support.

**The fix makes it stronger, not weaker.** Because every stencil component is in
{-1,0,1}, sign(v_i) = v_i, so **q(V) = (1,1,1).V exactly** for all 26 stencil vectors —
the charge functional is just the linear functional dotting with the body diagonal.
Hence, verified for every subset:

> **net charge = (1,1,1) . (vector sum)**

so **charge != 0 <=> the drift has a nonzero (1,1,1) component**. Charge IS the
(1,1,1)-component of the drift. That is an identity, not a survey, and it makes the no-go
airtight. Replace the asserted universal with this identity in the write-up.

### F9 — the two counterexamples ARE T^3_diamond (bridges §8.5 to §10.6 and board #37)
The two inversion-breaking zero-drift subsets are exactly:

```
[(-1,-1,-1), (-1,1,1), (1,-1,1), (1,1,-1)]
[(-1,-1,1), (-1,1,-1), (1,-1,-1), (1,1,1)]
```

— the two tetrahedra, i.e. **the four tetrahedral directions of T^3_diamond**. So the old
diamond architecture sits inside the cubic stencil as a drift-free, charge-free,
inversion-breaking subset of the body shell. This is a direct structural bridge to §10.6's
open "four or three?" question and to board #37 (the diamond degree-4 comparison): the
diamond is not an unrelated rival architecture, it is a distinguished subset of this one.
Worth pursuing — it may be the cleanest way to state what the switch actually changes.

### F10 [ACTION BEFORE COMMIT] — `lightcone.py` docstring carries the RETRACTED falsification claim
The script's docstring ends: *"Terrestrial Lorentz tests bound c-anisotropy at the 1e-18
level; sqrt(3) is a 73% effect. **This is a falsification-level problem, not a calibration
issue.**"*

That is exactly the claim V1.4's **10 Sept correction retracts**, and for exactly this
reason: it compares the single-tick reachable-set anisotropy against a bound derived from
long-wavelength optics — different quantities. The observable anisotropy is (ka)^2
suppressed, making the experiment a CONSTRAINT ON LATTICE SPACING (a < 0.34 fm at 1e-18),
not a refutation.

The computations are unaffected — the sqrt(3) is real, only the interpretation was
withdrawn. **But handoff #32 asks dcl-mathematics to commit these scripts as the
verification record.** Committed as-is, the repo carries a retracted falsification claim in
a docstring where it will later be read as current. Fix the docstring before committing.

**Inference for `confined.py`:** if `lightcone.py` predates the 10 Sept corrections, the
batch plausibly does too — and 10 Sept is also the date of the reflect-convention fix. So
`confined.py` may be the buggy-reflect version. Tell: whether the 5x5x5 reflecting row
gives E[T] = 4 with denominator-210 rationals, or something uglier with reflect
recommended against.

### F11 [MOST CONSEQUENTIAL] — §9.5's "Route to §9.4 option 3" nulls the WRONG moment
Document: *"Tuning shell weights to kill the degree-4 moment pushes the observable
anisotropy from O((a/lambda)^2) to O((a/lambda)^4)."* The moment in §9.5's table is computed
on **unit vectors**. The dispersion's quartic term is

> Sum p(s) (k.s)^4  =  Sum p(s) |s|^4 (k.n)^4

so each shell enters weighted by **|s|^4 = 1 : 4 : 9** (axial : face : body). Different
measure, different nulling mixture. Measured (`_check_deg4_dispersion.py`):

```
 a (axial share)       k=0.05       k=0.01      k=0.002   note
        1.000000    -0.055561    -0.055556    -0.055581   pure axial
        0.400000     0.080793     0.080808     0.080901   2/5 : 3/5  <- sec 9.5 nulling
        0.857143    -0.000009    -0.000001     0.000065   6/7 : 1/7  <- |s|^4-weighted
        0.000000     0.111094     0.111110     0.111119   pure body
```

**The 2/5 : 3/5 mixture makes the observable anisotropy WORSE** (+0.0808) than pure axial
(-0.0556), while nulling the unit-vector moment exactly. The correct dispersion-nulling
mixture solves (a + 3b)/(a + 9b) = 3/5, i.e. **a = 6b -> 6/7 axial : 1/7 body**, confirmed
numerically at ~1e-6 (float floor).

**Consequences:**
1. **Board #43 (E6) as written would solve for the wrong weights** — it says "solve shell
   weights for zero degree-4 moment (§9.5)". The |s|^4 weighting must be written into the
   issue, else the experiment increases the anisotropy it exists to reduce.
2. It partly answers §9.5's own "worth ten minutes" question about the 2/5 : 3/5
   shell-swap near-miss vs the §5.5 parity split: since 2/5 : 3/5 is not the operative
   mixture for the dispersion, the coincidence is much less meaningful than hoped.
3. The unit-vector moment is still the right quantity for the isotropy of the hop
   *landing-direction* distribution — this finding is specific to the dispersion route.

*Free cross-check:* pure axial gives -0.055556 = -1/18 exactly, independently reproducing
§9.4's k^2/18 law from a different code path.

### F12 — §5.6's -0.14% is measure-dependent and does NOT transfer to observables
Same weighting issue as F11. `confined.py::cubic4` normalises to unit vectors, but the hop
process's dispersion anisotropy is weighted by |v|^4, and these landing classes span
|v|^4 = 9 to 144. Measured (`_check_confined_measure.py`; k = 0.02 and 0.005 agree to 4
digits, confirming k^2 scaling):

```
construction               unit dev  |v|^4 dev  unit rel  |v|^4 rel  hop-dispersion aniso/k^2
5x5x5 reflect             -0.000818  +0.003154   -0.14%    +0.53%  -0.002061  -0.002058
5x5x5 stay/resample       +0.011111  -0.020323   +1.85%    -3.39%  +0.014720  +0.014725
3x3x3 (any conv.)         -0.266667  -0.266667  -44.44%   -44.44%  +0.111108  +0.111111
```

**The recommendation SURVIVES** — reflect wins under every measure, so board #39 (E2) can
proceed on the 5x5x5 reflecting cell as planned. But the headline -0.14% becomes **+0.53%**
under the |v|^4 measure (sign flipped, 4x larger) and the direct dispersion anisotropy is
**-0.00206/k^2**. Three numbers for three questions; do not quote -0.14% as an observable.

Margin over stay/resample: 13.6x (unit-vector), 6.4x (|v|^4), 7.1x (direct dispersion).

*Cross-check:* the 3x3x3 hop-dispersion anisotropy +0.111108 matches the pure-body figure
+0.111110 from the F11 run — as it must, since every non-null hop there is a body diagonal.
Two independent code paths agreeing.

*Method note:* my first version of this check excluded the null hop from the dispersion
sum, so E(0) != 0 and every row flattened to ~0. The null must be INCLUDED in the
dispersion (it contributes cos(0) = 1) while remaining EXCLUDED from the direction
moments. Corrected before the numbers above.

### F13 — "thirty times better" attaches the right number to the wrong baseline
§5.6: *"a degree-4 anisotropy of -0.14%, **thirty times better than stay or resample**."*
Stay/resample is +1.85%, so the ratio is **13.6x**, not thirty. The ~31x ratio is against
the **unconfined** case (+4.34%) — a different row of the same comparison table. Fix the
baseline or the multiplier.

### F14 [ACTION BEFORE COMMIT] — `reflbias.py` docstring contradicts its own output
Docstring: *"Both branches change the BCC coset by q(s) -- identical, since -s == s mod 2 --
so the induced quotient walk is untouched and **E[T] = 4 exactly**."* Its own run gives
**3.959316** (transmit) and **4.022767** (stay). The document has this RIGHT in §9.7:
"E[T] is no longer exactly 4: transmission adds an absorption channel and 'stay' adds a
quotient self-loop, so the §5.8 theorem (which covers pure reflection) does not extend to
either." The `r = 1.00, 1.00` row does give exactly 4.000000 — pure reflection, where the
theorem applies — so the theorem's scope is right and only the docstring overreaches.
Fix before committing (same class as F10).

### F15 — §5.8's "any shape" overstates the theorem; the L1-ball row uses an unstated rule
All 15 entries of §5.8's table reproduce (`_check_58_shapes.py`, `_check_l1_variants.py`).
But the reflect rule `s -> -s` is **undefined** on the L1 ball: at p = (-3,0,0) a transverse
step s = (0,1,0) gives |p+s|_1 = 4 AND |p-s|_1 = 4 — when a coordinate is zero, |0+1| =
|0-1|, so the L1 norm rises either way. 72 such cases across 48 transient sites.

Where the rule IS always legal: **exactly the boxes** (products of intervals) — verified
cube R=2 YES, slab YES, asymmetric box 3x5x7 YES; L1 ball NO, L2 ball NO.

The document's L1-ball reflect = 4.000000 **is** reproducible, under "reflect when legal,
stay when the reflection is also illegal" (exact rational 4). That fallback is never stated.
(`refl-else-resample` gives 42/11 = 3.818182.)

**Fix:** §5.6 is already careful ("for unit steps **in a box** this replacement is always
legal"). §5.8 must match it — say "any box", or "any region on which the reflection is
defined", not "any shape" — and state the L1-ball fallback. The argument that the
reflecting cell is *natural rather than tuned* survives for the recommended 5x5x5 cell,
since a cube is a box; only the sweeping generality needs narrowing.

## Open question for the user
Whether PM commits these to dcl-mathematics directly or hands them to that session
under handoff #32 (PM scope normally excludes direct edits to other repos). Unanswered.
