# dcl-cubic-taxicab

Verification scripts for the **cubic taxicab sphere** — the 3×3×3 stencil /
BCC 14-neighbourhood architecture proposed as a replacement for `T³_diamond`.

Source document: `dcl-mathematics` `notes/taxispacehops V1.4.md` (committed 5141d4a).
Figure: `dcl-mathematics` `figures/Cubic-Space.drawio` (403227e).
Tracked on project 6 as epic **#31** (board 030) with children **#32–#44**.

> **Status: validation, not adoption.** Nothing here is settled physics. The
> architecture is under an ordered experimental program (E0.1 … E7) that can
> still kill it — see board #40 (E3), the decision point for whether the geometry
> can carry relativistic kinematics at all.

## Why this repo exists

V1.4 §11.1 lists fourteen results as "already closed by exact computation; do not
re-derive." As of 2026-09-10 none of the scripts that produced them existed
anywhere on disk, so every one of those results was unreproducible. Board **#32**
recovered them. This repo is their home.

`dcl-mathematics` is deliberately paper + Lean only — its `CLAUDE.md` scopes the
no-code rule to Python numerics — so these live here instead, alongside the paper
they verify (author decision, 2026-09-11).

## Layout

```
src/utilities/     the eleven Appendix A scripts, as recovered
src/verification/  independent cross-checks written during the verification pass
VERIFICATION-LOG.md  what reproduced, what did not, and all fifteen findings
```

Run anything with the canonical interpreter (`C:\Users\jackd\.venv-win`,
Python 3.14, numpy + scipy). No script takes arguments.

## What each script establishes

| script | covers | verdict |
|---|---|---|
| `hop_fpd.py` | §5.1 rule-A first passage; the multinomial selection rule | exact |
| `branch68.py` | §5.3–5.5 first return, landing distribution | exact |
| `confined.py` | §5.6 confined cells, three boundary conventions | exact |
| `both68.py` | §6 BCC Voronoi cell + permutohedron | exact |
| `shells_6_12_8.py` | §7 Kac's lemma, branching, all three Bravais cells | exact |
| `greens.py` | §8.2 lattice Green's function, the 3/(2π) coefficient | exact |
| `charge.py` | §8.5 step-linear charge no-go; §8.3 coupling = 1/D | exact |
| `lightcone.py` | §9.2 cone shape and √3; §9.3 covariance isotropy | exact |
| `isotropy.py` | §9.5 direction counting, degree-2/4 moments | exact |
| `biased.py` | §9.6 exponential tilt, the time-dilation failure | exact |
| `reflbias.py` | §9.7 tilt covariance sign, wall-bias momentum capacity | exact |

**Every numerical result in all eleven scripts reproduced exactly.** All fifteen
findings are in the surrounding layer — prose, provenance tags, cross-references,
docstrings, and two moment-weighting errors. See `VERIFICATION-LOG.md`.

## The findings that change work, not just wording

- **F8** — §8.5's no-go has a false middle step ("every inversion-breaking subset
  carries a net drift"; two counterexamples). The conclusion survives and gets
  stronger: `q(V) = (1,1,1)·V` exactly, so charge **is** the (1,1,1)-component of
  the drift. An identity, not a survey.
- **F11** — §9.5's route to relaxing the spacing bound nulls the **wrong moment**.
  The dispersion's quartic term weights each shell by `|s|⁴` (1 : 4 : 9), so the
  nulling mixture is **6/7 axial : 1/7 body**, not 2/5 : 3/5. The documented
  mixture makes the observable anisotropy *worse* than pure axial. Posted to #43.
- **F15** — §5.8's theorem is stated for "a confining region of any shape," but
  the reflect rule `s → −s` is well-defined exactly on **boxes**; it is undefined
  on the L1 and L2 balls. The recommended 5×5×5 cell is a box, so the
  recommendation stands; only the generality needs narrowing.

Two findings make claims **stronger**: **F7** (the Green's coefficient is good to
six figures, not three — which settles §8.3's coupling-is-1/D) and **F9** (the
two drift-free inversion-breaking subsets of the body shell *are* `T³_diamond`'s
four tetrahedral directions — the old architecture sits inside this one).

## Corrections applied to recovered scripts

Two docstrings were stale and were fixed on landing; no numbers changed.

- `lightcone.py` asserted falsification at 1e-18 by Michelson–Morley-type data.
  **Retracted** — it compared a single-tick reachable set against a
  long-wavelength optics bound. The observable anisotropy is `(ka)²`-suppressed,
  making the experiment a *derived bound* on lattice spacing (a < 0.34 fm).
- `reflbias.py` asserted `E[T] = 4` exactly for both wall models. Its own output
  refutes that (3.959316 and 4.022767); the §5.8 theorem covers pure reflection
  only. The `r = 1.00/1.00` row does return exactly 4.

## Not here yet

No paper scaffold. `dcl-paper-experiment-template` still carries the defects its
open fix handoff (`2026-07-22-paper-template-build-discipline`) describes —
`build_dir := build`, `stage_dir := stage`, and a `VENV` pointing at a deleted
venv — so deriving from it now would import three known bugs. The scaffold lands
when that handoff is consumed and the paper is actually being written, which is
gated on the E-program returning verdicts.
