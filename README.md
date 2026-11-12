# canary-lower-bound

Exact computations behind

> V. V. Piterbarg, *Robust Lower Bound for a Canary Call*, SSRN working paper
> [5510958](https://papers.ssrn.com/abstract=5510958), revised version, November 2026,

and its companion paper

> V. V. Piterbarg, *Universal Minimal Measures for Canary Options: Explicit
> Families and Exact Thresholds*, SSRN working paper
> [7688081](https://ssrn.com/abstract=7688081), November 2026.

A Canary is a Bermudan option with two exercise dates. The first paper derives
the sharp model-independent lower bound on its price given the two marginal
distributions, and asks when one martingale measure attains the lower bounds
of all Canaries of a given pair of exercise slopes at once (a *universal
minimal measure*). This repository recomputes, in exact rational arithmetic,
every number in the two papers that is the result of a computation, and
regenerates the generated figures. Section and theorem numbers below refer to
the SSRN version of November 2026 of the first paper.

## Run

```
pip install -r requirements.txt     # numpy + matplotlib, for the figures only
python run_all.py                   # all checks (about 90 s) + all figures
python run_all.py --checks          # numbers only; exit code 0 iff all pass
python run_all.py --figures         # figures only, into figures/output/
```

The checks need nothing beyond the Python standard library. Every quantity is
a `fractions.Fraction`; floats never enter a decision. The code refuses float
input and refuses to return an LP value whose optimality it has not certified.

## What is computed

| Module | Content | Paper |
|---|---|---|
| `canary/exact.py` | finite measures, piece-wise linear functions, the **shadow measure** through the convex hull (with an exact certificate), a small exact simplex with certified optimality | Prop. 5.6 |
| `canary/lower_bound.py` | the sharp lower bound `L(k1, k2)` of a single Canary by the closed form over hold sets (the decision route) and by the dual's finite candidate set (an independent cross-check) | (5.6), Sect. 6 |
| `canary/decider.py` | **is there a universal minimal measure** for `(q1, q2, r)`? Portfolio LP, exact certificate `V(Q;.) == L(.)` on all of R², cutting-plane iteration, minimal conflicting set with its exact gap | Def. 8.1, Sect. 14 |
| `canary/threshold.py` | the greedy fibers and the marginals-only **threshold r\*** (n+1 convex hulls, no LP), the explicit families (dilations, binomial flow, convolutions) and the instances of the papers | Lemma 9.1, Thm. 9.4, Prop. 9.7 |

## What is checked (`checks/numcheck.py`, 51 checks)

| Block | Statement checked | Route |
|---|---|---|
| X | the counterexample of **Proposition 12.1**: two Canaries with no common minimal measure, `L_a = 1/2`, `L_b = 5/16`, `L_port = 5/6`, gap `1/48`, by a hand-checkable two-parameter description of all couplings; the decider returns the quoted pair and gap in two iterations, and the quoted pairs and gaps at dilation factors 3/2 and 5/2 and at `r = 21/20` | `checks/exhibit_certificate.py`, `checks/exhibit_rows.py` |
| H | the instance after the **Helly reduction** theorem: three Canaries on a pair of marginals with a two-dimensional polytope of couplings, every two of which have a common minimal measure while the three have none; lower bounds `63/16`, `57/16`, `45/64`, minimal portfolio value `14191/1728`, gap `1/108`. The paper prints the data and the result; the script is the derivation (the three segments of minimal measures in the polygon, their crossing points, the exact portfolio linear programs) | `checks/helly_instance.py` |
| R | every printed threshold constant: `(c+3)/4` and `(c+2)/3` for the three- and four-atom dilation families, `max(1,(c-1)/2)` for the uniform three-point law, `137/72` and `49/36` for Proposition 10.1, `16/15` for the exceptional convolution instance, the binomial-flow column `13/12, 17/16, 21/20, 25/24, 29/28, 31/30, 1153/1120` and the `(3,3)`, `(4,4)` cells `5/4`, `37/32` | `canary.threshold.selftest` |
| W | **Theorem 11.1** as an if-and-only-if: the three-atom family is universal at `c = 4r-3` and not at `c = 4r-3+1/50`; the four-atom family likewise at `c = 3r-2`, and universal on all of `(1,3)` for `r >= 5/3` | decider |
| C | **Theorem 11.2**: the binomial flow is universal at `r = r*(n)` and not at `r = 1` or `r*(n) - 1/1000`, `n = 2..5`; the `(3,3)` and `(4,4)` cells | decider |
| N | **Proposition 10.1**: instances that are universal at `r = 1`, not universal at several `r` in `(1, r*)`, and universal again at `r*`; the exceptional convolution instance universal at `r = 1` with `r* = 16/15` | decider |

A verdict is one of `universal` (a certified universal coupling is returned),
`not universal (conflicting pair)` (the pair and its exact portfolio gap),
`not universal (conflicting set)`, or `undecided` (iteration budget exhausted;
this never occurs on the instances of the papers).

## Figures (`figures/`)

| Script | Output | Paper |
|---|---|---|
| `shadow_and_chord.py` | `fig_shadow_hull_01.png`, `fig_chord_01.png` | Figures 2 and 7 |
| `windows_and_column.py` | `fig_windows_01.png`, `fig_column_01.png` | Figures 8 and 9 |
| `dens_smoothed.py` | `fig_dens_smoothed_01.png` | Figure 10 |
| `dual_step4.py` | `fig_step4_01.png` | Figure 11 |

Every plotted number is an exact rational; matplotlib only renders.

## A minimal example

```python
from fractions import Fraction as F
from canary.decider import UniversalityDecider
from canary.threshold import rstar, binom_centred, dilate

q1 = binom_centred(2)          # 1/4, 1/2, 1/4 on {-1, 0, 1}
q2 = dilate(q1, 2)             # the same masses on {-2, 0, 2}
print(rstar(q1, q2))           # 5/4 = (c+3)/4 at c = 2
d = UniversalityDecider(q1, q2, F(1)).decide()
print(d['verdict'], d['pair'], d['gap'])
# not universal (conflicting pair) [['-1/3', '0'], ['0', '1']] 1/48
d = UniversalityDecider(q1, q2, F(5, 4)).decide()
print(d['verdict'])            # universal
```

## Certificates

The companion paper proves its "only if" and uniqueness clauses (the
four-atom dilation window, the binomial-flow column `n = 2..8`, the cells
`(3,3)` and `(4,4)`) by finite exact computations whose certificates it
says are "recorded in the supplementary code". They are: `certificates/`
holds the eleven JSON certificate banks (exact rationals as strings) and a
[README](certificates/README.md) that maps each file to the theorem or lemma
it backs and describes its records.

## Provenance

The code was written by Anthropic's Claude models under the author's
direction, as described in the statement on the use of artificial
intelligence at the end of the first paper, and was used to verify the
papers' numbers before publication. This repository contains the kernel on
which the printed numbers and figures depend; exploratory research code is
not included.

## Licence

MIT, see `LICENSE`.
