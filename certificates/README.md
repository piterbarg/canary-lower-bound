# Certificate banks

The companion paper

> V. V. Piterbarg, *Universal Minimal Measures for Canary Options: Explicit
> Families and Exact Thresholds*, SSRN working paper, 2026,

proves several threshold statements by finite exact computations: linear
programs with rational data solved exactly, sign checks of polynomials at the
corners of boxes and triangles, and exact decisions of the universality
question at individual slope ratios. Its appendices say, of the four-atom
dilation family and of the binomial flow, that "every certificate is recorded
in the supplementary code". This folder is that record. Statement numbers below
refer to the companion paper (its label is given in parentheses).

A *certificate bank* is a JSON file written by the research code at the end of
a run. Every rational number in it is stored as a string (`"45/64"`,
`"-3/16"`, `"1"`), never as a float, so the files round-trip through
`fractions.Fraction` without loss. Booleans record which checks passed; lists
of records hold the certificates themselves. The files are copied byte for
byte from the author's working directory and are not edited for publication,
so they keep their original names, a few bookkeeping keys (timings, the name of
an internal planning note that is not included, short check labels) and the
index conventions of the code that wrote them, described below.

To load a bank:

```python
import json
from fractions import Fraction as F
bank = json.load(open("certificates/_r22_sandwich.json"))
seg = bank["segments"][0]
a, b = (F(x) for x in seg["dual_value_affine"])   # dual value a + b*r
print(a, b, a + b * F(5, 4))                       # 45/64 -3/16 15/32
```

The scripts that produced the banks are research code and are not included.
The thresholds and verdicts that the banks record (the window boundaries
`c = 3r - 2`, the column `13/12, ..., 1153/1120`, the cells `5/4` and
`37/32`, and the universal / not-universal verdicts at the slope ratios
tested) are recomputed independently, from the marginals alone, by
`python run_all.py` at the repository root; the banks add the pins,
witnesses and dual vectors behind the "only if" halves and the uniqueness
clauses, which `run_all.py` does not reproduce.

## Index of the banks

| File | Statement it backs | What each record certifies | Records |
|---|---|---|---|
| `_r14_g4win.json` | Theorem 3.4 (`thm:g4`), via Lemma A.19 (`lem:dil:g4pins`) | the hold sets that maximize the barrier value are prefixes: for each non-prefix hold set `H` and each cell of the strike plane, exact sign conditions at the corners (and along the recession directions) of the cell showing that a prefix barrier dominates, for every `r >= 1` and `c` in `(1, 3)` | 165 domination certificates, 8 check summaries |
| `_r14_g4win2.json` | Theorem 3.4 (`thm:g4`), Lemma A.19 (F1)-(F3) (`lem:dil:g4pins`), Lemma A.20 (`lem:dil:g4wit`) | the pinned Canaries of (F1) for the first two atoms on the whole window `c <= 3r - 2` and for the other atoms on sub-regions; the aggregate pin of the first three fibers on the serving wedge; 15 triangles of the `(c, r)` plane covering `{3r - 2 < c < 3, r >= 1}` on each of which a maximizer containing the second atom is exhibited (F2); the service modes of the chain coupling at 300 (cell, fiber) combinations; the margin identities of Lemma A.20(b) | 69 pins, 15 triangles, 300 service modes, 3 aggregate regions, 3 discovery instances, 16 check summaries |
| `_r14_perturb.json` | Theorem 3.4 (`thm:g4`), uniqueness clause (supporting evidence) | at six points of the serving wedge, the chain coupling serves and a small admissible perturbation of it (two sizes) does not, with the Canary that detects the perturbation | 12 rows |
| `_r18_window.json` | Theorem 4.2(b) (`thm:column`) at `n = 2, 3`; Lemma B.7 (`lem:fl:pins`), Lemma B.8 (`lem:fl:family`) | the boxes of strikes on which a prefix is a maximizer for every `r >= 1` (the pins), the rank of the pin system and the uniqueness of the pinned coupling, and the exact margins of the moving witness Canary at the breakpoints of `[1, r*(n)]` | 5 + 9 pin boxes, 9 + 9 margins |
| `_r19_column.json` | Theorem 4.2(b) (`thm:column`) at `n = 4, ..., 8`; the triangular pin lattice of Appendix B.4 (`app:flow:cert`) | the same three computations at the five remaining cells of the column; at `n = 7, 8` a second witness family anchored at the maximizing knot | 14 + 20 + 27 + 35 + 44 pin boxes; 45 to 3661 margins per cell |
| `_r21_window.json` | Theorem 4.3 (`thm:window33`), first two steps of the proof | at the cells `(3, 3)` and `(4, 4)`: the pin boxes for the pairs (prefix, atom), the rank of the pin system (22 of 28, 33 of 45) and an entry of the coupling that the pins leave free, the junction and knot attaining `r*`, and the margins of the moving Canary of Lemma B.8 along `[1, r*]` | 9 + 16 pin boxes, 7 + 5 margins |
| `_r21_window2.json` | Theorem 4.3 (`thm:window33`), the remark that searching all hold sets adds no rank | the pin pairs found when every subset `H` (not only prefixes) is tried: 18 and 44 pairs, same ranks 22 and 33, still not a point | 9 + 28 non-prefix pairs |
| `_r21_serve.json` | Theorem 4.3 (`thm:window33`): "the chain coupling serves exactly for `r >= r*`", "the tree" | exact decisions of the universality question at `(3, 3)`, `(4, 4)` and `(5, 5)` below and at `r*`; the composed tree never serves (15 probes) and violates 6 of 18, resp. 12 of 32, pin equalities; the minimum margin over the pinned polytope at nine slope ratios per cell: `(5/4 - r)/16` at `(3, 3)`, zero at `(4, 4)` | 9 decisions, 15 tree probes, 2 pin counts, 18 minimum margins |
| `_r22_sandwich.json` | Theorem 4.3 (`thm:window33`), "no universal measure for `r < 5/4` at `n = 3`"; Lemma B.9 (`lem:fl:sandwich`) | one dual vector, feasible for every `r`, whose value `45/64 - 3r/16` equals `L(c(r)) + (5/4 - r)/16` on the whole segment `[1, 5/4]`, and a coupling attaining it | 1 segment |
| `_r22_canary2.json` | Theorem 4.3 (`thm:window33`), "no universal measure for `r < 101/96` at `n = 4`" | the seven Canaries of the min-max programme, its exact values at nine slope ratios, the location of the edge `101/96` by a two-sided solve, and four consecutive intervals covering `[1, 101/96]` with a fixed dual line `101/768 - r/8` positive on each; an exact decision at the edge | 27 explored anchors, 7 Canaries, 9 + 9 probes, 4 certified pieces |
| `_r22_chart.json` | Theorem 4.3 (`thm:window33`), `n = 4` step (supporting evidence: how the five fixed Canaries with `k2 = -1` were chosen) | for three growing sets of moving Canaries and 21 slope ratios, a coupling in the pinned polytope with zero margin on that set, and a Canary on which it still fails | 21 rows |

Record counts are exact; "margins" counts the number of slope ratios at which
a margin is stored.

## Conventions shared by the files

* Rationals are strings; `true`/`false` are the outcomes of checks; `null`
  means "not applicable" or "not computed".
* **Dilation files** (`_r14_*`): the four atoms of `q1` are indexed `0..3`
  (the paper writes `x_1, ..., x_4`), so the hold set `[1]` is `{x_2}` and
  the prefixes are `[]`, `[0]`, `[0, 1]`, `[0, 1, 2]`, `[0, 1, 2, 3]`. The
  strike plane is cut into cells by the atoms of `q1` (index `t1seg`) and
  the levels of `q2` (index `t2seg`). Polynomials in `(c, r)` are stored as
  strings such as `"21/4 + -6*r + -21/4*c"`.
* **Flow files** (`_r18_*`, `_r19_*`, `_r21_*`, `_r22_*`): the atoms of `q1`
  are the integers `0..n`; a hold set `H` is a list of them, and the
  prefixes are `[0]`, `[0, 1]`, ... . `yb` is an atom of `q2` (an integer
  for two steps, a half-integer for three steps). A pin box is
  `[[k1_lo, k1_hi], [k2_lo, k2_hi]]`, a rectangle of strikes free of atoms
  with `[k2_lo, k2_hi]` inside `(yb - 1, yb)`; `extends_all_r` records the
  sign check of the `r`-coefficients described after Lemma B.7, which
  extends the certificate from `r = 1` to every `r >= 1`. The verdict strings
  `"UMM+"` and `"2UMM-"` of the exact decider mean *universal* and *not
  universal, a conflicting pair of Canaries exhibited*, in the vocabulary of
  `canary/decider.py`.
* **Moving Canary families** are stored as a pair `[X, S]` (or a triple
  `[X, S, k2]`), to be read, by comparison with the midpoint family of
  Lemma B.8, as `k1(r) = (X - S/r)/2` with `X = 2j + 1` and
  `S = C_j(k) + C_{j+1}(k)`; the stored values are consistent with this
  reading at every cell (for instance `X = 1, S = 13/12` at `n = 2`,
  anchored at the junction `0`, knot `0`). `margins` maps each slope ratio
  at which the margin was evaluated (breakpoints, endpoints and one point of
  each piece) to the exact margin
  `V(Q*; k(r)) - L(k(r))`; `margin_at_r1` and `margin_at_rstar` repeat the
  values at the two ends.

## The files one by one

### `_r14_g4win.json`

Top-level keys: `results`, `fails`, `certificates`, `hard`, `total_secs`,
and the planning-note pointer. `results` is a list of eight check summaries
(`check`, `tag`, `ok`, `detail`), among them: the shadow measures of all 16
hold sets are affine in `c` (384 checks), the coupling polytope has
dimension 6 at 26 values of `c`, the chain coupling (Lemma A.17) is
non-negative on `(1, 3)` symbolically, and 165 of 165 domination certificates
hold with an independent grid check of 4875 points finding no strict
non-prefix maximizer. `certificates` holds the 165 records with keys `H`
(the non-prefix hold set), `t1seg`, `t2seg` (the cell), `P` (the dominating
prefix), `P2` and `mu` (a second prefix and a weight: when `P2` is present,
the dominating barrier appears to be the combination of the two prefix
barriers with weight `mu`; 139 records have `P2 = null`, `mu = 1`), and
`conds`, the list of sign conditions, each a pair `[description, sign]` with
`description` either `"corner (t1, t2)"` or `"recession ... >= 0"` and
`sign` one of `"+"`, `"0"`. The exact meaning of the weight `mu` and of the
two recession conditions is not documented in the file beyond these strings.
`fails` and `hard` are empty.

### `_r14_g4win2.json`

Top-level keys: `results` (16 check summaries), `fails` (empty), `pins`,
`adaptive`, `qe`, `modes`, `resistant` (empty), `discovery`, `total_secs`,
and the planning-note pointer.

* `pins` (69 records: `row`, `t2`, `region`, `t1num`, `t1den`) are the
  Canaries of Lemma A.19 (F1): on the level `k2 = t2` the first strike
  `k1 = t1num / t1den` (both polynomials in `c, r`) makes the two prefixes
  ending at atom `row` the only maximizers, on the sub-region `region` of
  the `(c, r)` plane. `region` is one of `"FULL"` (the whole window
  `c <= 3r - 2`) or a half-plane `"c<=4r-3"`, `"c<=2r-1"`, `"c<=3/2r-1/2"`,
  `"c<=3r-2"`, `"c<=5r-4"`. Rows 0 and 1 carry pins with `region = "FULL"`,
  which is what Lemma A.20(a) uses; rows 2 and 3 are pinned on sub-regions
  only.
* `adaptive.agg` records the aggregate pin of the first three fibers on
  the serving wedge: `ok`, `regions` (3) and `leaves`, the three regions.
  The leaves are stored as lists of strings mixing plain numbers and
  Python `Fraction(...)` representations; which entries are the vertices
  and which are flags is not documented, and we do not interpret them.
* `qe` (`ok`, `triangles` = 15, `leaves`) is Lemma A.19 (F2): each leaf is a
  triangle `tri` of the `(c, r)` plane (three vertices `[c, r]`), a
  parameter `mu` and the maximizer `Hstar` that contains atom `1`
  (`[0, 1, 2]`, the prefix `P_3`, or `[0, 1, 2, 3]`, the whole support).
  `mu = 1` corresponds to the Canary `k^V` displayed after Lemma A.19; the other values
  (`1/4`, `1/2`, `-1/4`) select the "fixed rational convex combinations of
  `k^V` and `k^Z`" of (F2). The precise parametrization in `mu` is not
  written in the file.
* `modes` (300 records: `P`, `t1seg`, `t2seg`, `fiber`, `kind`, `mus`) are
  the service modes of the chain coupling: for the prefix `P` on the cell
  `(t1seg, t2seg)` and the fiber `fiber`, whether the hold/exercise
  consistency of the chain coupling is certified on the whole window
  (`kind = "window"`, 254 records) or on the serving wedge only
  (`kind = "wedge"`, 46 records). `mus` lists weights used by the
  certificate; its meaning is not documented.
* `discovery` records three instances above the boundary (`c = 3/2,
  r = 21/20`; `c = 2, r = 5/4`; `c = 5/2, r = 11/10`) at which the chain
  coupling fails, with the gap, the Canary `(t1, t2)`, the maximizers and
  the per-atom mode (hold versus exercise) at that Canary.
* `results` includes the identities of Lemma A.20(b) in the file's units
  (`4(c - 3r + 2)`, `3(c - 4r + 4)`), the slice identity `3(4r - 3 - c)`,
  and the summary "rows 0, 1 and the aggregate pinned".

### `_r14_perturb.json`

Keys `rows` (12 records) and `fails` (empty). Each row is a point `(c, r)`
on the serving wedge (six points, each with two perturbation sizes
`d3 = 1/1000` and `1/100000`): `wedge` and `chain_serves` are true, the
perturbed coupling is still non-negative (`pert_nonneg`) but does not serve
(`pert_serves = false`), and `witness` is the Canary `(k1, k2)` on which it
fails. The direction of the perturbation in the coupling polytope is not
recorded in the file.

### `_r18_window.json`

Keys `cell_2`, `cell_3` (the cells `n = 2, 3` of the two-step flow) and
`parametric_pin_standalone`. Each cell has: `n`, `rstar`; `certs`, the pin
boxes (`H`, `yb`, `box`, `extends_all_r`), 5 at `n = 2` and 9 at `n = 3`,
which are exactly the pairs (proper prefix, atom from the diagonal up);
`all_certs_extend_all_r`; `rank` of the pin system (14 and 21, the law
`(n^2 + 9n + 6)/2`); `unique_fixed_pins` and `unique` (the pinned polytope is
the single point `Q*`), `unique_witness` (`null`, no free entry found);
`family_X_S` and `family_witness` (the moving Canary of Lemma B.8 and its
margins: `1/96` and `1/64` at `r = 1`, zero at `r*`); `window_theorem`,
`covers_r1` (true: the argument covers `r = 1`), `r1_via_decide` (false: no
exact decision was needed at `r = 1`). `parametric_pin_used = false` and
`parametric_pin_ok = null` record that an alternative device, pins depending
on `r`, was not needed; `parametric_pin_standalone` (59 grid points, none
bad) is a stand-alone check of that device, not used in the proof.

### `_r19_column.json`

Keys `cells` (a map `"4"`, ..., `"8"`), `stopped_after` (8) and
`T1_ratified` (a flag recording that the `n = 4` cell reproduced figures
obtained in an earlier independent run). Each cell carries the same keys as
in `_r18_window.json` plus `wall_seconds`; the pin boxes number 14, 20, 27,
35, 44 (the triangular lattice `(n+1)(n+2)/2 - 1`), the ranks 29, 38, 48,
59, 71, uniqueness holds at every cell, and the margins of the moving Canary
anchored at the knot `0` of junction `n - 2` are `1/128, 1/256, 1/512,
1/1024, 1/2048` at `r = 1`. At `n = 7, 8` that family has
`ok_open_interval = false` (the knot `0` cannot reach past
`1 + 1/(4(n+1)) < r*(n)`), and the additional key `family2` (`jstar`,
`kstar`, `X_S`, `ok_open_interval`, `margins`, `margin_at_r1`,
`margin_at_rstar`) holds the family anchored at the maximizing knots
`(2, 2)` and `(3, 3)`, with margins `7/2560` and `33/10240` at `r = 1` and
zero at `r*`, as in the proof of Theorem 4.2(b).

### `_r21_window.json`

Keys `cell_3`, `cell_4` (the cells `(3, 3)` and `(4, 4)` of the flow with
`kappa = n`) and the planning-note pointer. Each cell has `n`, `k`
(= `kappa`), `rstar` (`5/4`, `37/32`), `argmax` (the junction and knot
attaining `r*`: `(2, 3/2)` and `(3, 3)`), `certs` (9 and 16 pin boxes, same
format as above), `n_certs`, `all_certs_extend_all_r`, `full_grid_count`
(18 and 32: the number of pin equalities, two per pair), `k2_lattice_count`
(9 and 14; its meaning is not inferable from the file), `rank` and `ncol`
(22 of 28, 33 of 45), `unique_fixed_pins = false`, `unique_witness` (a list
of four numbers: an entry index of the coupling array and three values,
evidently the minimum and maximum of that entry over the pinned polytope and
the chain coupling's value, in an order the file does not document),
`family_X_S` and `family_witness` (the moving Canary of Lemma B.8 anchored at
`argmax`: margins `1/64` and `5/1024` at `r = 1`, zero at `r*`), and the
flags `window_theorem = false`, `covers_r1 = false`, which record that at
these cells pins plus a single moving Canary do not by themselves close the
window; `_r22_sandwich.json` and `_r22_canary2.json` do.

### `_r21_window2.json`

Same two cells, when the search for pin pairs is extended from prefixes to
all subsets `H`: `n_certs` (18 and 44), `n_prefix` (9 and 16), `n_nonprefix`
(9 and 28), `nonprefix_list` (the extra pairs `H`, `yb`), and the unchanged
`rank`, `ncol`, `unique = false`, `unique_witness`. This is the sentence
"extending the search to all subsets adds certificates but no rank" of the
proof of Theorem 4.3.

### `_r21_serve.json`

Keys `bank` and `ok`. `bank` has five parts.

* `T4a` (9 records: `n`, `r`, `want`, `got`, `gap`, `pair`): exact decisions
  at the cells `(3, 3)`, `(4, 4)` and `(5, 5)` (the last with
  `r* = 681/560`, not stated in the paper) at `r = 1`, at one interior
  `r < r*`, and at `r*`: not universal below `r*` with the conflicting pair
  of Canaries and its exact portfolio gap, universal at `r*`.
* `T4b` (15 records: `n`, `r`, `below_rstar`, `served`): the composed tree
  serves at none of five slope ratios per cell, below and above `r*`.
* `T4c` (2 records: `n`, `violated_pin_rows`, `n_pin_rows`): the tree
  violates 6 of the 18 pin equalities at `(3, 3)` and 12 of the 32 at
  `(4, 4)`.
* `minmargin` (18 records: `n`, `r`, `min_margin`, `positive`,
  `matches_affine_profile`): the minimum over the pinned polytope of the
  margin at the moving Canary, at nine slope ratios per cell, several of
  them off the dyadic grid: at `(3, 3)` positive and equal to `(5/4 - r)/16`
  at every probe; at `(4, 4)` zero at every probe (the single moving Canary
  does not suffice there).
* `decide_seconds`: wall-clock timings of the decisions.

### `_r22_sandwich.json`

The smallest bank, and the certificate of the sentence "no universal measure
for `r < 5/4` at `n = 3`". Keys: `cell` (3), `rstar`, `k1_range`
(`[5/4, 3/2]`, the range of `k1(r)` on `[1, 5/4]`, crossing no atom),
`L_segments` (one segment `[1, 5/4]` on which the lower bound `L(c(r))` is a
single barrier line), `segments` (one record: `ra`, `rb` the ends,
`L_index`, `strong_duality_at_mid`, `dual_feasible_rfree` (the dual vector
is feasible independently of `r`, Lemma B.9(a)), `dual_value_affine`
(`45/64 - 3r/16`), `target_affine` (the same line, i.e. `L(c(r)) +
(5/4 - r)/16`), `dual_identity`, `primal_coverage` (a coupling in the
pinned polytope attains the bound), `extra_extractions` (empty),
`seg_ok`), and the verdicts `sandwich_certified` and `window33`. The dual
vector itself is not stored, only its value line and the feasibility flag.

### `_r22_canary2.json`

The certificate of "no universal measure for `r < 101/96` at `n = 4`".
Keys: `cell` (4), `rstar` (`37/32`), `anchor1` (`(3, 3)`, the anchor of the
first moving Canary), `exploration` (27 records `j`, `k`, `all_positive`,
`n_probed`, `worst`: candidate anchors for a second moving Canary, none of
which makes the min-margin positive at all probes), `registered_second`
(`null`) and `hf2_falsified_at_probes` (true: a single second moving Canary
does not suffice), `battery` (the seven Canaries of the proof as triples
`[X, S, k2]`: the family anchored at `(3, 3)`, the family anchored at
`(1, 1)`, and the five fixed Canaries `k1 = -1, -7/8, -3/4, -5/8, -1/2` on
the level `k2 = -1`), `battery_probes` and `sform_probe` (the exact min-max
value at nine slope ratios, positive at `1`, `641/640`, `133/128` and zero
from `69/64` on; the second list is the same programme in a second encoding
with the same zero pattern, `sform_same_zero_pattern`), `edge` (`101/96`),
`edge_steps` (the two-sided solve locating the edge in three steps),
`minmax_at_edge` (`0`), `breakpoints` (`1, 769/768, 337/326, 24/23,
101/96`), `pieces` (four records `lo`, `hi`, `depth`, `sigma`, `sigma_lo`,
`sigma_hi`, `minmax_at_mid`, `certified`: on each interval a fixed dual
vector has the affine value `sigma = 101/768 - r/8`, positive on the
interval except at the right end of the last one), `leftwedge` (true),
`wedge_profile_line` and `wedge_same_dual_line` (the same dual line serves
on all four pieces), `wedge_profile_upper_leg` and `wedge_profile_theorem`
(false: that the min-max *equals* this line was not certified and is not
claimed), and `edge_decide` (an exact decision at `r = 101/96`: not
universal). As in the previous file, the dual vectors are not stored, only
their value lines.

### `_r22_chart.json`

Key `rows` (21 records) and the planning-note pointer. Each row has `pair`
(the set of moving Canaries imposed: the family anchored at the maximizing
knot alone, or together with the family anchored at `(1, 1)` or at
`(2, 0)`), `r`, `survivor_served` (false in every row: a coupling in the
pinned polytope with zero margin on that set is still not universal),
`witness` (a Canary `(k1, k2)` on which it fails; 19 of the 21 lie on the
level `k2 = -1` with `k1` in `[-161/160, -2/5)`, the other two at
`k2 = -37/72` and `-8/5`), and `S_fit_per_junction` (`j0`..`j3`: a fitted
quantity per junction whose definition is not recoverable from the file).
This chart is how the five fixed Canaries on the level `k2 = -1` in the
`(4, 4)` argument were chosen; the proof itself rests on
`_r22_canary2.json`.
