"""Does a pair of marginals q1 <=cx q2 admit a universal minimal measure at
slope ratio r?  (Definition 8.1 of "Robust Lower Bound for a Canary Call",
SSRN version of November 2026; the algorithm is the one summarised in its
Section 15, "Universality".)

A coupling Q is universal iff V(Q; k) = L(k) for every Canary k = (k1, k2),
where V(Q; k) is the value of the Canary under Q and L(k) its sharp lower
bound.  The decision has three parts.

1. The portfolio LP (portfolio_lp): for a finite set of Canaries, the
   coupling minimising the sum of their values.  If the minimum exceeds the
   sum of the individual lower bounds, no universal measure exists and the
   set is a certificate.

2. The exact certificate (certify): given Q, decide V(Q; .) == L(.) on ALL
   of R^2 without sampling.  Outside a bounded box the identity holds
   structurally or by exact threshold inequalities (see certify); inside the
   box, on each rectangle of the grid supp q1 x supp q2, V and L are both
   piece-wise affine with kinks along finitely many known lines, so V - L
   vanishes identically iff it vanishes at every corner, edge-line and
   line-line intersection.  Any V < L would violate weak duality and raises
   InvariantError.

3. The cutting-plane iteration (decide): solve the portfolio LP on the
   Canaries found so far; if its value exceeds the sum of lower bounds,
   extract a minimal conflicting subset (in every instance of the paper a
   pair) with its exact gap; otherwise certify the solution on R^2 and, if
   it fails, add the Canary it fails on and repeat.

Verdicts: UNIVERSAL (a certified universal coupling is returned),
NOT_UNIVERSAL_PAIR (a conflicting pair and its exact gap), NOT_UNIVERSAL_SET
(a conflicting set of size >= 3 whose pairs all agree), UNDECIDED (the
iteration budget was exhausted; never silently treated as a verdict).
"""

from fractions import Fraction as F
from itertools import combinations

from .exact import (InvariantError, LPError, ScopeError, ambient_nodes,
                    cpot, cx_leq, measure, solve_lp, supp)
from .lower_bound import LowerBoundOracle

UNIVERSAL = 'universal'
NOT_UNIVERSAL_PAIR = 'not universal (conflicting pair)'
NOT_UNIVERSAL_SET = 'not universal (conflicting set)'
UNDECIDED = 'undecided'


# affine forms on a rectangle: f(k1,k2) = A*k1 + B*k2 + C
def _aff_eval(f, k1, k2):
    return f[0] * k1 + f[1] * k2 + f[2]


def _aff_sub(f, g):
    return (f[0] - g[0], f[1] - g[1], f[2] - g[2])


class UniversalityDecider:

    def __init__(self, q1_pairs, q2_pairs, r):
        if isinstance(r, float):
            raise TypeError('float r rejected')
        self.q1 = measure(q1_pairs)
        self.q2 = measure(q2_pairs)
        self.r = F(r)
        if self.r < 1:
            raise ScopeError('r >= 1 only (the paper''s standing convention)')
        nodes = ambient_nodes(self.q1, self.q2)
        ok, why = cx_leq(self.q1, self.q2, nodes)
        if not ok:
            raise ScopeError('q1 <=cx q2 fails: %s' % why)
        self.xs, self.ps = supp(self.q1), [m for _, m in self.q1]
        self.ys, self.ws = supp(self.q2), [m for _, m in self.q2]
        self.n1, self.n2 = len(self.xs), len(self.ys)
        self.oracle = LowerBoundOracle(self.q1, self.q2, self.r)
        # dimension of the polytope of martingale couplings:
        # n1*n2 minus the exact rank of the equality system
        self.dimM = self.n1 * self.n2 - self._eq_rank()

    # ------------------------------------------------------------ misc
    def _eq_rows(self, nvar_extra=0):
        nQ = self.n1 * self.n2
        rows, rhs = [], []
        for i in range(self.n1):
            a = [F(0)] * (nQ + nvar_extra)
            for j in range(self.n2):
                a[i * self.n2 + j] = F(1)
            rows.append(a); rhs.append(self.ps[i])
            a2 = [F(0)] * (nQ + nvar_extra)
            for j in range(self.n2):
                a2[i * self.n2 + j] = self.ys[j] - self.xs[i]
            rows.append(a2); rhs.append(F(0))
        for j in range(self.n2):
            a = [F(0)] * (nQ + nvar_extra)
            for i in range(self.n1):
                a[i * self.n2 + j] = F(1)
            rows.append(a); rhs.append(self.ws[j])
        return rows, rhs

    def _eq_rank(self):
        rows, _ = self._eq_rows()
        M = [list(rw) for rw in rows]
        rank, ncol = 0, len(M[0])
        for col in range(ncol):
            piv = next((i for i in range(rank, len(M)) if M[i][col] != 0),
                       None)
            if piv is None:
                continue
            M[rank], M[piv] = M[piv], M[rank]
            pv = M[rank][col]
            M[rank] = [v / pv for v in M[rank]]
            for i in range(len(M)):
                if i != rank and M[i][col] != 0:
                    f = M[i][col]
                    M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
            rank += 1
        return rank

    def V(self, Q, k1, k2):
        """the Canary value under the coupling Q: per fiber the larger of
        exercising at the first date and continuing."""
        k1, k2 = F(k1), F(k2)
        tot = F(0)
        for i in range(self.n1):
            ai = self.r * self.ps[i] * max(self.xs[i] - k1, F(0))
            ci = sum(Q[i][j] * (self.ys[j] - k2)
                     for j in range(self.n2) if self.ys[j] > k2)
            tot += max(ai, ci)
        return tot

    # ------------------------------------------------- the seed set
    def seed_canaries(self):
        """an initial Canary set: for each k2 at an atom of q2 (and one
        below all of them) the strike k1 at which r C_q1(k1) = C_q2(k2),
        plus one Canary below both supports."""
        K2 = list(self.ys) + [self.ys[0] - 1]
        out = []
        for k2 in K2:
            C = cpot(self.q2, k2)
            if C > 0:
                out.append((self._pin_k1(C), F(k2)))
        out += [(F(self.xs[0]) - 1, F(self.ys[0]) - 1)]
        return out

    def _pin_k1(self, C):
        m = sum(self.ps)
        lo, hi = F(self.xs[0]) - 1 - C / (self.r * m), F(self.xs[-1])
        pts = [F(x) for x in self.xs]
        f = lambda k: self.r * sum(p * (x - k) for x, p in
                                   zip(self.xs, self.ps) if x > k)
        cells = [lo] + pts
        for a, b in zip(cells, cells[1:]):
            va, vb = f(a), f(b)
            if va >= C >= vb:
                return a if va == vb else \
                    a + (va - C) * (b - a) / (va - vb)
        return hi

    def portfolio_lp(self, canaries):
        """min over couplings Q of sum_c V(Q; c); returns (Q, value)."""
        nQ = self.n1 * self.n2
        nT = len(canaries) * self.n1
        c = [F(0)] * nQ + [F(1)] * nT
        A_eq, b_eq = self._eq_rows(nT)
        A_ub, b_ub = [], []
        for ci, (k1, k2) in enumerate(canaries):
            for i in range(self.n1):
                ai = self.r * self.ps[i] * max(self.xs[i] - k1, F(0))
                t = nQ + ci * self.n1 + i
                row = [F(0)] * (nQ + nT); row[t] = F(-1)
                A_ub.append(row); b_ub.append(-ai)
                row = [F(0)] * (nQ + nT); row[t] = F(-1)
                for j in range(self.n2):
                    row[i * self.n2 + j] = max(self.ys[j] - k2, F(0))
                A_ub.append(row); b_ub.append(F(0))
        x, val = solve_lp(c, A_eq, b_eq, A_ub, b_ub,
                          [None] * (nQ + nT))
        Q = [[x[i * self.n2 + j] for j in range(self.n2)]
             for i in range(self.n1)]
        return Q, val

    def pair_gap(self, ca, cb):
        """exact two-Canary portfolio gap L_port - L(a) - L(b), each term a
        full coupling LP."""
        Qa, La = self.portfolio_lp([ca])
        Qb, Lb = self.portfolio_lp([cb])
        Qp, P = self.portfolio_lp([ca, cb])
        return P - La - Lb

    # ---------------------------------------------- the certificate
    def certify(self, Q):
        """exact decision of V(Q; .) == L(.) on R^2.  Returns (True, None)
        or (False, witness_canary).  Raises InvariantError on any V < L
        (a weak-duality breach, impossible for a martingale coupling).

        Coverage of R^2:
          - k2 >= max supp q2 or k1 >= max supp q1: V == L structurally;
          - k1 <= min supp q1 and k2 <= min supp q2: both equal
            sum_i p_i max(r (x_i-k1), x_i-k2) for every coupling;
          - the W tail (k1 below a solved threshold G1lo) and the S tail
            (k2 below a solved threshold G2lo): dominance inequalities
            checked exactly at the thresholds carry V == L in closed form;
          - the box [G1lo, max X] x [G2lo, max Y]: per grid rectangle,
            vertex enumeration of the piece-wise affine V - L."""
        r, xs, ps, ys = self.r, self.xs, self.ps, self.ys
        if any(isinstance(Q[i][j], float) for i in range(self.n1)
               for j in range(self.n2)):
            raise TypeError('float Q rejected -- exact rationals only')
        for i in range(self.n1):
            if sum(Q[i]) != ps[i] or \
               sum(Q[i][j] * ys[j] for j in range(self.n2)) \
               != ps[i] * xs[i]:
                raise InvariantError('Q row %d mass/mean off' % i)
        colsum_ok = all(sum(Q[i][j] for i in range(self.n1)) == self.ws[j]
                        for j in range(self.n2))
        if not colsum_ok:
            raise InvariantError('Q column sums off')
        if any(Q[i][j] < 0 for i in range(self.n1)
               for j in range(self.n2)):
            raise InvariantError('Q negative entry')
        # the within-box argument needs every shadow supported on supp q2;
        # the shadow certificate guarantees it, re-checked here
        yset = set(ys)
        for Hf, SH in self.oracle.barrier_forms():
            if any(z not in yset for z, _ in SH):
                raise InvariantError('shadow S_%s has off-grid support'
                                     % sorted(Hf))
        # solved thresholds for the W and S tails (the "-1" margins and the
        # xs[0]/ys[0] seeds guarantee G1lo < min X and G2lo < min Y, which
        # the tail arguments rest on)
        contmin = [sum(Q[i][j] * (ys[j] - ys[0])
                       for j in range(self.n2)) for i in range(self.n1)]
        psimin = {Hf: cpot(SH, ys[0])
                  for Hf, SH in self.oracle.barrier_forms() if Hf}
        cand = [F(xs[0])]
        for i in range(self.n1):
            if self.ps[i] > 0:
                cand.append(xs[i] - contmin[i] / (r * ps[i]))
        for Hf, _ in self.oracle.barrier_forms():
            if Hf:
                den = r * sum(ps[i] for i in Hf)
                cand.append(min(xs[i] for i in Hf) - psimin[Hf] / den)
        G1lo = min(cand) - 1
        cand2 = [F(ys[0])]
        for i in range(self.n1):
            cand2.append(xs[i] - r * (xs[i] - G1lo))
        G2lo = min(cand2) - 1
        for i in range(self.n1):
            if r * ps[i] * (xs[i] - G1lo) < contmin[i] or \
               ps[i] * (xs[i] - G2lo) < r * ps[i] * (xs[i] - G1lo):
                raise InvariantError('tail threshold fails at fiber %d'
                                     % i)
        for Hf, _ in self.oracle.barrier_forms():
            if Hf and r * sum(ps[i] * (xs[i] - G1lo) for i in Hf) \
               < psimin[Hf]:
                raise InvariantError('tail threshold fails at H=%s'
                                     % sorted(Hf))
        # the box: per-rectangle vertex enumeration
        B1 = sorted(set([G1lo] + [F(x) for x in xs]))
        B2 = sorted(set([G2lo] + [F(y) for y in ys]))
        forms = self.oracle.barrier_forms()
        for a1, b1 in zip(B1, B1[1:]):
            for a2, b2 in zip(B2, B2[1:]):
                bad = self._check_rect(Q, forms, a1, b1, a2, b2)
                if bad is not None:
                    return False, bad
        return True, None

    def _check_rect(self, Q, forms, a1, b1, a2, b2):
        r, xs, ps, ys = self.r, self.xs, self.ps, self.ys
        corners = [(a1, a2), (a1, b2), (b1, a2), (b1, b2)]
        # affine pieces of V on the rectangle
        ai_f, ci_f = [], []
        for i in range(self.n1):
            ai_f.append((-r * ps[i], F(0), r * ps[i] * xs[i])
                        if b1 <= xs[i] else (F(0), F(0), F(0)))
            wj = [Q[i][j] for j in range(self.n2) if ys[j] >= b2]
            yj = [ys[j] for j in range(self.n2) if ys[j] >= b2]
            ci_f.append((F(0), -sum(wj, F(0)),
                         sum(w * y for w, y in zip(wj, yj))))
        # affine barrier forms of L on the rectangle
        Lf = []
        for Hf, SH in forms:
            phiA = -r * sum(ps[i] for i in range(self.n1)
                            if i not in Hf and xs[i] >= b1)
            phiC = r * sum(ps[i] * xs[i] for i in range(self.n1)
                           if i not in Hf and xs[i] >= b1)
            sw = [(z, m) for z, m in SH if z >= b2]
            psiB = -sum((m for _, m in sw), F(0)) if sw else F(0)
            Lf.append((phiA, psiB,
                       phiC + sum(m * z for z, m in sw)))
        # prune dominated barrier forms (exact four-corner test)
        keep = []
        for f in Lf:
            dominated = False
            for g in Lf:
                if g is f or g == f:
                    continue
                d = [_aff_eval(_aff_sub(g, f), *c) for c in corners]
                if all(v >= 0 for v in d) and any(v > 0 for v in d):
                    dominated = True
                    break
            if not dominated:
                keep.append(f)
        keep = list(set(keep))
        # kink lines: fiber switches + pairs of kept barrier forms
        lines = [_aff_sub(ai_f[i], ci_f[i]) for i in range(self.n1)]
        lines += [_aff_sub(f, g) for f, g in combinations(keep, 2)]
        lines = [ln for ln in lines
                 if not (ln[0] == 0 and ln[1] == 0)]
        pts = set(corners)
        edges = [((a1, a2), (b1, a2)), ((a1, b2), (b1, b2)),
                 ((a1, a2), (a1, b2)), ((b1, a2), (b1, b2))]
        for ln in lines:
            A, B, C = ln
            for (p, q) in edges:
                if p[0] == q[0]:            # vertical edge k1 = p[0]
                    if B != 0:
                        k2 = -(A * p[0] + C) / B
                        if a2 <= k2 <= b2:
                            pts.add((p[0], k2))
                else:                       # horizontal edge k2 = p[1]
                    if A != 0:
                        k1 = -(B * p[1] + C) / A
                        if a1 <= k1 <= b1:
                            pts.add((k1, p[1]))
        for l1, l2 in combinations(lines, 2):
            det = l1[0] * l2[1] - l2[0] * l1[1]
            if det == 0:
                continue
            k1 = (-l1[2] * l2[1] + l2[2] * l1[1]) / det
            k2 = (-l1[0] * l2[2] + l2[0] * l1[2]) / det
            if a1 <= k1 <= b1 and a2 <= k2 <= b2:
                pts.add((k1, k2))
        for (k1, k2) in pts:
            v = self.V(Q, k1, k2)
            L = self.oracle.L(k1, k2)
            if v < L:
                raise InvariantError(
                    'V < L at (%s,%s): weak duality breach' % (k1, k2))
            if v > L:
                return (k1, k2)
        return None

    # --------------------------------------------------- the driver
    def decide(self, max_cuts=12):
        """the cutting-plane iteration.  Returns a dict with 'verdict' and,
        depending on the verdict, the certified coupling 'Q', or the
        conflicting 'pair' / 'core' with its exact 'gap'; 'iterations'
        counts the portfolio LPs solved."""
        C0 = self.seed_canaries()
        base = dict(n1=self.n1, n2=self.n2, r=str(self.r),
                    dimM=self.dimM)
        for it in range(max_cuts):
            Q, val = self.portfolio_lp(C0)
            Lsum = sum(self.oracle.L(k1, k2) for k1, k2 in C0)
            slack = val - Lsum
            if slack < 0:
                raise InvariantError('portfolio below sum of L')
            if slack > 0:
                return self._conflict(C0, slack, dict(base, iterations=it + 1))
            ok, witness = self.certify(Q)
            if ok:
                return dict(base, verdict=UNIVERSAL, iterations=it + 1,
                            Q=[[str(v) for v in row] for row in Q],
                            certificate='V == L on R^2 (exact)')
            C0 = C0 + [witness]
        return dict(base, verdict=UNDECIDED,
                    reason='cutting-plane budget exhausted (%d cuts)'
                    % max_cuts, n_seed=len(C0))

    # --------------------------------------------------- the conflict
    def _conflict(self, C0, slack, base):
        """a minimal conflicting subset of the current Canary set, by
        deletion, and its exact gap."""
        core = list(C0)
        i = 0
        while i < len(core):
            trial = core[:i] + core[i + 1:]
            if not trial:
                break
            if self._infeasible(trial):
                core = trial
            else:
                i += 1
        if not self._infeasible(core):
            raise InvariantError('slack > 0 but no infeasible subset')
        if len(core) == 2:
            gap = self.pair_gap(core[0], core[1])
            if gap <= 0:
                raise InvariantError('core pair with nonpositive gap')
            return dict(base, verdict=NOT_UNIVERSAL_PAIR,
                        pair=[[str(a), str(b)] for a, b in core],
                        gap=str(gap), gap_float=float(gap))
        for ca, cb in combinations(core, 2):
            gap = self.pair_gap(ca, cb)
            if gap > 0:
                return dict(base, verdict=NOT_UNIVERSAL_PAIR,
                            pair=[[str(ca[0]), str(ca[1])],
                                  [str(cb[0]), str(cb[1])]],
                            gap=str(gap), gap_float=float(gap),
                            via='pair inside a core of size %d' % len(core))
        return dict(base, verdict=NOT_UNIVERSAL_SET,
                    core=[[str(a), str(b)] for a, b in core],
                    core_size=len(core))

    def _infeasible(self, canaries):
        """is {Q : V(Q; c) <= L(c) for all c in canaries} empty?  Exact."""
        try:
            Q, val = self.portfolio_lp(canaries)
        except LPError:
            return True
        Lsum = sum(self.oracle.L(k1, k2) for k1, k2 in canaries)
        return val > Lsum
