"""The sharp lower bound L of a single Canary, from the marginals alone, by
two structurally different exact routes that can be cross-checked.

Conventions: the second exercise slope is normalised to s2 = 1 and the
first is s1 = r; the Canary has strikes (k1, k2).  Paper references are to
the SSRN version of November 2026 of "Robust Lower Bound for a Canary Call".

Route S (the decision route), the closed form (5.6):
    L(k1,k2) = max over H subset supp q1 of
               [ r sum_{i not in H} p_i (x_i-k1)^+  +  C_{S_H}(k2) ],
    S_H = the shadow of q1|_H in q2 (Proposition 5.6),
with the 2^n shadows computed once per instance and certified.

Route A (the cross-check), the dual in Section 6: L is the supremum over
affine mu(z) = c + theta z of
    F(mu) = sum_i p_i (g_i - mu(x_i))^+ + sum_j w_j min(mu(y_j), e_j),
    g_i = max(r (x_i-k1)^+, (x_i-k2)^+),   e_j = (y_j-k2)^+,
and the supremum is attained on the finite set of pairwise intersections of
the n+m lines {c + theta x_i = g_i}, {c + theta y_j = e_j}, together with
each line's points at theta = 0 and theta = 1.

A disagreement between the two routes raises OracleMismatch: it would be a
bug in this implementation, not a property of the instance.
"""

from fractions import Fraction as F
from itertools import combinations

from .exact import cpot, measure, shadow, supp


class OracleMismatch(Exception):
    """Route S and Route A disagree -- a bug in the implementation."""


def _subsets(n):
    for mask in range(1 << n):
        yield [i for i in range(n) if mask >> i & 1]


class LowerBoundOracle:
    """per-instance oracle: precomputes the 2^n certified shadows S_H."""

    def __init__(self, q1_pairs, q2_pairs, r):
        if isinstance(r, float):
            raise TypeError('float r rejected')
        self.q1 = measure(q1_pairs)
        self.q2 = measure(q2_pairs)
        self.r = F(r)
        self.n1 = len(self.q1)
        self.xs = [z for z, _ in self.q1]
        self.ps = [m for _, m in self.q1]
        self.ys = supp(self.q2)
        self.ws = [m for _, m in self.q2]
        self.SH = {}
        for H in _subsets(self.n1):
            nu = measure([(self.xs[i], self.ps[i]) for i in H])
            if not H:
                self.SH[frozenset()] = []
            else:
                self.SH[frozenset(H)] = shadow(nu, self.q2)

    # ---------------- Route S: the shadow closed form ----------------
    def closed_form(self, k1, k2):
        k1, k2 = F(k1), F(k2)
        best = None
        for H in _subsets(self.n1):
            Hf = frozenset(H)
            term = self.r * sum(self.ps[i] * (self.xs[i] - k1)
                                for i in range(self.n1)
                                if i not in Hf and self.xs[i] > k1)
            term += cpot(self.SH[Hf], k2) if H else F(0)
            if best is None or term > best:
                best = term
        return best

    def barrier_forms(self):
        """[(H frozenset, S_H)]: the hold set and its shadow, for every H."""
        return [(Hf, sh) for Hf, sh in self.SH.items()]

    # ---------------- Route A: the arrangement route ------------------
    def arrangement(self, k1, k2):
        k1, k2 = F(k1), F(k2)
        gs = [max(self.r * max(x - k1, F(0)), max(x - k2, F(0)))
              for x in self.xs]
        es = [max(y - k2, F(0)) for y in self.ys]
        lines = ([(x, g) for x, g in zip(self.xs, gs)]
                 + [(y, e) for y, e in zip(self.ys, es)])

        def Fval(c, th):
            v = sum(q * max(F(0), g - (c + th * x))
                    for (x, q), g in zip(self.q1, gs))
            v += sum(w * min(c + th * y, e)
                     for (y, w), e in zip(self.q2, es))
            return v

        best = None
        for (a, A), (b, B) in combinations(lines, 2):
            if a == b:
                continue
            th = (A - B) / (a - b)
            c = A - th * a
            v = Fval(c, th)
            if best is None or v > best:
                best = v
        for (a, A) in lines:
            for th in (F(0), F(1)):
                v = Fval(A - th * a, th)
                if best is None or v > best:
                    best = v
        return best

    # ---------------- the oracle ---------------------------------------
    def L(self, k1, k2, cross=False):
        """L(k1, k2) by the closed form; with cross=True the arrangement
        route is evaluated too and any mismatch raises."""
        vS = self.closed_form(k1, k2)
        if cross:
            vA = self.arrangement(k1, k2)
            if vS != vA:
                raise OracleMismatch(
                    'L mismatch at (%s,%s): closed form %s vs arrangement %s'
                    % (k1, k2, vS, vA))
        return vS
