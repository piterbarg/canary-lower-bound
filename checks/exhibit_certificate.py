"""The hand-checkable certificate for the counterexample of Proposition 13.1
("Robust Lower Bound for a Canary Call", SSRN version of November 2026):
two Canaries without a common minimal measure.

q1 = (1/4, 1/2, 1/4) on {-1, 0, 1},  q2 = the same masses on {-2, 0, 2},
s1 = s2 = 1,  a = (-1/3, 0),  b = (0, 1).

The martingale couplings form the two-parameter family
    fiber x=-1: (u + 1/8, 1/8 - 2u, u)
    fiber x= 0: (w, 1/4 + 2u + 2v, w),  w = 1/8 - u - v
    fiber x= 1: (v, 1/8 - 2v, v + 1/8),  0 <= u, v <= 1/16.
Checks: (H1) every member is a martingale coupling and the square is tight;
(H2) the polytope has dimension 2; (H3) the closed forms
V_a = 2u + max(1/6, 1/4 - 2u - 2v) + max(1/3, 2v + 1/4) and V_b = 3/8 - v
match the exact values on a rational grid; (H4) min (V_a + V_b) = 5/6,
attained only at (u, v) = (0, 1/24), while min V_a = 1/2 and min V_b = 5/16,
so the portfolio gap is 1/48.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from canary.decider import UniversalityDecider  # noqa: E402


def g3(s):
    s = F(s)
    return [(-s, F(1, 4)), (F(0), F(1, 2)), (s, F(1, 4))]


def run(verbose=True):
    dec = UniversalityDecider(g3(1), g3(2), F(1))
    ca, cb = (F(-1, 3), F(0)), (F(0), F(1))

    def Qof(u, v):
        w = F(1, 8) - u - v
        return [[u + F(1, 8), F(1, 8) - 2 * u, u],
                [w, F(1, 4) + 2 * u + 2 * v, w],
                [v, F(1, 8) - 2 * v, v + F(1, 8)]]

    def V(Q, k1, k2):
        tot = F(0)
        for i, (x, p) in enumerate(zip(dec.xs, dec.ps)):
            cont = sum(Q[i][j] * max(y - k2, F(0)) for j, y in enumerate(dec.ys))
            tot += max(p * max(x - k1, F(0)), cont)
        return tot

    ok = [True]

    def check(name, cond):
        ok[0] &= bool(cond)
        if verbose:
            print(('PASS ' if cond else 'FAIL ') + name)

    def is_coupling(Q):
        for i, (x, p) in enumerate(zip(dec.xs, dec.ps)):
            if sum(Q[i]) != p or sum(Q[i][j] * y for j, y in enumerate(dec.ys)) != p * x:
                return False
            if any(v < 0 for v in Q[i]):
                return False
        return all(sum(Q[i][j] for i in range(3)) == dec.ws[j] for j in range(3))

    grid = [F(k, 192) for k in range(0, 13)]   # 0 .. 1/16
    check('H1 family in M on the 13x13 grid', all(is_coupling(Qof(u, v)) for u in grid for v in grid))
    check('H1b corners outside the square violate positivity',
          not is_coupling(Qof(F(1, 16) + F(1, 192), F(0))) and not is_coupling(Qof(F(0), F(-1, 192))))
    check('H2 dim M == 2', dec.dimM == 2)
    Va = lambda u, v: 2 * u + max(F(1, 6), F(1, 4) - 2 * u - 2 * v) + max(F(1, 3), 2 * v + F(1, 4))
    Vb = lambda u, v: F(3, 8) - v
    check('H3 V_a closed form == exact V on grid', all(V(Qof(u, v), *ca) == Va(u, v) for u in grid for v in grid))
    check('H3 V_b closed form == exact V on grid', all(V(Qof(u, v), *cb) == Vb(u, v) for u in grid for v in grid))
    Qa, La = dec.portfolio_lp([ca]); Qb, Lb = dec.portfolio_lp([cb]); Qp, P = dec.portfolio_lp([ca, cb])
    check('H4 LP values L_a=1/2, L_b=5/16, L_port=5/6', (La, Lb, P) == (F(1, 2), F(5, 16), F(5, 6)))
    # for u >= 0 the sum is >= f(v) := max(1/6, 1/4-2v) + max(1/3, 2v+1/4) + 3/8 - v,
    # since 2u + max(1/6, 1/4-2u-2v) >= max(1/6, 1/4-2v); f is PL with its kink at v = 1/24, f = 5/6
    f = lambda v: max(F(1, 6), F(1, 4) - 2 * v) + max(F(1, 3), 2 * v + F(1, 4)) + F(3, 8) - v
    check('H4 f(1/24) = 5/6 and f >= 5/6 on [0,1/16]',
          f(F(1, 24)) == F(5, 6) and all(f(F(k, 960)) >= F(5, 6) for k in range(0, 61)))
    check('H4 2u + max(1/6, 1/4-2u-2v) >= max(1/6, 1/4-2v) on grid',
          all(2 * u + max(F(1, 6), F(1, 4) - 2 * u - 2 * v) >= max(F(1, 6), F(1, 4) - 2 * v)
              for u in grid for v in grid))
    check('H4 attained: V(Q(0,1/24)) sums to 5/6',
          V(Qof(F(0), F(1, 24)), *ca) + V(Qof(F(0), F(1, 24)), *cb) == F(5, 6))
    check('H4 the portfolio-LP coupling is Q(0,1/24)', Qp == Qof(F(0), F(1, 24)))
    check('H4 min V_b = 5/16 at v=1/16; min V_a = 1/2 at (0,1/24)',
          Vb(F(0), F(1, 16)) == F(5, 16) and Va(F(0), F(1, 24)) == F(1, 2))
    return ok[0]


if __name__ == '__main__':
    print('ALL PASS' if run() else 'SOME FAIL')
