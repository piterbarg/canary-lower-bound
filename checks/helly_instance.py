"""Three Canaries that are pairwise universal but not jointly universal: the
instance printed after the Helly reduction theorem of "Robust Lower Bound for
a Canary Call" (appendix "Three Canaries, Pairwise Universal but Not Jointly").
The paper prints the data and the result; this script is the derivation.

q1 = (1/3, 2/3) on {0, 1},  q2 = (1/12, 5/12, 1/4, 1/4) on {-7/2, -13/4, 11/4, 13/2},
s1 = s2 = 1, and the Canaries
    k^a = (-157/48, -13/4),   k^b = (-139/48, -5/2),   k^c = (-7/192, 59/16).

The mechanism.  Let u(Q) be the continuation value of the fiber at x1 = 0.  The
martingale constraints give the continuation value of the other fiber as
C_{q2}(k2) - u(Q), so for k1 < 0
    V(Q; k) = max(-p1 k1, u(Q)) + max(p2 (1 - k1), C_{q2}(k2) - u(Q)).
With k1 = p2 - C_{q2}(k2) the two exercise values add up to C_{q2}(k2), and with
u* = -p1 k1 this becomes V(Q; k) = C_{q2}(k2) + |u(Q) - u*|.  So L(k) = C_{q2}(k2)
as soon as some coupling has u(Q) = u*, and the minimal measures of k are the
couplings on the hyperplane u(Q) = u*.

The couplings.  Writing (a, b) for the masses of the fiber at x1 = 0 on the two
upper atoms of q2, the fiber is (1/3 - a - b - c, c, a, b) with c = 14/3 - 25a - 40b
(zero mean) and the other fiber is q2 minus it: the couplings form a polygon in
the (a, b)-plane (dimension 2).  The minimal measures of the three Canaries are
the segments
    6 a + 39/4 b = 157/144,    21/4 a + 9 b = 139/144,    b = 7/1620
of the polygon.  They cross pairwise at three distinct points, so each pair has
a common minimal measure and the three have none.

Checks: (H1) the parametrization gives martingale couplings and the three
crossing points lie in the polygon; (H2) k1 = p2 - C_{q2}(k2) and the lower
bounds are 63/16, 57/16, 45/64; (H3) V(Q; k) = C_{q2}(k2) + |u(Q) - u*| on a
rational grid of the polygon; (H4) the coupling at each crossing point is
minimal for exactly two of the three Canaries; (H5) the exact portfolio linear
programs: every pair has gap 0, and the three together have the minimal
portfolio value 14191/1728 = 63/16 + 57/16 + 45/64 + 1/108.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from canary.decider import UniversalityDecider  # noqa: E402

Q1 = [(F(0), F(1, 3)), (F(1), F(2, 3))]
Q2 = [(F(-7, 2), F(1, 12)), (F(-13, 4), F(5, 12)), (F(11, 4), F(1, 4)), (F(13, 2), F(1, 4))]
KA, KB, KC = (F(-157, 48), F(-13, 4)), (F(-139, 48), F(-5, 2)), (F(-7, 192), F(59, 16))
CROSS = {'ab': (F(77, 540), F(13, 540)), 'ac': (F(283, 1620), F(7, 1620)), 'bc': (F(667, 3780), F(7, 1620))}


def run(verbose=True):
    dec = UniversalityDecider(Q1, Q2, F(1))
    xs, ps, ys = dec.xs, dec.ps, dec.ys
    w2 = [w for _, w in Q2]
    ok = [True]

    def check(name, cond):
        ok[0] &= bool(cond)
        if verbose:
            print(('PASS ' if cond else 'FAIL ') + name)

    def Qof(a, b):
        c = F(14, 3) - 25 * a - 40 * b
        f1 = [F(1, 3) - a - b - c, c, a, b]
        return [f1, [w - v for w, v in zip(w2, f1)]]

    def inside(Q):
        return all(v >= 0 for row in Q for v in row)

    def is_coupling(Q):
        return all(sum(Q[i]) == ps[i] and sum(Q[i][j] * ys[j] for j in range(4)) == ps[i] * xs[i]
                   for i in range(2)) and all(Q[0][j] + Q[1][j] == w2[j] for j in range(4))

    def V(Q, k):
        k1, k2 = k
        return sum(max(ps[i] * max(xs[i] - k1, F(0)),
                       sum(Q[i][j] * max(ys[j] - k2, F(0)) for j in range(4))) for i in range(2))

    def cq2(k2):
        return sum(w * max(y - k2, F(0)) for y, w in Q2)

    def u(Q, k2):
        return sum(Q[0][j] * max(ys[j] - k2, F(0)) for j in range(4))

    cans = {'a': KA, 'b': KB, 'c': KC}
    L = {t: dec.oracle.L(*k) for t, k in cans.items()}

    # H1
    pts = {t: Qof(*ab) for t, ab in CROSS.items()}
    check('H1 the parametrization gives martingale couplings; the three crossing points are distinct '
          'points of the polygon',
          all(is_coupling(Q) and inside(Q) for Q in pts.values()) and len(set(CROSS.values())) == 3)
    # H2
    check('H2 k1 = p2 - C_q2(k2) for the three Canaries; L = 63/16, 57/16, 45/64',
          all(k[0] == ps[1] - cq2(k[1]) for k in cans.values())
          and [L['a'], L['b'], L['c']] == [F(63, 16), F(57, 16), F(45, 64)]
          and all(L[t] == cq2(cans[t][1]) for t in cans))
    # H3
    grid = [Qof(F(i, 400), F(j, 4000)) for i in range(0, 120) for j in range(0, 400, 7)]
    grid = [Q for Q in grid if inside(Q)]
    check('H3 V(Q; k) = C_q2(k2) + |u(Q) - u*| on %d couplings of the polygon' % len(grid),
          len(grid) > 50 and all(V(Q, k) == cq2(k[1]) + abs(u(Q, k[1]) + ps[0] * k[0])
                                 for Q in grid for k in cans.values()))
    # the three segments in the coordinates (a, b)
    seg = {'a': lambda a, b: 6 * a + F(39, 4) * b - F(157, 144),
           'b': lambda a, b: F(21, 4) * a + 9 * b - F(139, 144),
           'c': lambda a, b: b - F(7, 1620)}
    check('H3 the hyperplanes u(Q) = u* are the three printed segments',
          all((u(Qof(a, b), cans[t][1]) + ps[0] * cans[t][0] == 0) == (seg[t](a, b) == 0)
              for t in cans for a, b in list(CROSS.values()) + [(F(1, 10), F(1, 100))])
          and all(seg[t](*CROSS[pair]) == 0 for pair in CROSS for t in pair))
    # H4
    good = True
    for pair, Q in pts.items():
        other = (set('abc') - set(pair)).pop()
        good &= all(V(Q, cans[t]) == L[t] for t in pair) and V(Q, cans[other]) > L[other]
    check('H4 the coupling at each crossing point is minimal for exactly two of the three Canaries', good)
    check('H4 the coupling at the crossing of a and c misses b by exactly 1/108',
          V(pts['ac'], KB) - L['b'] == F(1, 108))
    # H5
    gaps = [dec.portfolio_lp([cans[s], cans[t]])[1] - L[s] - L[t] for s, t in ('ab', 'ac', 'bc')]
    check('H5 portfolio LP: every pair has gap 0', gaps == [0, 0, 0])
    val = dec.portfolio_lp([KA, KB, KC])[1]
    check('H5 portfolio LP: L_port = 14191/1728, the sum of the three lower bounds plus 1/108',
          val == F(14191, 1728) and val - sum(L.values()) == F(1, 108))
    return ok[0]


if __name__ == '__main__':
    good = run(verbose=True)
    print('HELLY INSTANCE: ALL PASS' if good else 'HELLY INSTANCE: SOME FAIL')
    sys.exit(0 if good else 1)
