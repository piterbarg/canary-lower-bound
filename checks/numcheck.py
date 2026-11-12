"""Recompute every number in "Robust Lower Bound for a Canary Call" (SSRN
version of November 2026) and its companion paper that is the result of a
computation, from the marginals alone, in exact rational arithmetic.

Two independent routes are used where both apply: the universality decider
(canary.decider) and the marginals-only threshold calculator
(canary.threshold), which share no code.  Quantifiers are asserted as the
papers state them: at a window boundary, "universal at the boundary and not
universal just beyond it".

Blocks:
  X  the counterexample of Proposition 12.1 (checks/exhibit_certificate.py)
     and the decider on the dilation family (checks/exhibit_rows.py)
  H  the three Canaries that are pairwise universal but not jointly, the
     instance after the Helly reduction theorem (checks/helly_instance.py)
  R  every threshold constant (canary.threshold.selftest)
  W  the dilation windows of Theorem 11.1 as if-and-only-if statements
  C  the binomial-flow thresholds of Theorem 11.2 and the (3,3), (4,4) cells
  N  Proposition 10.1: the sufficient condition is not necessary
Run time: one to two minutes (exact LPs).
"""
import os
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from canary.decider import UniversalityDecider, UNIVERSAL, NOT_UNIVERSAL_PAIR  # noqa: E402
from canary import threshold as T  # noqa: E402
import exhibit_certificate  # noqa: E402
import exhibit_rows  # noqa: E402
import helly_instance  # noqa: E402


def run(verbose=True):
    ok = [True]

    def check(name, cond, detail=''):
        ok[0] &= bool(cond)
        if verbose:
            print(('PASS ' if cond else 'FAIL ') + name, detail)

    def verdict(q1, q2, r):
        return UniversalityDecider(T.measure(q1), T.measure(q2), F(r)).decide()['verdict']

    t0 = time.time()
    # ---- X: the counterexample ------------------------------------------
    check('X exhibit hand certificate', exhibit_certificate.run(verbose=False))
    check('X decider on the dilation family (five quoted points)', exhibit_rows.run(verbose=False))

    # ---- H: the Helly instance ------------------------------------------
    check('H three Canaries, pairwise universal but not jointly (gap 1/108)',
          helly_instance.run(verbose=False))

    # ---- R: thresholds via the marginals-only calculator -----------------
    check('R threshold constants', T.selftest(verbose=False))

    # ---- W: the dilation windows as iffs (Theorem 11.1) -----------------
    G3 = T.binom_centred(2); G4 = T.binom_centred(3)
    eps = F(1, 50)
    for r in (F(1), F(9, 8), F(5, 4), F(7, 5)):          # three atoms: r in [1, 3/2)
        cb = 4 * r - 3
        if 1 < cb < 3:
            check('W three-atom r=%s: universal at c=c*=%s' % (r, cb),
                  verdict(G3, T.dilate(G3, cb), r) == UNIVERSAL)
            check('W three-atom r=%s: not universal at c=c*+eps' % r,
                  verdict(G3, T.dilate(G3, cb + eps), r) == NOT_UNIVERSAL_PAIR)
        if cb <= 1:
            check('W three-atom r=%s: not universal at c=1+eps (empty window)' % r,
                  verdict(G3, T.dilate(G3, 1 + eps), r) == NOT_UNIVERSAL_PAIR)
    for r in (F(1), F(4, 3), F(3, 2), F(5, 3), F(2)):     # four atoms: all r >= 1
        cb = 3 * r - 2
        if 1 < cb < 3:
            check('W four-atom r=%s: universal at c=c*=%s' % (r, cb),
                  verdict(G4, T.dilate(G4, cb), r) == UNIVERSAL)
            check('W four-atom r=%s: not universal at c=c*+eps' % r,
                  verdict(G4, T.dilate(G4, cb + eps), r) == NOT_UNIVERSAL_PAIR)
        if cb <= 1:
            check('W four-atom r=%s: not universal at c=1+eps' % r,
                  verdict(G4, T.dilate(G4, 1 + eps), r) == NOT_UNIVERSAL_PAIR)
        if cb >= 3:
            check('W four-atom r=%s: universal at c=3-eps' % r,
                  verdict(G4, T.dilate(G4, 3 - eps), r) == UNIVERSAL)
    U3 = [(-1, F(1, 3)), (0, F(1, 3)), (1, F(1, 3))]
    check('W uniform3 x dilation 3: universal at r=1 (r*=1)', verdict(U3, T.dilate(U3, 3), 1) == UNIVERSAL)

    # ---- C: the binomial flow (Theorem 11.2) and the diagonal cells ------
    col = {2: F(13, 12), 3: F(17, 16), 4: F(21, 20), 5: F(25, 24)}
    for n, rs in col.items():
        q1, q2 = T.binom_flow(n, 2)
        check('C flow n=%d: universal at r=r*=%s' % (n, rs), verdict(q1, q2, rs) == UNIVERSAL)
        check('C flow n=%d: not universal at r=1' % n, verdict(q1, q2, 1) == NOT_UNIVERSAL_PAIR)
        check('C flow n=%d: not universal at r=r*-1/1000' % n,
              verdict(q1, q2, rs - F(1, 1000)) == NOT_UNIVERSAL_PAIR)
    q1, q2 = T.binom_flow(3, 3)
    check('C (3,3): universal at r=5/4', verdict(q1, q2, F(5, 4)) == UNIVERSAL)
    check('C (3,3): not universal at r=1', verdict(q1, q2, 1) == NOT_UNIVERSAL_PAIR)
    check('C (3,3): not universal at r=5/4-1/1000', verdict(q1, q2, F(5, 4) - F(1, 1000)) == NOT_UNIVERSAL_PAIR)
    q1, q2 = T.binom_flow(4, 4)
    check('C (4,4): universal at r=37/32', verdict(q1, q2, F(37, 32)) == UNIVERSAL)
    check('C (4,4): not universal at r=1', verdict(q1, q2, 1) == NOT_UNIVERSAL_PAIR)

    # ---- N: the sufficient condition is not necessary (Proposition 10.1) -
    qa1, qa2 = T.CONVERSE_A
    check('N (a): universal at r=1', verdict(qa1, qa2, 1) == UNIVERSAL)
    for r in (F(11, 10), F(5, 4), F(3, 2), F(7, 4), F(9, 5)):
        check('N (a): not universal at r=%s' % r, verdict(qa1, qa2, r) == NOT_UNIVERSAL_PAIR)
    check('N (a): universal at r=r*=137/72', verdict(qa1, qa2, F(137, 72)) == UNIVERSAL)
    qb1, qb2 = T.CONVERSE_B
    check('N (b): universal at r=1', verdict(qb1, qb2, 1) == UNIVERSAL)
    for j in (1, 3, 6, 9):
        check('N (b): not universal at r=1+10^-%d' % j, verdict(qb1, qb2, 1 + F(1, 10 ** j)) == NOT_UNIVERSAL_PAIR)
    check('N (b): not universal at r=49/36-10^-9', verdict(qb1, qb2, F(49, 36) - F(1, 10 ** 9)) == NOT_UNIVERSAL_PAIR)
    check('N (b): universal at r=49/36', verdict(qb1, qb2, F(49, 36)) == UNIVERSAL)
    qe1, kern = T.EXCEPTIONAL_Q1, T.EXCEPTIONAL_KERNEL
    check('N exceptional convolution: universal at r=1 with r*=16/15',
          verdict(qe1, T.convolve(qe1, kern), 1) == UNIVERSAL)
    if verbose:
        print('elapsed %.0fs' % (time.time() - t0))
    return ok[0]


if __name__ == '__main__':
    print('NUMCHECK ALL PASS' if run() else 'NUMCHECK: SOME FAIL')
