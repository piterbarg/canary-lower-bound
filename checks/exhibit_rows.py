"""The cutting-plane decider on the three-atom dilation family at the five
(r, c) points quoted in the Numerics sections of the two papers:

  r = 1,     c = 2:    conflicting pair (-1/3, 0), (0, 1), gap 1/48, two iterations
  r = 1,     c = 3/2:  conflicting pair (-1/6, 0), (0, 1/2), gap 1/72
  r = 1,     c = 5/2:  conflicting pair (-1/2, 0), (0, 3/2), gap 1/40
  r = 21/20, c = 11/10: universal (c <= 4r - 3 = 6/5)
  r = 21/20, c = 2:    not universal (a conflicting pair)

q1 = (1/4, 1/2, 1/4) on {-1, 0, 1}, q2 = its dilation by c.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from canary.decider import UniversalityDecider, UNIVERSAL, NOT_UNIVERSAL_PAIR  # noqa: E402
from canary.threshold import binom_centred, dilate  # noqa: E402

EXPECTED = [
    (F(1), F(2), NOT_UNIVERSAL_PAIR, [['-1/3', '0'], ['0', '1']], '1/48', 2),
    (F(1), F(3, 2), NOT_UNIVERSAL_PAIR, [['-1/6', '0'], ['0', '1/2']], '1/72', None),
    (F(1), F(5, 2), NOT_UNIVERSAL_PAIR, [['-1/2', '0'], ['0', '3/2']], '1/40', None),
    (F(21, 20), F(11, 10), UNIVERSAL, None, None, None),
    (F(21, 20), F(2), NOT_UNIVERSAL_PAIR, None, None, None),
]


def run(verbose=True):
    G3 = binom_centred(2)
    ok = True
    for r, c, verdict, pair, gap, iters in EXPECTED:
        d = UniversalityDecider(G3, dilate(G3, c), r).decide()
        good = d['verdict'] == verdict
        if pair is not None:
            good &= d.get('pair') == pair
        if gap is not None:
            good &= d.get('gap') == gap
        if iters is not None:
            good &= d.get('iterations') == iters
        ok &= good
        if verbose:
            print(('PASS ' if good else 'FAIL ') + 'r=%s c=%s: %s' % (r, c, d['verdict']),
                  'pair', d.get('pair'), 'gap', d.get('gap'), 'iterations', d.get('iterations'))
    return ok


if __name__ == '__main__':
    print('ALL PASS' if run() else 'SOME FAIL')
