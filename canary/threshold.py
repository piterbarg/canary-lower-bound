"""The greedy fibers and the marginals-only threshold r* of Theorem 9.4 of
"Robust Lower Bound for a Canary Call" (SSRN version of November 2026), in
exact rationals, plus the explicit families of marginals used in the paper
and its companion.

W_j = the shadow of p_j delta_{x_j} in the remainder R_{j-1} = q2 - S_{j-1},
    which is the quantile window of R_{j-1} of mass p_j and barycentre x_j
    (Lemma 9.1(a)); it is found by sliding a window of mass p_j along the
    quantile axis of the remainder (the barycentre is continuous and
    non-decreasing in the position) and solving exactly.
C_j(k) = int (y-k)^+ W_j(dy) / p_j;
r*     = max_j sup_k (C_{j+1} - C_j)(k) / (x_{j+1} - x_j),
    the sup being attained at an atom of q2 or on the left ray (Theorem
    9.4), so it is evaluated at the atoms of q2 and one point left of them.

selftest() recomputes the threshold constants printed in the two papers:
  three-atom symmetric law x dilation by c:  r* = (c+3)/4
  four-atom symmetric law x dilation by c:   r* = (c+2)/3
  uniform on {-1,0,1} x dilation by c:        r* = max(1, (c-1)/2)
  Proposition 10.1 (a): 137/72;  (b): 49/36
  the exceptional convolution instance of Section 10: 16/15
  binomial flow, two steps, p = 1/2, n = 2..8:
      13/12, 17/16, 21/20, 25/24, 29/28, 31/30, 1153/1120
  binomial flow (3,3): 5/4;  (4,4): 37/32
  uniform lattice with kernel (2/3 at -1, 1/3 at 2):
      1 + sigma^2 (n-2)/(2n), sigma^2 = 2, for n large enough
"""
from fractions import Fraction as F
from math import comb


def measure(pairs):
    acc = {}
    for z, m in pairs:
        if m:
            acc[F(z)] = acc.get(F(z), F(0)) + F(m)
    return sorted((z, m) for z, m in acc.items() if m > 0)


def call(meas, k):
    return sum(m * max(y - k, F(0)) for y, m in meas)


def window(R, mass, mean):
    """the quantile window of the atomic measure R with the given mass and
    barycentre; returns the sub-measure as a sorted list.  Exact."""
    ys = [y for y, _ in R]; ms = [m for _, m in R]
    tot = sum(ms)
    assert mass <= tot, (mass, tot)

    def win(t):
        """sub-measure of R between quantile levels t and t+mass."""
        out = []; acc = F(0)
        for y, m in R:
            lo, hi = acc, acc + m
            a, b = max(lo, t), min(hi, t + mass)
            if b > a:
                out.append((y, b - a))
            acc = hi
        return out

    def bary(t):
        w = win(t)
        return sum(y * m for y, m in w) / mass

    # breakpoints of t -> bary(t): where t or t+mass crosses a cumulative level
    cum = [F(0)]
    for m in ms:
        cum.append(cum[-1] + m)
    cands = set()
    for c in cum:
        for t in (c, c - mass):
            if F(0) <= t <= tot - mass:
                cands.add(t)
    ts = sorted(cands)
    vals = [bary(t) for t in ts]
    if len(ts) == 1:                      # the window is the whole remainder
        w = win(ts[0])
        assert sum(y * m for y, m in w) == mass * mean
        return w
    assert vals[0] <= mean <= vals[-1], (vals[0], mean, vals[-1])
    for (t0, b0), (t1, b1) in zip(zip(ts, vals), zip(ts[1:], vals[1:])):
        if b0 <= mean <= b1:
            if b1 == b0:
                t = t0
            else:
                t = t0 + (t1 - t0) * (mean - b0) / (b1 - b0)   # bary is affine between breakpoints
            w = win(t)
            assert sum(m for _, m in w) == mass
            assert sum(y * m for y, m in w) == mass * mean
            return w
    raise ValueError('no window')


def greedy_fibers(q1, q2):
    """the greedy fibers W_1, ..., W_n (Section 9.1)."""
    R = list(q2); W = []
    for x, p in q1:
        w = window(R, p, x)
        W.append(w)
        Rd = {y: m for y, m in R}
        for y, m in w:
            Rd[y] -= m
        R = sorted((y, m) for y, m in Rd.items() if m > 0)
    assert not R or sum(m for _, m in R) == 0
    return W


def rstar(q1, q2):
    """the threshold r*(q1, q2) of Theorem 9.4."""
    q1 = measure(q1); q2 = measure(q2)
    n = len(q1)
    if n < 2:
        return F(1)
    W = greedy_fibers(q1, q2)
    ks = [y for y, _ in q2] + [q2[0][0] - 1]
    best = F(0)
    for j in range(n - 1):
        pj, pj1 = q1[j][1], q1[j + 1][1]
        gj = q1[j + 1][0] - q1[j][0]
        for k in ks:
            v = (call(W[j + 1], k) / pj1 - call(W[j], k) / pj) / gj
            best = max(best, v)
    return best


# ---------------------------------------------------------------- families
def binom_centred(N, step=F(1), shift=None):
    """the centred symmetric binomial law on N+1 atoms (N=2: the three-atom
    law 1/4, 1/2, 1/4 on {-1,0,1}; N=3: the four-atom law on half-integers)."""
    sh = F(N, 2) if shift is None else shift
    return measure([(step * (k - sh), F(comb(N, k), 2 ** N)) for k in range(N + 1)])


def dilate(q, c):
    """the dilation of q by the factor c about its mean."""
    q = measure(q); m = sum(x * p for x, p in q)
    return measure([(m + F(c) * (x - m), p) for x, p in q])


def binom_flow(n, kappa, p=F(1, 2)):
    """the binomial flow: q1 = Binomial(n, p), q2 = Binomial(n + kappa, p)
    shifted to the same mean."""
    def binom(N):
        return measure([(F(k), F(comb(N, k)) * p ** k * (1 - p) ** (N - k)) for k in range(N + 1)])
    q1 = binom(n)
    q2 = measure([(y - kappa * p, m) for y, m in binom(n + kappa)])
    return q1, q2


def convolve(q1, kern):
    """q2 = q1 * kern for a mean-zero atomic kernel."""
    return measure([(x + z, p * w) for x, p in q1 for z, w in kern])


# --------------------------------------------------------------- instances
# Proposition 10.1 (a) and (b), and the exceptional convolution instance
CONVERSE_A = ([(-3, F(1, 2)), (2, F(3, 7)), (3, F(1, 14))],
              [(-3, F(1, 2)), (F(-3, 2), F(3, 56)), (1, F(8, 35)),
               (F(3, 2), F(3, 70)), (F(9, 2), F(9, 56)), (F(19, 2), F(1, 70))])
CONVERSE_B = ([(3, F(2, 5)), (4, F(3, 5))],
              [(z, F(1, 10)) for z in (0, 1, 2, 3, 5, 6, 7)] + [(4, F(3, 10))])
EXCEPTIONAL_Q1 = [(0, F(1, 2)), (F(1, 2), F(1, 4)), (1, F(1, 4))]
EXCEPTIONAL_KERNEL = [(F(-5, 12), F(2, 3)), (F(10, 12), F(1, 3))]
FLOW_COLUMN = [F(13, 12), F(17, 16), F(21, 20), F(25, 24), F(29, 28), F(31, 30), F(1153, 1120)]


def selftest(verbose=True):
    """recompute every threshold constant printed in the papers; returns True
    iff all agree."""
    ok = [True]

    def check(name, got, exp):
        good = (got == exp); ok[0] &= good
        if verbose:
            print(('PASS ' if good else 'FAIL ') + name, got, 'expected', exp)

    G3 = binom_centred(2); G4 = binom_centred(3)
    for c in (F(5, 4), F(3, 2), F(2), F(5, 2)):
        check('three-atom r*(c=%s) = (c+3)/4' % c, rstar(G3, dilate(G3, c)), (c + 3) / 4)
        check('four-atom r*(c=%s) = (c+2)/3' % c, rstar(G4, dilate(G4, c)), (c + 2) / 3)
    U3 = [(-1, F(1, 3)), (0, F(1, 3)), (1, F(1, 3))]
    for c in (F(2), F(3), F(4), F(5)):
        check('uniform3 r*(c=%s) = max(1,(c-1)/2)' % c, rstar(U3, dilate(U3, c)), max(F(1), (c - 1) / 2))
    check('Proposition 10.1(a) r* = 137/72', rstar(*CONVERSE_A), F(137, 72))
    check('Proposition 10.1(b) r* = 49/36', rstar(*CONVERSE_B), F(49, 36))
    check('exceptional convolution r* = 16/15',
          rstar(EXCEPTIONAL_Q1, convolve(EXCEPTIONAL_Q1, EXCEPTIONAL_KERNEL)), F(16, 15))
    for n, exp in zip(range(2, 9), FLOW_COLUMN):
        q1, q2 = binom_flow(n, 2)
        check('binomial flow n=%d, two steps, r*' % n, rstar(q1, q2), exp)
    for n, exp in ((3, F(5, 4)), (4, F(37, 32))):
        q1, q2 = binom_flow(n, n)
        check('binomial flow (%d,%d) r*' % (n, n), rstar(q1, q2), exp)
    # uniform lattice, kernel (2/3 at -1, 1/3 at 2): sigma^2 = 2/3 + 4/3 = 2
    kl = [(-1, F(2, 3)), (2, F(1, 3))]
    for n in (10, 12, 16):
        q1 = [(k, F(1, n)) for k in range(n)]
        got = rstar(q1, convolve(q1, kl)); exp = 1 + F(2) * (n - 2) / (2 * n)
        if verbose:
            print(('PASS ' if got == exp else 'INFO ') + 'uniform lattice n=%d r*' % n, got,
                  'closed form', exp, '(the closed form holds for n large enough)')
    return ok[0]


if __name__ == '__main__':
    print('ALL PASS' if selftest() else 'SOME FAIL')
