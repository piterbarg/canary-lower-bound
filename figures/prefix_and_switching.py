"""Figures 8 and 9 of "Robust Lower Bound for a Canary Call", on the marginals
of Example 9.5: q1 = (1/4, 1/2, 1/4) on {-1, 0, 1} and q2 its dilation by 2.
  fig_prefix_01.png    -- the slice k1 -> L(k1, k2) at k2 = 1/2 and equal
      slopes as the upper envelope of the minimal barrier values of the four
      prefixes; the four other hold sets lie on or below it (Theorem 9.4).
  fig_switching_01.png -- the switching values of the three atoms as functions
      of the second strike, at slope ratios 1 and 5/4 (Section 9.5).  They are
      read off the slice: k1 -> L(k1, k2) is convex and piece-wise linear, at
      each kink the largest maximizing hold set grows by a group of atoms, and
      the kink is the switching value of those atoms (Theorems 9.2 and 9.3).
Every plotted number is an exact rational; floats only at plotting.
Output: figures/output/
"""
import itertools
import os
import sys
from fractions import Fraction as F
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from canary.exact import cpot  # noqa: E402
from canary.lower_bound import LowerBoundOracle  # noqa: E402
from canary import threshold as T  # noqa: E402

OUT = os.path.join(HERE, 'output')
INK, GRAY, GRID = '#0b0b0b', '#8a8984', '#e1e0d9'
FIB = ['#2a78d6', '#eb6834', '#1baf7a']
RC = {'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
      'axes.edgecolor': '#c3c2b7', 'axes.labelcolor': INK, 'xtick.color': INK,
      'ytick.color': INK, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
      'legend.frameon': False, 'mathtext.default': 'regular'}

E = T.binom_centred(2)
Q2 = T.dilate(E, 2)


def setname(xs, H):
    return '{' + ', '.join(str(xs[i]) for i in sorted(H)) + '}' if H else 'nothing'


def switching_values(O, k2):
    """the switching values (kappa-hat_i(k2)) of the atoms of q1, exactly: the
    kinks of the slice k1 -> L(k1, k2), each attributed to the atoms by which
    the largest maximizing hold set grows there."""
    xs, ps, r = O.xs, O.ps, O.r
    n = len(xs)
    Hs = list(O.SH.keys())
    cont = {H: (cpot(O.SH[H], k2) if H else F(0)) for H in Hs}
    M = (xs[-1] - xs[0]) + 10
    bounds = [xs[0] - M] + list(xs) + [xs[-1] + M]
    segments = []                       # (a, b, the hold sets attaining the maximum on (a, b))
    for m in range(n + 1):              # between consecutive atoms every barrier value is affine
        a, b = bounds[m], bounds[m + 1]
        aff = {}
        for H in Hs:
            idx = [i for i in range(n) if i not in H and i >= m]
            aff[H] = (-r * sum(ps[i] for i in idx), r * sum(ps[i] * xs[i] for i in idx) + cont[H])
        cands = {a, b}
        for H, G in itertools.combinations(Hs, 2):
            (s1, c1), (s2, c2) = aff[H], aff[G]
            if s1 != s2:
                k = (c2 - c1) / (s1 - s2)
                if a < k < b:
                    cands.add(k)
        cands = sorted(cands)
        for u, v in zip(cands, cands[1:]):
            mid = (u + v) / 2
            vals = {H: s * mid + c for H, (s, c) in aff.items()}
            top = max(vals.values())
            act = frozenset(H for H, val in vals.items() if val == top)
            if segments and segments[-1][2] == act and segments[-1][1] == u:
                segments[-1] = (segments[-1][0], v, act)
            else:
                segments.append((u, v, act))
    kappa = {}
    for (_, b1, A), (_, _, B) in zip(segments, segments[1:]):
        HL, HR = frozenset().union(*A), frozenset().union(*B)      # the largest maximizers
        assert HL in A and HR in B
        if HL < HR:
            for i in HR - HL:
                assert i not in kappa
                kappa[i] = b1
    return [kappa[i] for i in range(n)]


def panel_slice(ax, r, k2):
    O = LowerBoundOracle(E, Q2, r)
    xs, ps = O.xs, O.ps
    n = len(xs)
    grid = [F(-22, 10) + F(j, 100) * F(35, 10) for j in range(101)]

    def barrier(H, k1):
        return O.r * sum(ps[i] * max(xs[i] - k1, F(0)) for i in range(n) if i not in H) + \
            (cpot(O.SH[H], k2) if H else F(0))
    prefixes = [frozenset(range(j)) for j in range(n + 1)]
    for H in O.SH:
        if H in prefixes:
            continue
        ax.plot([float(k) for k in grid], [float(barrier(H, k)) for k in grid], color=GRAY, lw=1.0, ls='--', zorder=1)
    env = [float(max(barrier(H, k) for H in O.SH)) for k in grid]
    ax.plot([float(k) for k in grid], env, color=FIB[0], lw=5, alpha=0.30, solid_capstyle='butt', zorder=2)
    for H in prefixes:
        ys = [float(barrier(H, k)) for k in grid]
        ax.plot([float(k) for k in grid], ys, color=INK, lw=1.1, zorder=3)
        ax.annotate('hold ' + setname(xs, H), (float(grid[0]), ys[0]), xytext=(4, 3), textcoords='offset points',
                    fontsize=7.5, color=INK, zorder=6,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
    ax.set_xlabel('first strike $k_1$')
    ax.set_ylabel('barrier value')
    ax.set_title('lower bound at $k_2$ = %s, $s_1/s_2$ = %s' % (k2, r), loc='left', fontsize=9)
    ax.legend(handles=[Line2D([], [], color=FIB[0], lw=5, alpha=0.30, label='lower bound'),
                       Line2D([], [], color=INK, lw=1.1, label='prefix hold sets'),
                       Line2D([], [], color=GRAY, lw=1.0, ls='--', label='other hold sets')],
              loc='upper right', fontsize=7.5)


def panel_kappa(ax, r, title):
    O = LowerBoundOracle(E, Q2, r)
    n = len(O.xs)
    lo, hi = F(-27, 10), F(27, 10)
    ks = [lo + (hi - lo) * F(2 * j + 1, 482) for j in range(241)]
    kap = [[] for _ in range(n)]
    for k2 in ks:
        vals = switching_values(O, k2)
        for i in range(n):
            kap[i].append(float(vals[i]))
    for y in O.ys:
        ax.axvline(float(y), color=INK, lw=0.6, ls=':', alpha=0.6, zorder=0)
    for i in range(n):
        ax.plot([float(k) for k in ks], kap[i], color=FIB[i], lw=2, zorder=2 + i)
        ax.annotate('x = %s' % O.xs[i], (float(ks[-1]), kap[i][-1]), xytext=(4, -3), textcoords='offset points',
                    fontsize=8, color=INK)
    ax.set_xlabel('second strike $k_2$')
    ax.set_ylabel('switching value')
    ax.set_xlim(float(lo), float(hi) + 0.7)
    ax.set_title(title, loc='left', fontsize=9)


def main():
    os.makedirs(OUT, exist_ok=True)
    matplotlib.rcdefaults()
    with plt.rc_context(RC):
        fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.5), constrained_layout=True)
        panel_slice(ax, F(1), F(1, 2))
        fig.savefig(os.path.join(OUT, 'fig_prefix_01.png'), dpi=200)
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), constrained_layout=True)
        panel_kappa(axes[0], F(1), '(a) switching values, $s_1/s_2$ = 1')
        panel_kappa(axes[1], F(5, 4), '(b) switching values, $s_1/s_2$ = 5/4')
        fig.savefig(os.path.join(OUT, 'fig_switching_01.png'), dpi=200)
        plt.close(fig)
    print('figures written to', OUT)


if __name__ == '__main__':
    main()
