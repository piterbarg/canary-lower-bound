"""Figures 10 and 11 of "Robust Lower Bound for a Canary Call":
  fig_windows_01.png -- universality regions of the three- and four-atom
      dilation families in the (c, r) plane: boundaries c = 4r-3 (three
      atoms, r in [1,3/2)) and c = 3r-2 (four atoms, all r >= 1), c in (1,3)
      (Theorem 12.1).
  fig_column_01.png  -- the binomial-flow thresholds r*(n), n = 2..8, two
      steps, p = 1/2, with the lower bound 1 + 1/(4(n+1)) (Theorem 12.2).
Every number is a theorem constant; nothing is fitted.
Output: figures/output/
"""
import os
from fractions import Fraction as F
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output')
BLUE, ORANGE, GRAY, INK = '#2a78d6', '#eb6834', '#8a8984', '#0b0b0b'
RC = {'font.size': 10, 'axes.spines.top': False,
      'axes.spines.right': False, 'axes.edgecolor': GRAY,
      'xtick.color': GRAY, 'ytick.color': GRAY}


def main():
    os.makedirs(OUT, exist_ok=True)
    matplotlib.rcdefaults()
    plt.rcParams.update(RC)
    # ---------------------------------------------------------------- windows
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    r = np.linspace(1, 2, 401)
    c3 = 4 * r - 3          # three-atom boundary, theorem range r in [1, 3/2)
    c4 = 3 * r - 2          # four-atom boundary, all r >= 1
    m3 = r < 1.5
    ax.fill_betweenx(r[m3], 1, np.minimum(c3[m3], 3), color=BLUE, alpha=0.18, lw=0)
    ax.fill_betweenx(r, 1, np.minimum(c4, 3), color=ORANGE, alpha=0.18, lw=0)
    ax.plot(np.minimum(c3[m3], 3), r[m3], color=BLUE, lw=2, label=r'three atoms: $c=4r-3$')
    ax.plot(np.minimum(c4, 3), r, color=ORANGE, lw=2, label=r'four atoms: $c=3r-2$')
    ax.axvline(3, color=GRAY, lw=0.8, ls=':')
    ax.axvline(1, color=GRAY, lw=0.8, ls=':')
    ax.set_xlim(0.9, 3.1); ax.set_ylim(1, 2)
    ax.set_xlabel('dilation factor $c$'); ax.set_ylabel('slope ratio $r=s_1/s_2$')
    ax.text(1.35, 1.75, 'universal\n(chain coupling)', color=INK, ha='center')
    ax.text(2.55, 1.12, 'not universal', color=INK, ha='center')
    ax.legend(frameon=False, loc='upper left')
    ax.grid(axis='y', color='#e6e5e1', lw=0.6)
    fig.savefig(os.path.join(OUT, 'fig_windows_01.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ---------------------------------------------------------------- column
    ns = list(range(2, 9))
    rstar = [F(13, 12), F(17, 16), F(21, 20), F(25, 24), F(29, 28), F(31, 30), F(1153, 1120)]
    lb = [1 + F(1, 4 * (n + 1)) for n in ns]
    for n, a, b in zip(ns, rstar, lb):
        assert a >= b, (n, a, b)
        if n <= 6:
            assert a == b, (n, a, b)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(ns, [float(v) for v in lb], color=GRAY, lw=1.6, ls='--', marker='o', ms=5,
            mfc='white', label=r'lower bound $1+1/(4(n+1))$')
    ax.plot(ns, [float(v) for v in rstar], color=BLUE, lw=2, marker='o', ms=7,
            label=r'exact threshold $r^*(n)$')
    for n, v in zip(ns, rstar):
        ax.annotate(f'{v.numerator}/{v.denominator}', (n, float(v)), textcoords='offset points',
                    xytext=(0, 8), ha='center', fontsize=8.5, color=INK)
    ax.set_xlabel('$n$ (steps before the first exercise date)')
    ax.set_ylabel('$r^*$')
    ax.set_xticks(ns)
    ax.set_ylim(1.0, 1.11)
    ax.legend(frameon=False, loc='upper right')
    ax.grid(axis='y', color='#e6e5e1', lw=0.6)
    fig.savefig(os.path.join(OUT, 'fig_column_01.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('figures written to', OUT)


if __name__ == '__main__':
    main()
