"""Figure 10 of "Robust Lower Bound for a Canary Call": Step 4 of the proof
of Proposition 6.2 (Appendix C.1) -- replacing a maximizer rho_0 by the node
interpolant of rho_1 = max(rho_0, e_2) with controlled half-line slopes.
  fig_step4_01.png
Every plotted number is an exact rational; matplotlib only renders.
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
      'xtick.color': GRAY, 'ytick.color': GRAY,
      'axes.labelcolor': INK}


def main():
    os.makedirs(OUT, exist_ok=True)
    matplotlib.rcdefaults()
    plt.rcParams.update(RC)
    # instance: nodes t_1 < ... < t_5 (atoms of q1 and q2), the call e_2 with
    # slope s2 = 1 and strike k2 = 2/5, and a convex piece-wise linear rho_0
    # with kinks at the nodes 0 and 1 that dips below e_2
    nodes = [F(-2), F(-1), F(0), F(1), F(2)]
    s2, k2 = F(1), F(2, 5)

    def e2(y):
        return s2 * max(y - k2, F(0))

    def rho0(y):
        return F(-3, 25) - F(1, 5) * y + F(9, 20) * max(y, F(0)) + F(3, 4) * max(y - F(1), F(0))

    def rho1(y):
        return max(rho0(y), e2(y))

    # rho_0 crosses the flat part of e_2 at -3/5 and stays below the call part
    # (slope 1, parallel) for y > 1; the off-node kinks of rho_1 are -3/5 and k_2
    cross0, cross = F(-3, 5), k2
    assert rho0(cross0) == e2(cross0) and all(rho0(y) < e2(y) for y in [F(1), F(2), F(3)])
    vals = [rho1(t) for t in nodes]
    sig = [(vals[i + 1] - vals[i]) / (nodes[i + 1] - nodes[i]) for i in range(len(nodes) - 1)]
    sl_left, sl_right = min(F(0), sig[0]), max(s2, sig[-1])

    def rhobar(y):
        if y <= nodes[0]:
            return vals[0] + sl_left * (y - nodes[0])
        if y >= nodes[-1]:
            return vals[-1] + sl_right * (y - nodes[-1])
        for i in range(len(nodes) - 1):
            if nodes[i] <= y <= nodes[i + 1]:
                lam = (y - nodes[i]) / (nodes[i + 1] - nodes[i])
                return vals[i] + lam * (vals[i + 1] - vals[i])

    grid = [F(-5, 2) + F(i, 40) for i in range(0, 201)]
    for y in grid:
        assert rhobar(y) >= rho1(y) >= e2(y), y
    for t in nodes:
        assert rhobar(t) == rho1(t)

    xs = np.array([float(y) for y in grid])
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.plot(xs, [float(e2(y)) for y in grid], color=GRAY, lw=1.4, label=r'$e_2$')
    ax.plot(xs, [float(rho0(y)) for y in grid], color=BLUE, lw=1.6, ls='--', label=r'$\rho_0$ (maximizer)')
    ax.plot(xs, [float(rho1(y)) for y in grid], color=BLUE, lw=1.0, alpha=0.55, label=r'$\rho_1=\max(\rho_0,e_2)$')
    ax.plot(xs, [float(rhobar(y)) for y in grid], color=ORANGE, lw=2.0, label=r'$\bar\rho$ (node interpolant)')
    for i, t in enumerate(nodes):
        ax.plot([float(t)], [float(vals[i])], 'o', ms=5, color=ORANGE, mec='white', mew=1.0, zorder=5)
        ax.annotate(r'$t_%d$' % (i + 1), (float(t), -0.22), ha='center', va='top', fontsize=9, color=INK,
                    annotation_clip=False)
    ax.set_xticks([float(t) for t in nodes]); ax.set_xticklabels([])
    for y in [cross0, k2]:
        ax.plot([float(y)], [float(rho1(y))], 'o', ms=9, mfc='none', mec=INK, mew=1.1, zorder=6)
    ax.annotate(r'kinks of $\rho_1$ off the nodes', (float(cross) + 0.08, float(rho1(cross)) - 0.05),
                ha='left', va='top', fontsize=8.5, color=INK)
    ax.axvline(float(k2), color=GRAY, lw=0.6, ls=':')
    ax.annotate(r'$k_2$', (float(k2), 1.95), ha='center', va='bottom', fontsize=9, color=GRAY)
    ax.set_xlim(-2.5, 2.5); ax.set_ylim(-0.15, 2.1)
    ax.set_yticks([0, 0.5, 1, 1.5, 2])
    ax.grid(axis='y', color='#e5e5e5', lw=0.6)
    ax.legend(loc='upper left', frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_step4_01.png'), dpi=200)
    print('figure written to', OUT, '; crossings', cross0, cross, '; half-line slopes', sl_left, sl_right)


if __name__ == '__main__':
    main()
