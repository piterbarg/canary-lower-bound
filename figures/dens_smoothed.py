"""Figure 12 of "Robust Lower Bound for a Canary Call" (Theorem 13.2)
(the counterexample for densities): the smoothed marginals
  q1^{eps,delta} = (1-delta) q1 * law(eps U) + delta m,
  q2^{eps,delta} = (1-delta) q2 * law(eps (U+U')) + delta m,
U triangular on [-1,1], m triangular on [-3,3], drawn at eps = 0.3, delta = 0.05,
with the atoms of q1 (at -1, 0, 1) and q2 (at -2, 0, 2) as faint bars.
Every curve is an explicit closed-form density; nothing is computed.
  fig_dens_smoothed_01.png
Output: figures/output/
"""
import os, sys
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

Q1 = [(-1, .25), (0, .5), (1, .25)]
Q2 = [(-2, .25), (0, .5), (2, .25)]
EPS, DELTA = 0.3, 0.05


def tri(z, h):
    """the triangular density on [-h, h]"""
    return np.maximum(1 - np.abs(z) / h, 0) / h


def tri2(z, h):
    """the density of h (U + U'), U, U' independent triangular on [-1, 1]"""
    t = np.abs(z) / h
    out = np.zeros_like(t)
    m1 = t <= 1
    m2 = (t > 1) & (t <= 2)
    out[m1] = (4 - 6 * t[m1] ** 2 + 3 * t[m1] ** 3) / 6
    out[m2] = (2 - t[m2]) ** 3 / 6
    return out / h


def dens_q1(y, eps=EPS, delta=DELTA):
    return (1 - delta) * sum(p * tri(y - x, eps) for x, p in Q1) + delta * tri(y, 3)


def dens_q2(y, eps=EPS, delta=DELTA):
    return (1 - delta) * sum(p * tri2(y - x, eps) for x, p in Q2) + delta * tri(y, 3)


def main(out=OUT):
    os.makedirs(out, exist_ok=True)
    matplotlib.rcdefaults()
    plt.rcParams.update(RC)
    y = np.linspace(-3.2, 3.2, 2001)
    # sanity: unit mass and zero mean of both densities
    for f in (dens_q1, dens_q2):
        v = f(y)
        dy = y[1] - y[0]
        assert abs(np.sum(v) * dy - 1) < 1e-3 and abs(np.sum(v * y) * dy) < 1e-9
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for x, p in Q1:
        ax.vlines(x, 0, p, color=BLUE, lw=6, alpha=.22)
    for x, p in Q2:
        ax.vlines(x, 0, p, color=ORANGE, lw=6, alpha=.22)
    ax.plot(y, dens_q1(y), color=BLUE, lw=2, label=r'$q_1^{\varepsilon,\delta}$')
    ax.plot(y, dens_q2(y), color=ORANGE, lw=2, label=r'$q_2^{\varepsilon,\delta}$')
    ax.plot(y, DELTA * tri(y, 3), color=GRAY, lw=1, ls=':', label=r'$\delta\,m$')
    ax.set_xlabel('$y$')
    ax.set_ylabel('density')
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(0, 1.7)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.legend(frameon=False, loc='upper right')
    ax.grid(axis='y', color='#e6e5e1', lw=0.6)
    fig.savefig(os.path.join(out, 'fig_dens_smoothed_01.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('figure written to', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else OUT)
