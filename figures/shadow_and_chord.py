"""Figures 2 and 7 of "Robust Lower Bound for a Canary Call":
  fig_shadow_hull_01.png -- the shadow through the hull (Example in Section 5.3)
  fig_chord_01.png       -- the chord sub-hedge (Theorem 6.5)
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


def call(meas, k):
    return sum(m * max(y - k, F(0)) for y, m in meas)


def hull_pl(xs, ys):
    """lower convex hull of the points (xs, ys) (xs sorted); returns the
    indices of the hull vertices (Andrew's monotone chain, lower part)."""
    pts = list(zip(xs, ys))
    hull = []
    for i, p in enumerate(pts):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = pts[hull[-2]], pts[hull[-1]]
            if (y2 - y1) * (p[0] - x1) - (p[1] - y1) * (x2 - x1) >= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    return hull


def main():
    os.makedirs(OUT, exist_ok=True)
    matplotlib.rcdefaults()
    plt.rcParams.update(RC)
    # ------------------------------------------------------------ shadow
    q2 = [(F(-2), F(1, 4)), (F(0), F(1, 2)), (F(2), F(1, 4))]
    nu = [(F(-1), F(1, 4))]
    nodes = [F(-3), F(-2), F(-1), F(0), F(2), F(3)]
    g = [call(q2, k) - call(nu, k) for k in nodes]
    hidx = hull_pl(nodes, g)
    hull_pts = [(nodes[i], g[i]) for i in hidx]
    sh = {}
    for a, b, c in zip(hull_pts, hull_pts[1:], hull_pts[2:]):
        sl1 = (b[1] - a[1]) / (b[0] - a[0]); sl2 = (c[1] - b[1]) / (c[0] - b[0])
        if sl2 - sl1 > 0:
            sh[b[0]] = sl2 - sl1
    eta = {y: m - sh.get(y, F(0)) for y, m in q2}
    assert eta == {F(-2): F(1, 8), F(0): F(1, 8), F(2): F(0)}, eta
    assert {k: v for k, v in sh.items()} == {F(-2): F(1, 8), F(0): F(3, 8), F(2): F(1, 4)}, sh

    ks = np.linspace(-3, 3, 601)
    Cq2 = [float(call(q2, F(k).limit_denominator(10**6))) for k in ks]
    Cnu = [float(call(nu, F(k).limit_denominator(10**6))) for k in ks]
    gg = [a - b for a, b in zip(Cq2, Cnu)]
    hx = [float(x) for x, _ in hull_pts]; hy = [float(y) for _, y in hull_pts]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                                   gridspec_kw={'hspace': 0.12})
    ax1.plot(ks, Cq2, color=BLUE, lw=2, label=r'$C_{q_2}$')
    ax1.plot(ks, Cnu, color=ORANGE, lw=2, label=r'$C_{\nu}$')
    ax1.set_ylabel('call function')
    ax1.legend(frameon=False, loc='upper right')
    ax1.grid(axis='y', color='#e6e5e1', lw=0.6)
    ax2.plot(ks, gg, color=INK, lw=2, label=r'$g=C_{q_2}-C_{\nu}$')
    ax2.plot(hx, hy, color=GRAY, lw=2, ls='--', label=r'$\hat g=\lfloor g\rfloor$')
    for x in (-2.0, 0.0, 2.0):
        ax2.plot([x], [float(call(q2, F(x)) - call(nu, F(x)))], 'o', color=GRAY,
                 ms=7, mfc='white', mew=2)
    ax2.annotate('contact set $K$', xy=(0, float(g[3])), xytext=(0.35, 0.55),
                 color=GRAY, arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.8))
    ax2.set_xlabel('strike $k$')
    ax2.set_ylabel('difference')
    ax2.legend(frameon=False, loc='upper right')
    ax2.grid(axis='y', color='#e6e5e1', lw=0.6)
    ax2.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    fig.savefig(os.path.join(OUT, 'fig_shadow_hull_01.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------------ chord
    # e2(y) = (y - k2)^+ with k2 inside a gap (alpha, beta) of the agreement
    # set Z; alpha, beta, k2 are chosen, rho = max(e2, mu) is computed exactly
    alpha, beta, k2 = F(-1), F(3, 2), F(1, 2)
    s2 = F(1)
    e2 = lambda y: s2 * max(y - k2, F(0))
    slope = (e2(beta) - e2(alpha)) / (beta - alpha)
    mu = lambda y: e2(alpha) + slope * (y - alpha)
    assert F(0) <= slope <= s2
    ys = np.linspace(-2.5, 3, 551)
    E2 = [float(e2(F(y).limit_denominator(10**6))) for y in ys]
    MU = [float(mu(F(y).limit_denominator(10**6))) for y in ys]
    RHO = [max(a, b) for a, b in zip(E2, MU)]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(ys, E2, color=BLUE, lw=2, label=r'$e_2$')
    ax.plot(ys, MU, color=GRAY, lw=1.6, ls='--', label=r'$\mu$ (chord)')
    ax.plot(ys, RHO, color=ORANGE, lw=2.6, alpha=0.9, label=r'$\rho=\max(e_2,\mu)$')
    ax.plot([-2.5, float(alpha)], [-0.08, -0.08], color=INK, lw=4, solid_capstyle='butt')
    ax.plot([float(beta), 3], [-0.08, -0.08], color=INK, lw=4, solid_capstyle='butt')
    ax.text(-2.4, -0.3, r'$Z$', color=INK)
    ax.text(2.7, -0.3, r'$Z$', color=INK)
    for x, lab in ((alpha, r'$\alpha$'), (beta, r'$\beta$'), (k2, r'$k_2$')):
        ax.axvline(float(x), color='#e6e5e1', lw=0.8)
        ax.text(float(x), 2.45, lab, ha='center', color=INK)
    ax.set_ylim(-0.4, 2.6)
    ax.set_xlabel('$y$')
    ax.set_yticks([0, 1, 2])
    ax.legend(frameon=False, loc='upper left')
    fig.savefig(os.path.join(OUT, 'fig_chord_01.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('figures written to', OUT, '; shadow =', {str(k): str(v) for k, v in eta.items()},
          'chord slope =', slope)


if __name__ == '__main__':
    main()
