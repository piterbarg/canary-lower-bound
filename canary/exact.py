"""Exact-rational building blocks: finite measures, piece-wise linear
functions, the shadow measure computed through the convex hull, and a small
exact simplex whose optimality is certified before a value is returned.

Everything is fractions.Fraction end to end; floats never enter any
decision path.  Paper references are to the SSRN version of November 2026
of "Robust Lower Bound for a Canary Call":

  shadow()   Proposition 5.6 (the shadow through the hull): the shadow of
             nu in q2 is the measure whose call function is C_q2 minus the
             largest convex minorant of g = C_q2 - C_nu, with the left ray
             slope fixed to -(mass q2 - mass nu) and the right ray slope 0.

A note on node sets: the call function of a finitely supported measure is
piece-wise linear with kinks only at its atoms, so exact comparisons of two
call functions need only be made at the union of the two supports plus the
two ray slopes.  Every function below that takes a node set says so.
"""

from fractions import Fraction as F


class InvariantError(Exception):
    """A state the theory proves unreachable was reached.  This is a bug in
    the implementation, never a verdict about the instance."""


class ScopeError(Exception):
    """Input outside the scope of the routine (not in convex order, slope
    ratio below one, wrong number of atoms)."""


class ShadowCertError(Exception):
    """The computed shadow failed its exact certificate."""


class LPError(Exception):
    pass


# =========================================================================
# measures: sorted list of (point, mass) pairs, Fractions, masses > 0
# =========================================================================

def measure(pairs):
    """Normalise to sorted [(z, m), ...] with zero masses dropped.
    Raises on floats, negative mass or duplicate points."""
    acc = {}
    for z, m in pairs:
        if isinstance(z, float) or isinstance(m, float):
            raise TypeError('float rejected in measure -- pass exact '
                            'rationals')
        z, m = F(z), F(m)
        if m < 0:
            raise ValueError('negative mass %s at %s' % (m, z))
        if z in acc:
            raise ValueError('duplicate point %s' % z)
        acc[z] = m
    return sorted((z, m) for z, m in acc.items() if m > 0)


def mass(mu):
    return sum(m for _, m in mu)


def mean_sum(mu):
    """sum of z*m (the un-normalised first moment)."""
    return sum(z * m for z, m in mu)


def supp(mu):
    return [z for z, _ in mu]


def cpot(mu, k):
    """the call function C_mu(k) = sum m (z-k)^+ (Section 5.2)."""
    k = F(k)
    return sum(m * (z - k) for z, m in mu if z > k)


def cx_leq(q1, q2, nodes):
    """exact convex-order test q1 <=cx q2 on the declared node set `nodes`
    (must be supp q1 union supp q2: both call functions are piece-wise
    linear with kinks only there, and their rays agree once mass and mean
    agree, so node checks decide the order on all of R).
    Returns (ok, why)."""
    if mass(q1) != mass(q2):
        return False, 'mass %s != %s' % (mass(q1), mass(q2))
    if mean_sum(q1) != mean_sum(q2):
        return False, 'mean-sum %s != %s' % (mean_sum(q1), mean_sum(q2))
    for z in nodes:
        a, b = cpot(q1, z), cpot(q2, z)
        if a > b:
            return False, 'C_q1(%s)=%s > C_q2(%s)=%s' % (z, a, z, b)
    return True, ''


def ambient_nodes(q1, q2):
    """supp q1 union supp q2, sorted."""
    return sorted(set(supp(q1)) | set(supp(q2)))


# =========================================================================
# piecewise-linear functions on R:
#   breakpoints xs (sorted), values ys at xs, left/right ray slopes.
#   Affine between consecutive breakpoints; at least one breakpoint.
# =========================================================================

class PL:
    __slots__ = ('xs', 'ys', 'sL', 'sR')

    def __init__(self, xs, ys, sL, sR):
        if not (len(xs) == len(ys) and xs):
            raise ValueError('PL needs >= 1 breakpoint')
        if not all(xs[i] < xs[i + 1] for i in range(len(xs) - 1)):
            raise ValueError('PL breakpoints must be strictly sorted')
        self.xs = [F(x) for x in xs]
        self.ys = [F(y) for y in ys]
        self.sL = F(sL)
        self.sR = F(sR)

    @staticmethod
    def const(c, at=F(0)):
        return PL([F(at)], [F(c)], F(0), F(0))

    @staticmethod
    def from_calls(mu, extra_nodes=()):
        """C_mu as a PL: breakpoints at supp mu (plus any extra nodes);
        sL = -mass, sR = 0."""
        pts = sorted(set(supp(mu)) | set(F(z) for z in extra_nodes))
        if not pts:
            return PL.const(0)
        return PL(pts, [cpot(mu, z) for z in pts], -mass(mu), F(0))

    def __call__(self, k):
        k = F(k)
        xs, ys = self.xs, self.ys
        if k <= xs[0]:
            return ys[0] + self.sL * (k - xs[0])
        if k >= xs[-1]:
            return ys[-1] + self.sR * (k - xs[-1])
        lo, hi = 0, len(xs) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= k:
                lo = mid
            else:
                hi = mid
        t = (k - xs[lo]) / (xs[hi] - xs[lo])
        return ys[lo] + t * (ys[hi] - ys[lo])

    def slopes(self):
        """the slope sequence: [sL, cell slopes..., sR]."""
        out = [self.sL]
        for i in range(len(self.xs) - 1):
            out.append((self.ys[i + 1] - self.ys[i])
                       / (self.xs[i + 1] - self.xs[i]))
        out.append(self.sR)
        return out

    def slope_left_of(self, x):
        x = F(x)
        sl = self.slopes()
        for i, bp in enumerate(self.xs):
            if bp >= x:
                return sl[i]
        return self.sR

    def slope_right_of(self, x):
        x = F(x)
        sl = self.slopes()
        for i in range(len(self.xs) - 1, -1, -1):
            if self.xs[i] <= x:
                return sl[i + 1]
        return self.sL

    def kink_at(self, x):
        """slope jump (right minus left) at x; 0 if affine through x."""
        return self.slope_right_of(x) - self.slope_left_of(x)

    def refit(self, new_xs):
        """same function on a superset of breakpoints."""
        new_xs = sorted(set([F(x) for x in new_xs]) | set(self.xs))
        return PL(new_xs, [self(x) for x in new_xs], self.sL, self.sR)

    def _zip(self, other):
        return sorted(set(self.xs) | set(other.xs))

    def add(self, other):
        xs = self._zip(other)
        return PL(xs, [self(x) + other(x) for x in xs],
                  self.sL + other.sL, self.sR + other.sR)

    def sub(self, other):
        xs = self._zip(other)
        return PL(xs, [self(x) - other(x) for x in xs],
                  self.sL - other.sL, self.sR - other.sR)

    def add_const(self, c):
        return PL(self.xs, [y + F(c) for y in self.ys], self.sL, self.sR)

    def scale(self, c):
        c = F(c)
        if c < 0:
            raise ValueError('scale by nonnegative only')
        return PL(self.xs, [y * c for y in self.ys], self.sL * c,
                  self.sR * c)

    def max_with(self, other):
        """pointwise max, exact: merges breakpoints, then inserts every
        interior-cell crossing and any ray crossing as new breakpoints."""
        xs = self._zip(other)
        d = self.sub(other).refit(xs)
        cross = []
        for i in range(len(d.xs) - 1):
            a, b = d.xs[i], d.xs[i + 1]
            va, vb = d.ys[i], d.ys[i + 1]
            if (va < 0 < vb) or (vb < 0 < va):
                cross.append(a + (b - a) * (-va) / (vb - va))
        if d.sL != 0:
            t = d.xs[0] - d.ys[0] / d.sL
            if t < d.xs[0]:
                cross.append(t)
        if d.sR != 0:
            t = d.xs[-1] - d.ys[-1] / d.sR
            if t > d.xs[-1]:
                cross.append(t)
        xs2 = sorted(set(xs) | set(cross))
        ys2 = [max(self(x), other(x)) for x in xs2]
        dL = self.sL - other.sL
        dR = self.sR - other.sR
        vL = self(xs2[0]) - other(xs2[0])
        if vL > 0 or (vL == 0 and dL < 0):
            sL = self.sL
        elif vL < 0 or (vL == 0 and dL > 0):
            sL = other.sL
        else:
            sL = self.sL
        vR = self(xs2[-1]) - other(xs2[-1])
        if vR > 0 or (vR == 0 and dR > 0):
            sR = self.sR
        elif vR < 0 or (vR == 0 and dR < 0):
            sR = other.sR
        else:
            sR = self.sR
        return PL(xs2, ys2, sL, sR)

    def le(self, other):
        """exact f <= g on all of R: node values on the merged set plus
        the two ray-slope conditions."""
        d = other.sub(self)
        if any(y < 0 for y in d.ys):
            return False
        return d.sL <= 0 and d.sR >= 0

    def eq(self, other):
        d = other.sub(self)
        return (all(y == 0 for y in d.ys) and d.sL == 0 and d.sR == 0)

    def breakpoints_with_kinks(self):
        """[(x, jump)] at breakpoints where the slope actually jumps."""
        sl = self.slopes()
        out = []
        for i, x in enumerate(self.xs):
            j = sl[i + 1] - sl[i]
            if j != 0:
                out.append((x, j))
        return out


# =========================================================================
# the shadow through the hull (Proposition 5.6): eta = S^{q2}(nu),
# C_eta = C_q2 - ghat, ghat the largest convex minorant of g = C_q2 - C_nu
# on R (left ray slope -d, d = mass(q2) - mass(nu); right ray slope 0).
# =========================================================================

def shadow(nu, q2):
    """the shadow of nu in q2, computed on the node set supp(nu) union
    supp(q2) (the kink set of g, which suffices for the hull on R).

    Returns eta as a measure.  Raises ShadowCertError unless eta passes the
    exact certificate: 0 <= eta <= q2, mass(eta) = mass(nu),
    mean-sum(eta) = mean-sum(nu), nu <=cx eta, and mass is removed from q2
    only at nodes where the hull touches g."""
    nodes = sorted(set(supp(nu)) | set(supp(q2)))
    q2d = dict(q2)
    g = [cpot(q2, z) - cpot(nu, z) for z in nodes]
    d = mass(q2) - mass(nu)
    if d < 0:
        raise ShadowCertError('mass(nu)=%s > mass(q2)=%s -- no shadow'
                              % (mass(nu), mass(q2)))
    n = len(nodes)
    if n == 1:
        raw = [(nodes[0], mass(nu))]
        _certify_shadow(raw, nu, q2)
        return measure(raw)
    # candidate slopes for the sup of affine minorants: chord slopes
    # clipped to the admissible ray range [-d, 0], plus the two ray slopes
    chords = set()
    for j in range(n):
        for l in range(j + 1, n):
            chords.add((g[l] - g[j]) / (nodes[l] - nodes[j]))
    cands = sorted(set([max(-d, min(F(0), s)) for s in chords] + [-d, F(0)]))
    ghat = []
    for i in range(n):
        best = None
        for b in cands:
            v = min(g[j] + b * (nodes[i] - nodes[j]) for j in range(n))
            if best is None or v > best:
                best = v
        ghat.append(best)
    # second differences of the hull, with ray slopes -d (left) and 0 (right)
    jumps = []
    for i in range(n):
        sl = ((ghat[i] - ghat[i - 1]) / (nodes[i] - nodes[i - 1])
              if i > 0 else -d)
        sr = ((ghat[i + 1] - ghat[i]) / (nodes[i + 1] - nodes[i])
              if i < n - 1 else F(0))
        jumps.append(sr - sl)
    raw = [(z, q2d.get(z, F(0)) - j) for z, j in zip(nodes, jumps)]
    _certify_shadow(raw, nu, q2)
    # contact check: mass is removed at a node only where the hull touches g
    for i, z in enumerate(nodes):
        if raw[i][1] < q2d.get(z, F(0)) and ghat[i] != g[i]:
            raise ShadowCertError(
                'mass removed at %s where ghat(%s)=%s != g(%s)=%s'
                % (z, z, ghat[i], z, g[i]))
    return measure(raw)


def _certify_shadow(raw, nu, q2):
    """the exact shadow certificate on the raw node/mass list:
    0 <= eta <= q2, mass, mean, and nu <=cx eta at every node."""
    q2d = dict(q2)
    for z, m in raw:
        if m < 0:
            raise ShadowCertError('eta < 0 at %s' % z)
        if m > q2d.get(z, F(0)):
            raise ShadowCertError(
                'eta(%s)=%s > q2(%s)=%s' % (z, m, z, q2d.get(z, F(0))))
    tot = sum(m for _, m in raw)
    if tot != mass(nu):
        raise ShadowCertError(
            'mass(eta)=%s != mass(nu)=%s' % (tot, mass(nu)))
    ms = sum(z * m for z, m in raw)
    if ms != mean_sum(nu):
        raise ShadowCertError(
            'mean-sum(eta)=%s != mean-sum(nu)=%s' % (ms, mean_sum(nu)))
    for z, _ in raw:
        pe = sum(mm * (zz - z) for zz, mm in raw if zz > z)
        if pe < cpot(nu, z):
            raise ShadowCertError(
                'nu <=cx eta fails at %s: %s < %s'
                % (z, pe, cpot(nu, z)))


# =========================================================================
# a small exact simplex: min c.x  s.t.  A_eq x = b_eq, A_ub x <= b_ub,
# 0 <= x <= ub (ub entries may be None).  Dense, Fractions, Bland's rule.
# =========================================================================

def solve_lp(c, A_eq, b_eq, A_ub, b_ub, ub):
    """returns (x, val).  Exact.  Raises LPError on infeasible/unbounded.
    Optimality is certified before returning: primal feasibility is
    re-checked exactly, every nonbasic reduced cost is >= 0, and the
    objective equals the tableau's tracked value (zero duality gap)."""
    c = [F(v) for v in c]
    n = len(c)
    rows = []
    rhs = []
    for a, b in zip(A_eq, b_eq):
        rows.append([F(v) for v in a] + []), rhs.append(F(b))
    nub = 0
    for a, b in zip(A_ub, b_ub):
        rows.append([F(v) for v in a]), rhs.append(F(b))
        nub += 1
    for j, u in enumerate(ub):
        if u is not None:
            r = [F(0)] * n
            r[j] = F(1)
            rows.append(r), rhs.append(F(u))
            nub += 1
    m = len(rows)
    nslack = nub
    # tableau columns: x (n) | slacks (nslack) | artificials (m)
    ncol = n + nslack + m
    T = []
    si = 0
    neq = len(list(A_eq))
    for i in range(m):
        row = rows[i] + [F(0)] * (nslack + m)
        if i >= neq:                       # inequality row: slack
            row[n + si] = F(1)
            si += 1
        T.append(row + [rhs[i]])
    for i in range(m):
        if T[i][-1] < 0:
            T[i] = [-v for v in T[i]]
    for i in range(m):
        T[i][n + nslack + i] = F(1)
    basis = [n + nslack + i for i in range(m)]

    def pivot(T, basis, obj, banned=()):
        banned = set(banned)
        while True:
            # Bland's rule: smallest-index entering with negative reduced cost
            enter = -1
            for j in range(ncol):
                if j not in banned and obj[j] < 0:
                    enter = j
                    break
            if enter < 0:
                return
            leave, best = -1, None
            for i in range(m):
                if T[i][enter] > 0:
                    ratio = T[i][-1] / T[i][enter]
                    if best is None or ratio < best or \
                       (ratio == best and basis[i] < basis[leave]):
                        best, leave = ratio, i
            if leave < 0:
                raise LPError('unbounded')
            piv = T[leave][enter]
            T[leave] = [v / piv for v in T[leave]]
            for i in range(m):
                if i != leave and T[i][enter] != 0:
                    f = T[i][enter]
                    T[i] = [a - f * b for a, b in zip(T[i], T[leave])]
            f = obj[enter]
            obj[:] = [a - f * b for a, b in zip(obj, T[leave])]
            basis[leave] = enter

    # phase 1
    obj1 = [F(0)] * (ncol + 1)
    for j in range(n + nslack, ncol):
        obj1[j] = F(1)
    for i in range(m):
        obj1 = [a - b for a, b in zip(obj1, T[i])]
    pivot(T, basis, obj1)
    if any(basis[i] >= n + nslack and T[i][-1] != 0 for i in range(m)):
        raise LPError('infeasible')
    # drive residual artificials out of the basis where possible
    for i in range(m):
        if basis[i] >= n + nslack:
            for j in range(n + nslack):
                if T[i][j] != 0:
                    piv = T[i][j]
                    T[i] = [v / piv for v in T[i]]
                    for k2 in range(m):
                        if k2 != i and T[k2][j] != 0:
                            f = T[k2][j]
                            T[k2] = [a - f * b
                                     for a, b in zip(T[k2], T[i])]
                    basis[i] = j
                    break
    # phase 2: no big-M; artificial columns are simply barred from entering
    art0 = n + nslack
    obj2 = [F(v) for v in c] + [F(0)] * (nslack + m) + [F(0)]
    for i in range(m):
        if obj2[basis[i]] != 0:
            f = obj2[basis[i]]
            obj2 = [a - f * b for a, b in zip(obj2, T[i])]
    pivot(T, basis, obj2, banned=range(art0, ncol))
    x = [F(0)] * n
    for i in range(m):
        if basis[i] < n:
            x[basis[i]] = T[i][-1]
    # exact certification: primal feasibility + reduced costs + duality gap
    for a, b in zip(A_eq, b_eq):
        if sum(F(ai) * xi for ai, xi in zip(a, x)) != F(b):
            raise LPError('certify: eq violated')
    for a, b in zip(A_ub, b_ub):
        if sum(F(ai) * xi for ai, xi in zip(a, x)) > F(b):
            raise LPError('certify: ub violated')
    for j, u in enumerate(ub):
        if x[j] < 0 or (u is not None and x[j] > F(u)):
            raise LPError('certify: bound violated')
    for j in range(art0):
        if obj2[j] < 0:
            raise LPError('certify: negative reduced cost survives')
    val = sum(ci * xi for ci, xi in zip(c, x))
    if val != -obj2[-1]:
        raise LPError('certify: duality gap %s' % (val + obj2[-1]))
    return x, val
