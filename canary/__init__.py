"""Exact computations for "Robust Lower Bound for a Canary Call"
(V. V. Piterbarg, SSRN 5510958) and its companion paper on explicit families.

Everything is exact rational arithmetic (fractions.Fraction); floats never
enter a decision.  Modules:

  exact        measures, piece-wise linear functions, the shadow measure
               through the convex hull, a small certified exact simplex
  lower_bound  the sharp lower bound L of a single Canary, by two routes
  decider      is there a universal minimal measure for a pair of marginals
               at a given slope ratio?  (cutting-plane iteration + exact
               certificate)
  threshold    the greedy fibers and the marginals-only threshold r*
"""
