# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Exact algebra for sharpness, non-CP witnesses, and the loss counterexample.

These checks are not a formalization of the all-dimensional main theorem.
"""
from __future__ import annotations

import sympy as sp


def exact_report(max_dimension: int = 6) -> dict:
    if not __debug__:
        raise RuntimeError("Run without -O: verification assertions must be enabled")
    if not isinstance(max_dimension, int) or max_dimension < 2:
        raise ValueError("max_dimension must be an integer >=2")
    d = sp.Symbol("d", integer=True, positive=True)
    c, a = sp.symbols("c a", real=True)
    multiplicities = [(d-1)**2, 2*(d-1), 1]
    trace_g = -sp.Rational(1, 2)*multiplicities[1] - multiplicities[2]
    identities = {
        "multiplicity_sum_minus_d_squared": sp.expand(sum(multiplicities)-d*d),
        "trace_G_plus_d": sp.expand(trace_g+d),
        "shifted_equality_slack": sp.expand(d*(c-1)+(d*d-d)*c-(c*d*d-d)),
        "replacement_rhs_minus_trace": sp.expand(a*(c-1)+(d*d-a)*c-(c*d*d-d)),
    }
    assert identities["replacement_rhs_minus_trace"] == d-a
    assert all(v == 0 for k, v in identities.items() if k != "replacement_rhs_minus_trace")
    checked = []
    for n in range(2, max_dimension+1):
        def unit(i, j):
            x = sp.zeros(n); x[i, j] = 1
            return x
        v, p = unit(0, 1), unit(1, 1)
        def action(x):
            return v*x*v.T - (p*x+x*p)/2
        pairs = [(unit(i, j), sp.Integer(0)) for j in range(n) for i in range(n)
                 if i != 1 and j != 1]
        pairs += [(unit(1, j), -sp.Rational(1, 2)) for j in range(n) if j != 1]
        pairs += [(unit(j, 1), -sp.Rational(1, 2)) for j in range(n) if j != 1]
        pairs += [(unit(1, 1)-unit(0, 0), -sp.Integer(1))]
        assert len(pairs) == n*n
        assert all(action(x) == lam*x for x, lam in pairs)
        basis = sp.Matrix.hstack(*[sp.Matrix([x[i,j] for j in range(n) for i in range(n)])
                                  for x, _ in pairs])
        assert basis.rank() == n*n
        checked.append(n)
    q = sp.Symbol("q", nonnegative=True)
    p = sp.diag(0, 1)
    v = sp.Matrix([[0, 1], [0, 0]])
    k0 = sp.eye(2)-p+q*p
    # q=exp(-t/2) in [0,1], so K1^dag K1=(1-q^2)P exactly.
    assert sp.simplify(k0.T*k0+(1-q*q)*(v.T*v)-sp.eye(2)) == sp.zeros(2)
    # Pure state-dependent loss K=diag(0,2), eigenvalue on E_ij=-(k_i+k_j)/2.
    loss_spectrum = [sp.Integer(0), -sp.Integer(1), -sp.Integer(1), -sp.Integer(2)]
    return {
        "all_checks_passed": True,
        "scope": "Exact algebraic identities; finite-dimension eigenbasis checks, not a formal proof of the theorem.",
        "identities": {k: str(v) for k, v in identities.items()},
        "eigenbasis_dimensions_checked": checked,
        "kraus_completeness_in_q": "zero matrix",
        "state_dependent_loss_spectrum": [str(v) for v in loss_spectrum],
        "reduction_choi_eigenvalue_at_d3_alpha_half": str(1-sp.Rational(1, 2)*3),
    }
