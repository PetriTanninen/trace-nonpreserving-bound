# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
import numpy as np
import pytest

from spectral_bounds import bound_report, congruence
from spectral_bounds.examples import nilpotent_cp_map, non_cp_non_tp_map, random_cp_map
from spectral_bounds.normalization import (NormalizationError,
    normalize_with_eigenmatrix, regularized_normalization)


@pytest.mark.parametrize("d", range(2,7))
@pytest.mark.parametrize("epsilon", [1e-2,1e-5,1e-8])
def test_random_normalization(d,epsilon):
    p = random_cp_map(d,np.random.default_rng(d))
    n = regularized_normalization(p,epsilon)
    assert n.minimum_eigenvalue_h > 0
    assert n.trace_preservation_residual < 1e-8
    assert n.eigenmatrix_residual < 1e-10
    assert n.similarity_residual < 1e-10
    assert n.matched_spectrum_error < 1e-8
    before, after = bound_report(n.perturbed), bound_report(n.normalized)
    assert before.maximum_real == pytest.approx(n.r)
    assert before.slack/n.r == pytest.approx(after.slack)
    assert after.maximum_real == pytest.approx(1)


@pytest.mark.parametrize("epsilon", [1e-1,1e-3,1e-5,1e-8])
@pytest.mark.parametrize("family", ["zero","nilpotent","noncp"])
def test_singular_and_noncp_limits(family,epsilon):
    p = {"zero": np.zeros((9,9)), "nilpotent": nilpotent_cp_map(3),
         "noncp": non_cp_non_tp_map(3)}[family]
    n = regularized_normalization(p,epsilon)
    assert n.minimum_eigenvalue_h > 0
    assert bound_report(n.perturbed).holds_within_tolerance
    assert bound_report(n.normalized).holds_within_tolerance
    assert n.trace_preservation_residual < 1e-8


def test_complex_faithful_eigenmatrix():
    h = np.array([[2, .2+.3j],[.2-.3j,1]],complex)
    r = 1.7
    # Identity superoperator has a degenerate eigenvalue but any H>0 works.
    n = normalize_with_eigenmatrix(r*np.eye(4),h,r)
    assert np.allclose(n.normalized,np.eye(4))
    assert np.allclose(n.eigenmatrix,h/np.trace(h))


def test_failures_are_explicit():
    with pytest.raises(ValueError): regularized_normalization(np.eye(4),0)
    with pytest.raises(ValueError): regularized_normalization(np.eye(4),-1)
    with pytest.raises(ValueError): normalize_with_eigenmatrix(np.eye(4),np.eye(2),0)
    with pytest.raises(NormalizationError):
        normalize_with_eigenmatrix(np.eye(4),np.diag([1,0]),1)
    with pytest.raises(NormalizationError):
        normalize_with_eigenmatrix(np.eye(4),np.array([[1,1],[0,1]]),1)
    with pytest.raises(NormalizationError):
        normalize_with_eigenmatrix(2*np.eye(4),np.eye(2),1)
