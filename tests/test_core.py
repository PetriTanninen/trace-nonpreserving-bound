# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
import numpy as np
import pytest
from scipy.linalg import expm

from spectral_bounds.core import (apply, bound_report, choi, choi_min_eigenvalue,
    congruence, depolarizing, from_action, kraus_map, matrix_unit, spectrum_distance,
    superoperator, trace_preservation_residual, unvec, vec)
from spectral_bounds.examples import (amplitude_damping, amplitude_damping_kraus,
    lindblad_generator, nilpotent_cp_map, non_cp_non_tp_map, random_cp_map,
    random_tilted_generator, reduction_map, state_dependent_loss, transpose_map,
    two_positive_non_cp_generator)


@pytest.mark.parametrize("d", range(1, 7))
def test_vec_and_action(d):
    rng = np.random.default_rng(101+d)
    a, x, b = [rng.normal(size=(d,d)) + 1j*rng.normal(size=(d,d)) for _ in range(3)]
    assert np.allclose(unvec(vec(x), d), x)
    assert np.allclose(vec(a@x@b), np.kron(b.T, a)@vec(x))
    s = from_action(lambda y: a@y@a.conj().T, d)
    assert np.allclose(s, congruence(a))
    assert np.allclose(apply(s,x), a@x@a.conj().T)
    assert np.allclose(apply(depolarizing(d),x), np.trace(x)*np.eye(d))
    # Hilbert-Schmidt adjoint and Choi ordering: no hidden transpose convention.
    assert np.allclose(np.vdot(a,apply(s,x)), np.vdot(apply(s.conj().T,a),x))
    assert np.allclose(choi(s), np.outer(vec(a), vec(a).conj()))
    assert np.allclose(choi(depolarizing(d)), np.eye(d*d))


@pytest.mark.parametrize("d", range(2, 7))
@pytest.mark.parametrize("seed", range(5))
def test_random_cp_and_generators(d,seed):
    rng = np.random.default_rng(seed+100*d)
    p = random_cp_map(d,rng)
    assert choi_min_eigenvalue(p) > -1e-10
    assert bound_report(p).holds_within_tolerance
    g = random_tilted_generator(d,rng)
    assert bound_report(g).holds_within_tolerance
    assert choi_min_eigenvalue(expm(.2*g)) > -1e-10


@pytest.mark.parametrize("d", range(2, 8))
@pytest.mark.parametrize("c", [-2.0, 0.0, .3, 2.0])
def test_sharpness(d,c):
    g = amplitude_damping(d,c)
    r = bound_report(g)
    assert r.trace == pytest.approx(c*d*d-d)
    assert r.minimum_real == pytest.approx(c-1)
    assert r.maximum_real == pytest.approx(c)
    assert abs(r.slack) <= r.tolerance
    for t in (0.0, .03, .5, 2.0):
        channel = kraus_map(amplitude_damping_kraus(d,t,c))
        assert np.allclose(expm(t*g),channel,rtol=1e-11,atol=1e-11)
        assert trace_preservation_residual(channel) == pytest.approx(abs(np.expm1(c*t)),abs=1e-10)
    a = d+.3
    assert a*r.minimum_real+(d*d-a)*r.maximum_real < r.trace-0.2


@pytest.mark.parametrize("d", range(3, 7))
def test_noncp_families(d):
    r = reduction_map(d)
    assert choi_min_eigenvalue(r) == pytest.approx(1-d/2)
    phi = non_cp_non_tp_map(d)
    assert choi_min_eigenvalue(phi) < 0
    assert trace_preservation_residual(phi) > 0.1
    assert bound_report(phi).holds_within_tolerance
    g = two_positive_non_cp_generator(d)
    assert bound_report(g).holds_within_tolerance
    assert choi_min_eigenvalue(expm(1e-4*g)) < -1e-7
    # Positivity at k=2 is established analytically in MATHEMATICS.md,
    # not inferred from random state samples or the tested inequality.


def test_controls_and_state_dependent_loss():
    r = bound_report(transpose_map(2))
    assert r.slack == pytest.approx(-2)
    assert not r.holds_within_tolerance
    g = state_dependent_loss([0,2])
    assert sorted(np.linalg.eigvals(g).real) == [-2,-1,-1,0]
    assert bound_report(g).slack == pytest.approx(0)
    d = np.diag([1,np.exp(-.4)])
    assert np.allclose(expm(.4*g),congruence(d))
    assert trace_preservation_residual(expm(.4*g)) > .1


@pytest.mark.parametrize("d", range(1,7))
def test_zero_nilpotent_and_scalar(d):
    for p in [np.zeros((d*d,d*d)), nilpotent_cp_map(d), -.7*np.eye(d*d)]:
        # A negative scalar identity is a CP-semigroup generator, not a positive map.
        assert bound_report(p).slack == pytest.approx(0)
    p = nilpotent_cp_map(d)
    assert np.allclose(np.linalg.matrix_power(p,d),0)


def test_nonreal_spectrum_and_spectral_matching():
    h = np.diag([0,2])
    g = lindblad_generator(h,[])
    r = bound_report(g)
    assert abs(np.linalg.eigvals(g).imag).max() == pytest.approx(2)
    assert r.minimum_real == pytest.approx(0)
    assert r.maximum_real == pytest.approx(0)
    assert r.slack == pytest.approx(0)
    assert spectrum_distance([1,1j,-1j],[1j,1,-1j]) == 0


@pytest.mark.parametrize("bad", [np.zeros((3,3)), np.zeros((2,3)), np.eye(4)*np.nan, [], [1,2]])
def test_reject_invalid_superoperators(bad):
    with pytest.raises(ValueError):
        superoperator(bad)


def test_validation_and_non_hp_rejection():
    with pytest.raises(ValueError): bound_report(1j*np.eye(4))
    with pytest.raises(ValueError): bound_report(np.eye(4),rtol=-1)
    with pytest.raises(ValueError): choi_min_eigenvalue(1j*np.eye(4))
    with pytest.raises(ValueError): kraus_map([])
    with pytest.raises(ValueError): apply(np.eye(4),np.eye(3))
    with pytest.raises(ValueError): unvec(np.ones(5),2)
    with pytest.raises(ValueError): matrix_unit(2,-1,0)
    with pytest.raises(ValueError): depolarizing(True)
    with pytest.raises(ValueError): lindblad_generator([[0,1],[0,0]],[])
    with pytest.raises(ValueError): state_dependent_loss([-1,2])
