# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Analytically admissible maps and generators used in the accompanying tests."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .core import (Matrix, congruence, depolarizing, dimension, from_action,
                   kraus_map, matrix_unit, square_matrix)


def cp_generator(drift: ArrayLike, jumps: list[ArrayLike]) -> Matrix:
    """L(X)=A X+X A^dagger+sum K X K^dagger; exp(tL) is CP.

    No trace-preservation condition is imposed. Lie product approximation of
    exp(t(A*+*A^dagger)) and exp(t Phi) establishes the CP-semigroup property.
    """
    a = square_matrix(drift)
    d = a.shape[0]
    return np.kron(np.eye(d), a) + np.kron(a.conj(), np.eye(d)) + kraus_map(jumps, d=d)


def lindblad_generator(hamiltonian: ArrayLike, jumps: list[ArrayLike]) -> Matrix:
    h = square_matrix(hamiltonian, name="Hamiltonian")
    if not np.allclose(h, h.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError("Hamiltonian must be Hermitian")
    ks = [square_matrix(k) for k in jumps]
    if any(k.shape != h.shape for k in ks):
        raise ValueError("jump dimensions do not match the Hamiltonian")
    loss = sum((k.conj().T @ k for k in ks), np.zeros_like(h))
    return cp_generator(-1j*h - loss/2, ks)


def amplitude_damping(d: int, c: float = 0.0) -> Matrix:
    d = dimension(d)
    if d < 2 or not np.isfinite(c):
        raise ValueError("amplitude damping requires d>=2 and a finite real shift")
    v = matrix_unit(d, 0, 1)
    p = matrix_unit(d, 1, 1)
    return cp_generator(-p/2, [v]) + float(c)*np.eye(d*d)


def amplitude_damping_kraus(d: int, t: float, c: float = 0.0) -> list[Matrix]:
    d = dimension(d)
    if d < 2 or t < 0 or not np.isfinite(t + c):
        raise ValueError("requires d>=2, t>=0, and finite t,c")
    p, v = matrix_unit(d, 1, 1), matrix_unit(d, 0, 1)
    with np.errstate(over="raise", invalid="raise"):
        factor = np.exp(c*t/2)
        return [factor*(np.eye(d) - p + np.exp(-t/2)*p),
                factor*np.sqrt(-np.expm1(-t))*v]


def reduction_map(d: int, alpha: float = 0.5) -> Matrix:
    """R_alpha(X)=tr(X)I-alpha X.

    For 0<=alpha<=1/2 this is 2-positive. At alpha=1/2 and d>=3 it
    is NOT CP: its Choi matrix is I-alpha|Omega><Omega|, with eigenvalue
    1-alpha*d. The full elementary argument is in docs/MATHEMATICS.md.
    Other alpha are allowed to support negative controls.
    """
    d = dimension(d)
    if not np.isfinite(alpha):
        raise ValueError("alpha must be finite")
    return depolarizing(d) - float(alpha)*np.eye(d*d)


def non_cp_non_tp_map(d: int) -> Matrix:
    """An explicit non-TP, 2-positive, non-CP map for d>=3."""
    d = dimension(d)
    if d < 3:
        raise ValueError("requires d>=3")
    root = np.diag(np.sqrt(np.arange(1, d+1, dtype=float)))
    inv = np.diag(1 / np.diag(root))
    channel = reduction_map(d) / (d - 0.5)
    return 1.7 * congruence(inv) @ channel @ congruence(root)


def transpose_map(d: int) -> Matrix:
    return from_action(lambda x: x.T, dimension(d))


def nilpotent_cp_map(d: int) -> Matrix:
    d = dimension(d)
    shift = np.diag(np.ones(d-1), k=1)
    return kraus_map([shift])


def state_dependent_loss(rates: ArrayLike) -> Matrix:
    r = np.asarray(rates, dtype=float)
    if r.ndim != 1 or r.size == 0 or not np.isfinite(r).all() or (r < 0).any():
        raise ValueError("rates must be a nonempty finite nonnegative vector")
    return cp_generator(-np.diag(r)/2, [])


def random_kraus(d: int, rank: int, rng: np.random.Generator) -> list[Matrix]:
    d, rank = dimension(d), dimension(rank)
    return [(rng.standard_normal((d, d)) + 1j*rng.standard_normal((d, d)))
            / np.sqrt(2*d*rank) for _ in range(rank)]


def random_cp_map(d: int, rng: np.random.Generator) -> Matrix:
    return kraus_map(random_kraus(d, d+1, rng))


def random_tilted_generator(d: int, rng: np.random.Generator) -> Matrix:
    """Random GKSL jumps with positive counting weights, original loss retained."""
    d = dimension(d)
    ks = random_kraus(d, d+1, rng)
    b = rng.standard_normal((d, d)) + 1j*rng.standard_normal((d, d))
    h = (b + b.conj().T)/(2*np.sqrt(d))
    q = sum((k.conj().T @ k for k in ks), np.zeros((d, d), complex))
    weights = np.exp(rng.uniform(-1.0, 1.0, len(ks)))
    return cp_generator(-1j*h - q/2, [np.sqrt(w)*k for w, k in zip(weights, ks)])


def twisted_reduction_map(d: int) -> Matrix:
    """A TP, 2-positive, non-CP map with a traceless unitary Choi witness."""
    d = dimension(d)
    if d < 3:
        raise ValueError("requires d>=3")
    u = np.roll(np.eye(d), 1, axis=0)  # cyclic shift, so tr(U)=0
    return (depolarizing(d) - 0.5*congruence(u)) / (d - 0.5)


def two_positive_non_cp_generator(d: int, c: float = 0.3) -> Matrix:
    """exp(tL) is 2-positive but not CP at sufficiently small positive t.

    L=T-id+c*id, with T=twisted_reduction_map(d); use the exponential series.
    Its Choi derivative has a negative expectation along vec(U) orthogonal to
    vec(I), so its semigroup is genuinely outside the CP-semigroup class.
    """
    t = twisted_reduction_map(d)
    if not np.isfinite(c):
        raise ValueError("c must be finite")
    return t + (float(c)-1)*np.eye(d*d)
