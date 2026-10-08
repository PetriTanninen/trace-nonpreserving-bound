# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Numerical realization of manuscript Lemma 5 and its epsilon perturbation.

Brouwer's theorem supplies existence in the proof; this module uses an eigensolver
instead. Positivity / 2-positivity must be established independently for input.
No clipping of negative eigenvalues, and no claimed certificate, is performed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.linalg import eig, eigvals

from .core import (Matrix, congruence, depolarizing, hermiticity_residual,
                   spectrum_distance, square_matrix, superoperator,
                   trace_preservation_residual, unvec, vec)


class NormalizationError(RuntimeError):
    """No numerically reliable faithful normalization was obtained."""


@dataclass(frozen=True)
class Normalization:
    perturbed: Matrix
    normalized: Matrix
    eigenmatrix: Matrix
    r: float
    epsilon: float
    minimum_eigenvalue_h: float
    condition_h: float
    eigenmatrix_residual: float
    trace_preservation_residual: float
    similarity_residual: float
    matched_spectrum_error: float

    def diagnostics(self) -> dict[str, float]:
        keys = ("r", "epsilon", "minimum_eigenvalue_h", "condition_h",
                "eigenmatrix_residual", "trace_preservation_residual",
                "similarity_residual", "matched_spectrum_error")
        return {k: float(getattr(self, k)) for k in keys}


def normalize_with_eigenmatrix(a: ArrayLike, h: ArrayLike, r: float, *,
                               tol: float = 1e-8, epsilon: float = 0.0) -> Normalization:
    """Implement C_H A C_H^{-1}/r and verify numerical residuals.

    a must already include any perturbation. epsilon is diagnostic metadata.
    """
    a, d = superoperator(a)
    h = square_matrix(h, name="eigenmatrix")
    if h.shape != (d, d):
        raise ValueError("eigenmatrix has the wrong dimension")
    if tol <= 0 or not np.isfinite(tol) or r <= 0 or not np.isfinite(r):
        raise ValueError("r and tol must be finite and strictly positive")
    if epsilon < 0 or not np.isfinite(epsilon):
        raise ValueError("epsilon must be finite and nonnegative")
    if hermiticity_residual(a) > tol:
        raise NormalizationError("input does not preserve Hermiticity")
    if np.linalg.norm(h-h.conj().T) > tol * max(1.0, np.linalg.norm(h)):
        raise NormalizationError("candidate eigenmatrix is not Hermitian")
    h = (h+h.conj().T)/2  # only remove verified roundoff-level skew-Hermitian part
    tr = float(np.trace(h).real)
    if tr <= 0:
        raise NormalizationError("candidate eigenmatrix has nonpositive trace")
    h = h/tr
    w, u = np.linalg.eigh(h)
    if w[0] <= 0 or w[0] <= 100*np.finfo(float).eps*w[-1]:
        raise NormalizationError("eigenmatrix is singular or numerically ill-conditioned; increase epsilon")
    root = (u*np.sqrt(w)) @ u.conj().T
    inv = (u/np.sqrt(w)) @ u.conj().T
    c, ci = congruence(root), congruence(inv)
    t = (c @ a @ ci)/r
    eigen_residual = float(np.linalg.norm(a.conj().T @ vec(h)-r*vec(h)) /
                           max(np.finfo(float).tiny, np.linalg.norm(a)*np.linalg.norm(h),
                               r*np.linalg.norm(h)))
    tp = trace_preservation_residual(t)
    sim = float(np.linalg.norm(t @ c - c @ a/r) /
                max(1.0, np.linalg.norm(c @ a/r)))
    if max(eigen_residual, tp, sim) > tol:
        raise NormalizationError(f"normalization residual too large: eigen={eigen_residual:g}, TP={tp:g}, similarity={sim:g}")
    match = spectrum_distance(eigvals(t), eigvals(a)/r)
    return Normalization(a.copy(), t, h, float(r), float(epsilon), float(w[0]),
                         float(w[-1]/w[0]), eigen_residual, tp, sim, match)


def regularized_normalization(a: ArrayLike, epsilon: float = 1e-6, *,
                              tol: float = 1e-8) -> Normalization:
    """Perturb by epsilon*Delta, solve the adjoint eigenproblem, normalize.

    Requires epsilon>0. For a known positive input the mathematical construction
    has a faithful positive eigenmatrix; numerical success is separately checked.
    This function does NOT decide positivity or 2-positivity of a general map.
    """
    a, d = superoperator(a)
    if epsilon <= 0 or not np.isfinite(epsilon):
        raise ValueError("epsilon must be finite and strictly positive")
    if tol <= 0 or not np.isfinite(tol):
        raise ValueError("tol must be finite and positive")
    p = a + epsilon*depolarizing(d)
    values, vectors = eig(p.conj().T)
    k = int(np.argmax(values.real))
    lam = values[k]
    if lam.real <= 0 or abs(lam.imag) > tol*max(1.0, abs(lam)):
        raise NormalizationError("dominant eigenvalue is not positive real within tolerance")
    h = unvec(vectors[:, k], d)
    tr = np.trace(h)
    if abs(tr) < 100*np.finfo(float).eps*np.linalg.norm(h):
        raise NormalizationError("dominant eigenvector has numerically zero trace")
    h = h/tr  # fixes the arbitrary complex phase without taking entrywise absolute values
    return normalize_with_eigenmatrix(p, h, float(lam.real), tol=tol, epsilon=epsilon)
