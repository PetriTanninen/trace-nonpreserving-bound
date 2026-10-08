# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Finite-dimensional superoperators, with column-stacked vectorization.

vec(A X B) = (B.T kron A) vec(X).  A superoperator is NEVER assumed
Hermitian as a d**2 by d**2 matrix, even when it preserves Hermiticity.
The checks here are floating-point diagnostics, not positivity certificates.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from math import isqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment

Matrix = NDArray[np.complex128]


def dimension(d: int) -> int:
    if isinstance(d, (bool, np.bool_)) or not isinstance(d, (int, np.integer)) or d < 1:
        raise ValueError("d must be a positive integer")
    return int(d)


def square_matrix(a: ArrayLike, *, name: str = "matrix") -> Matrix:
    a = np.asarray(a, dtype=np.complex128)
    if a.ndim != 2 or a.shape[0] == 0 or a.shape[0] != a.shape[1]:
        raise ValueError(f"{name} must be a nonempty square matrix")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains a nonfinite entry")
    return a


def superoperator(a: ArrayLike) -> tuple[Matrix, int]:
    a = square_matrix(a, name="superoperator")
    d = isqrt(a.shape[0])
    if d * d != a.shape[0]:
        raise ValueError("superoperator dimension must equal d**2")
    return a, d


def vec(x: ArrayLike) -> Matrix:
    """Column-stack a square matrix into a one-dimensional complex array."""
    return square_matrix(x).reshape(-1, order="F")


def unvec(x: ArrayLike, d: int) -> Matrix:
    d = dimension(d)
    a = np.asarray(x, dtype=np.complex128)
    if a.ndim != 1 or a.size != d * d or not np.all(np.isfinite(a)):
        raise ValueError("vector must be finite, one-dimensional, and have length d**2")
    return a.reshape((d, d), order="F")


def matrix_unit(d: int, i: int, j: int) -> Matrix:
    d = dimension(d)
    if not (isinstance(i, (int, np.integer)) and isinstance(j, (int, np.integer))
            and 0 <= i < d and 0 <= j < d):
        raise ValueError("matrix-unit indices must be integers in range(d)")
    e = np.zeros((d, d), dtype=np.complex128)
    e[i, j] = 1
    return e


def from_action(action: Callable[[Matrix], ArrayLike], d: int) -> Matrix:
    d = dimension(d)
    columns = []
    for j in range(d):
        for i in range(d):
            image = square_matrix(action(matrix_unit(d, i, j)), name="image")
            if image.shape != (d, d):
                raise ValueError("action must map M_d to M_d")
            columns.append(vec(image))
    return np.column_stack(columns)


def apply(a: ArrayLike, x: ArrayLike) -> Matrix:
    a, d = superoperator(a)
    x = square_matrix(x)
    if x.shape != (d, d):
        raise ValueError("matrix and superoperator dimensions do not match")
    return unvec(a @ vec(x), d)


def kraus_map(operators: Iterable[ArrayLike], *, d: int | None = None) -> Matrix:
    ks = [square_matrix(k, name="Kraus operator") for k in operators]
    if d is None:
        if not ks:
            raise ValueError("supply d for an empty Kraus family")
        d = ks[0].shape[0]
    d = dimension(d)
    if any(k.shape != (d, d) for k in ks):
        raise ValueError("all Kraus operators must have shape (d,d)")
    result = np.zeros((d * d, d * d), dtype=np.complex128)
    for k in ks:
        result += np.kron(k.conj(), k)
    return result


def depolarizing(d: int) -> Matrix:
    """The UNNORMALIZED map Delta(X)=tr(X) I, not the channel Delta/d."""
    d = dimension(d)
    v = vec(np.eye(d))
    return np.outer(v, v.conj())


def congruence(s: ArrayLike) -> Matrix:
    s = square_matrix(s)
    return np.kron(s.conj(), s)


def choi(a: ArrayLike) -> Matrix:
    """J(A)=sum_ij E_ij tensor A(E_ij), using an unnormalized Omega."""
    a, d = superoperator(a)
    jmat = np.zeros_like(a)
    for i in range(d):
        for j in range(d):
            jmat[i*d:(i+1)*d, j*d:(j+1)*d] = unvec(a[:, i + j*d], d)
    return jmat


def hermiticity_residual(a: ArrayLike) -> float:
    j = choi(a)
    return float(np.linalg.norm(j - j.conj().T) / max(1.0, np.linalg.norm(j)))


def choi_min_eigenvalue(a: ArrayLike, *, tol: float = 1e-10) -> float:
    """CP diagnostic only. Nonnegative Choi spectrum does not test mere 2-positivity."""
    j = choi(a)
    if tol <= 0 or not np.isfinite(tol):
        raise ValueError("tol must be finite and positive")
    if hermiticity_residual(a) > tol:
        raise ValueError("map does not preserve Hermiticity within tolerance")
    return float(np.linalg.eigvalsh((j + j.conj().T) / 2)[0])


def trace_preservation_residual(a: ArrayLike, *, generator: bool = False) -> float:
    a, d = superoperator(a)
    ident = vec(np.eye(d))
    target = np.zeros_like(ident) if generator else ident
    return float(np.linalg.norm(a.conj().T @ ident - target) / np.linalg.norm(ident))


def spectrum_distance(x: ArrayLike, y: ArrayLike) -> float:
    """Maximum distance in a minimum-total-cost matching, not sorted-list comparison.

    This is a diagnostic matching, NOT the optimum bottleneck matching.
    """
    x, y = np.asarray(x).reshape(-1), np.asarray(y).reshape(-1)
    if x.size == 0 or x.size != y.size or not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("spectra must be finite, nonempty, and equally sized")
    costs = np.abs(x[:, None] - y[None, :])
    rows, cols = linear_sum_assignment(costs)
    return float(costs[rows, cols].max())


@dataclass(frozen=True)
class BoundReport:
    d: int
    trace: float
    trace_imaginary: float
    minimum_real: float
    maximum_real: float
    spectral_radius: float
    rhs: float
    slack: float
    tolerance: float
    holds_within_tolerance: bool
    hp_residual: float

    def as_dict(self) -> dict:
        return asdict(self)


def bound_report(a: ArrayLike, *, atol: float = 1e-10, rtol: float = 1e-9) -> BoundReport:
    """Evaluate d*m+(d**2-d)*s-Tr(A). Assumes neither positivity nor TP.

    Rejects maps that fail a Hermiticity-preservation check. A nonnegative slack
    is only a necessary spectral condition, not evidence of 2-positivity.
    """
    a, d = superoperator(a)
    if min(atol, rtol) < 0 or not np.isfinite(atol + rtol):
        raise ValueError("tolerances must be finite and nonnegative")
    hp = hermiticity_residual(a)
    if hp > max(1e-12, atol + rtol):
        raise ValueError("the theorem requires a Hermiticity-preserving map or generator")
    ev = eigvals(a)
    m, s = float(ev.real.min()), float(ev.real.max())
    tr = np.trace(a)
    rhs = d*m + (d*d - d)*s
    scale = max(1.0, abs(tr), abs(d*m) + abs((d*d - d)*s))
    tolerance = atol + rtol * scale
    if abs(tr.imag) > tolerance:
        raise ValueError("superoperator trace is not real within tolerance")
    slack = float(rhs - tr.real)
    return BoundReport(d, float(tr.real), float(tr.imag), m, s,
                       float(np.abs(ev).max()), float(rhs), slack, float(tolerance),
                       bool(slack >= -tolerance), hp)
