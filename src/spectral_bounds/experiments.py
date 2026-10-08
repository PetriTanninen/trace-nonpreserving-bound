# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Deterministic experiment runner. Every reported assertion is checked."""
from __future__ import annotations

import csv
import json
import os
import platform
from pathlib import Path

import numpy as np
import scipy
import sympy
from scipy.linalg import expm

from . import __version__
from .core import (bound_report, choi_min_eigenvalue, kraus_map,
                   trace_preservation_residual)
from .examples import (amplitude_damping, amplitude_damping_kraus,
                       nilpotent_cp_map, non_cp_non_tp_map, random_cp_map,
                       random_tilted_generator, state_dependent_loss,
                       transpose_map, two_positive_non_cp_generator)
from .normalization import regularized_normalization
from .symbolic import exact_report


def _csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("cannot write an empty experiment table")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")


def reproduce(output: str | Path, *, seed: int = 20261007, samples: int = 40) -> dict:
    """Run seeded families for d=2..6 and exact algebra; fail on unmet checks.

    `samples` is the number of random CP maps AND tilted generators per dimension.
    The runner writes no network resources and performs no publication actions.
    """
    if not __debug__:
        raise RuntimeError("Run without -O: verification assertions must be enabled")
    if not isinstance(samples, int) or samples < 1:
        raise ValueError("samples must be a positive integer")
    if not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    map_rows, generator_rows, normalization_rows, sharp_rows = [], [], [], []
    for d in range(2, 7):
        for j in range(samples):
            p = random_cp_map(d, rng)
            report = bound_report(p)
            assert report.holds_within_tolerance, ("CP map", d, j, report)
            map_rows.append({"family": "random_cp", "d": d, "sample": j,
                             **report.as_dict()})
            g = random_tilted_generator(d, rng)
            report = bound_report(g)
            assert report.holds_within_tolerance, ("tilted generator", d, j, report)
            cp_min = choi_min_eigenvalue(expm(0.2*g))
            assert cp_min >= -1e-9
            generator_rows.append({"family": "random_tilted_cp_generator", "d": d,
                                   "sample": j, **report.as_dict(),
                                   "generator_tp_residual": trace_preservation_residual(g, generator=True),
                                   "exp_0_2_choi_min": cp_min})
            if j < 3:
                n = regularized_normalization(p, 1e-5)
                normalization_rows.append({"family": "random_cp", "d": d, "sample": j,
                                           **n.diagnostics()})
        for c in (-1.0, 0.0, 0.3, 2.0):
            g = amplitude_damping(d, c)
            report = bound_report(g)
            residual = float(np.linalg.norm(expm(0.37*g)-kraus_map(amplitude_damping_kraus(d, 0.37, c))))
            assert abs(report.slack) <= report.tolerance
            assert residual < 1e-10
            sharp_rows.append({"d": d, "c": c, **report.as_dict(),
                               "kraus_exponential_residual_at_0_37": residual})
    outside_cp_rows = []
    for d in range(3, 7):
        p = non_cp_non_tp_map(d)
        report = bound_report(p)
        cp_min = choi_min_eigenvalue(p)
        assert report.holds_within_tolerance and cp_min < -1e-4
        assert trace_preservation_residual(p) > 0.1
        outside_cp_rows.append({"family": "non_cp_non_tp_2positive_map", "d": d,
                                **report.as_dict(), "choi_min": cp_min,
                                "small_time_semigroup_choi_min": "not_applicable"})
        n = regularized_normalization(p, 1e-6)
        normalization_rows.append({"family": "non_cp_non_tp_2positive_map", "d": d,
                                   "sample": 0, **n.diagnostics()})
        g = two_positive_non_cp_generator(d)
        report = bound_report(g)
        exp_cp_min = choi_min_eigenvalue(expm(1e-4*g))
        assert report.holds_within_tolerance and exp_cp_min < -1e-7
        outside_cp_rows.append({"family": "genuinely_2positive_non_cp_semigroup_generator", "d": d,
                                **report.as_dict(), "choi_min": choi_min_eigenvalue(g),
                                "small_time_semigroup_choi_min": exp_cp_min})
    epsilon_rows = []
    for name, p in [("zero", np.zeros((9, 9))), ("nilpotent_shift", nilpotent_cp_map(3))]:
        exact_limit = bound_report(p)
        assert abs(exact_limit.slack) < 1e-12
        for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            n = regularized_normalization(p, eps)
            r = bound_report(n.perturbed)
            assert r.holds_within_tolerance
            epsilon_rows.append({"family": name, **n.diagnostics(),
                                 "m_perturbed": r.minimum_real, "s_perturbed": r.maximum_real,
                                 "slack_perturbed": r.slack,
                                 "slack_limit": exact_limit.slack})
    # A Hamiltonian has a non-Hermitian superoperator spectrum: tests of min/max
    # near t=0 must use real parts, not magnitudes or Hermitian eigenvalue solvers.
    from .examples import lindblad_generator, matrix_unit
    h = np.array([[0.0, 0.7], [0.7, 0.2]])
    g = lindblad_generator(h, [matrix_unit(2, 0, 1)]) + 0.17*np.eye(4)
    base = bound_report(g)
    derivative_rows = []
    for t in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
        exp_report = bound_report(expm(t*g))
        assert exp_report.holds_within_tolerance
        derivative_rows.append({"t": t, "slack_exp_over_t": exp_report.slack/t,
                                "slack_generator": base.slack,
                                "absolute_error": abs(exp_report.slack/t-base.slack)})
    assert derivative_rows[-1]["absolute_error"] < 1e-3
    transpose = bound_report(transpose_map(2))
    assert not transpose.holds_within_tolerance and abs(transpose.slack + 2) < 1e-12
    loss = bound_report(state_dependent_loss([0, 2]))
    assert abs(loss.slack) < 1e-12
    exact = exact_report()
    _csv(output/"random_cp_maps.csv", map_rows)
    _csv(output/"random_tilted_generators.csv", generator_rows)
    _csv(output/"normalization.csv", normalization_rows)
    _csv(output/"sharpness.csv", sharp_rows)
    _csv(output/"outside_complete_positivity.csv", outside_cp_rows)
    _csv(output/"epsilon_limits.csv", epsilon_rows)
    _csv(output/"generator_derivative.csv", derivative_rows)
    _json(output/"symbolic_checks.json", exact)
    _json(output/"negative_controls.json", {
        "transpose_d2": transpose.as_dict(),
        "state_dependent_loss": loss.as_dict(),
        "loss_spectrum": [0, -1, -1, -2],
        "note": "Transpose is positive but not 2-positive; state-dependent loss disproves only the broad spectral-shift sentence, not the theorem."
    })
    summary = {
        "seed": seed, "samples_per_dimension_per_random_family": samples,
        "dimensions": [2, 3, 4, 5, 6],
        "counts": {"random_cp_maps": len(map_rows), "random_tilted_cp_generators": len(generator_rows),
                   "sharpness_cases": len(sharp_rows), "outside_cp_cases": len(outside_cp_rows),
                   "normalizations": len(normalization_rows), "epsilon_limit_cases": len(epsilon_rows),
                   "generator_derivative_points": len(derivative_rows)},
        "all_internal_checks_passed": True,
        "maximum_tp_residual_in_normalization": max(x["trace_preservation_residual"] for x in normalization_rows+epsilon_rows),
        "maximum_eigenmatrix_residual": max(x["eigenmatrix_residual"] for x in normalization_rows+epsilon_rows),
        "maximum_sharpness_absolute_slack": max(abs(x["slack"]) for x in sharp_rows),
        "minimum_random_cp_slack": min(x["slack"] for x in map_rows),
        "minimum_random_tilted_generator_slack": min(x["slack"] for x in generator_rows),
        "negative_control_transpose_slack": transpose.slack,
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__, "sympy": sympy.__version__,
                        "package": __version__,
                        "thread_environment": {k: os.environ.get(k, "unset") for k in
                                               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}},
        "limitations": "Finite seeded checks, explicit non-CP families, and exact identities; not exhaustive sampling, interval arithmetic, a positivity decision procedure, or a machine-checked proof."
    }
    _json(output/"validation_summary.json", summary)
    return summary
