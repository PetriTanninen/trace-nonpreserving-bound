# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Reproducibility support, not computer-assisted proof, for the spectral bound."""
from .core import (BoundReport, apply, bound_report, choi, choi_min_eigenvalue,
                   congruence, depolarizing, from_action, kraus_map,
                   trace_preservation_residual, unvec, vec)
from .normalization import (Normalization, NormalizationError,
                            normalize_with_eigenmatrix, regularized_normalization)

__version__ = "0.1.0"
__author__ = "Petri Tanninen, Sentient Machine Corporation"
__copyright__ = "Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation"
__license__ = "MIT"
__all__ = ["BoundReport", "Normalization", "NormalizationError", "apply", "bound_report",
           "choi", "choi_min_eigenvalue", "congruence", "depolarizing", "from_action",
           "kraus_map", "normalize_with_eigenmatrix", "regularized_normalization",
           "trace_preservation_residual", "unvec", "vec"]
