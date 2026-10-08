# Removing trace preservation from a universal spectral bound

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23225994.svg)](https://doi.org/10.5281/zenodo.23225994)

Numerical and symbolic support for **Removing Trace Preservation from a Universal Spectral Bound for 2-Positive Maps**, Petri Tanninen, Sentient Machine Corporation, October 2026.

For a superoperator on `M_d`, the paper proves, using its cited trace-preserving theorem,

$$\operatorname{Tr}A\le d\,m(A)+(d^2-d)\,s(A),$$

where `m` and `s` are the smallest and largest real parts of its eigenvalues. The claim applies to 2-positive maps and generators of 2-positive semigroups, without trace preservation.

**Status:** this repository supplies reproducibility diagnostics, not a formal proof certificate, independent peer review, or an originality claim. The mathematical proof is in the paper. The review found no gap in the main proof conditional on the cited Theorem A; its corrections to two explanatory passages, a LaTeX cross-reference issue and one bibliography entry are incorporated in the final manuscript. The paper remains the author's responsibility.

## Start here

[Correctness review](docs/REVIEW.md) · [Mathematical conventions and example proofs](docs/MATHEMATICS.md) · [Reproduction guide](docs/REPRODUCIBILITY.md) · [Revision history](docs/CHANGES.md) · [Publication checklist](docs/PUBLICATION.md)

The author-approved final manuscript (LaTeX source and PDF) is in `paper/`, with checksums in `paper/SHA256SUMS`. Superseded drafts are not included in this release; the changes from the pre-revision draft are listed in `docs/CHANGES.md`. The manuscript does not contain an invented repository URL or DOI.

## Install and run

Python 3.11 or later is required. The supplied reference run used Python 3.13.5.

```bash
python -m venv .venv
# macOS / Linux:
source .venv/bin/activate
# Windows PowerShell, instead: .venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m pytest -q
spectral-bound reproduce --output results/reproduced --seed 20261007 --samples 40
spectral-bound symbolic
```

The same command is available as `python -m spectral_bounds reproduce`. Computation after installation is offline. Do not run with `python -O`: the experiment runner deliberately refuses to run with assertion checking disabled.

For the specific numerical dependency versions used here, install `requirements-reference.txt` on Python 3.13, then install this project with `--no-deps`. See the reproduction guide for exact commands and environment limitations.

## What is checked

The reference run includes 200 random Kraus maps, 200 random tilted completely positive semigroup generators, 20 exact-sharpness numerical cases, eight explicit cases outside complete positivity, 19 normalization checks, and 16 perturbation-limit cases. Nine time steps illustrate the generator limit. `results/reference/` contains CSV data, exact-algebra output and the JSON environment/validation report. The test suite has 105 passing tests in the recorded local run.

The explicit non-CP examples matter: random Kraus sampling covers **only** completely positive maps, not the entire 2-positive cone. The repository also constructs 2-positive, non-CP maps and generators whose semigroups are not CP at small positive times, with elementary admissibility proofs in the mathematics guide. The qubit transpose is a negative control: it is positive but not 2-positive and violates the bound.

## API example

```python
from spectral_bounds import bound_report, regularized_normalization
from spectral_bounds.examples import non_cp_non_tp_map

phi = non_cp_non_tp_map(3)  # analytically 2-positive, non-CP and non-TP
report = bound_report(phi)
print(report.slack, report.holds_within_tolerance)

normalization = regularized_normalization(phi, epsilon=1e-6)
print(normalization.r)
print(normalization.trace_preservation_residual)
```

All superoperators use column-stacked vectorization. `Tr(A)` is the trace of the `d**2` by `d**2` superoperator, not `tr(A(I))`. Eigenvalues are computed with a general nonsymmetric eigensolver. The Choi matrix is used to check complete positivity and Hermiticity preservation; a positive Choi matrix is **not** a characterization of mere 2-positivity. No routine decides arbitrary 2-positivity.

The normalization routine verifies a faithful adjoint eigenmatrix, trace preservation and similarity residuals. It refuses singular or numerically ill-conditioned eigenmatrices instead of clipping negative eigenvalues. Residuals are not interval bounds; ill-conditioned eigenproblems can remain inaccurate despite small residuals.

## Repository layout

```text
src/spectral_bounds/       Reusable implementation and CLI
tests/                    Unit, symbolic and end-to-end tests
docs/                     Review, mathematics, reproducibility, publication notes
paper/                    Final manuscript: LaTeX source, PDF, SHA256SUMS
results/reference/       Recorded local numerical/symbolic results and test logs
scripts/build_paper.py     Two-pass pdfLaTeX build with warning checks
.github/workflows/        Tests and builds on GitHub Actions; not yet run remotely
CITATION.cff              Software metadata and preferred manuscript citation
```

## Attribution, license and publication

The external trace-preserving theorem and conjecture are from F. vom Ende, D. Chruściński, G. Kimura and P. Muratore-Ginanneschi, *Universal bound on the eigenvalues of 2-positive trace-preserving maps*, Linear Algebra and its Applications 730 (2026), 262–275, [arXiv:2506.02145](https://arxiv.org/abs/2506.02145), [DOI](https://doi.org/10.1016/j.laa.2025.10.022). The repository does not reproduce their paper. Other source checks appear in the review.

Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation. The software is released under the MIT License; see `LICENSE`. The manuscript in `paper/` is the author's scholarly work and is distributed under the terms the author chooses for the paper (for example, the arXiv license selected at submission). No GitHub repository, public release or archival DOI was created by preparing this package.
