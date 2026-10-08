# Reproducing the support computations

## Recorded local run

The shipped run used seed `20261007`, `--samples 40`, Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0. BLAS/OpenMP thread counts were set to one. Pytest 9.0.2 executed 105 tests successfully. See `results/reference/validation_summary.json` and `pytest.txt` for the executed results, not just intended tests.

The reference run provides finite numerical evidence. It does not replace the manuscript's proof. The configured GitHub Actions matrix can be used to execute these tests remotely.

## Installation

A flexible supported install is:

```bash
python -m pip install -e ".[dev]"
```

To match the recorded scientific dependency versions, use Python 3.13 and:

```bash
python -m pip install -r requirements-reference.txt
python -m pip install --no-deps -e .
```

The reference requirements pin the scientific/test packages and their ordinary Python dependencies. They do not lock an operating-system image, BLAS build, CPU architecture, TeX distribution, or package download hashes. This provides dependency pinning rather than strict bitwise cross-platform reproducibility.

## Run commands

```bash
# POSIX shells; Windows users can set equivalent environment variables.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
python -m pytest -q
python -m spectral_bounds reproduce --output results/reproduced --seed 20261007 --samples 40
python -m spectral_bounds symbolic --max-dimension 6
```

`--samples` controls samples per dimension per random family, not total samples. Dimensions are 2 through 6. Exact outputs should agree; floating-point values can differ in their last digits across machines and dependency versions. Each experiment checks numerical tolerances before reporting success. Commands fail with nonzero exit status on an unmet assertion or failed normalization. Optimized Python execution (`-O`) is rejected so assertions cannot silently disappear.

## Results schema

`random_cp_maps.csv` and `random_tilted_generators.csv` contain superoperator trace, extremal real eigenvalues, spectral radius, right-hand side, slack, tolerance and the tolerance-qualified outcome. Tilted generators also record the generator trace-preservation residual and the Choi minimum of their exponential at t=0.2.

`outside_complete_positivity.csv` distinguishes non-CP maps from generators of genuinely non-CP, 2-positive semigroups. In this file a generator's own Choi minimum is not a test of its exponential; the latter is recorded separately at t=1e-4.

`normalization.csv` records the perturbed eigenvalue r, minimum eigenvalue and condition number of H, relative adjoint-eigenmatrix residual, trace-preservation residual, similarity residual and matched spectrum diagnostic. `epsilon_limits.csv` records these quantities as epsilon tends toward zero for the zero and nilpotent maps. Spectrum matching minimizes total absolute matching cost and reports the largest matched distance; it does not represent an optimal bottleneck distance or a certified eigenvalue error bound.

`sharpness.csv` checks equality and the Kraus/exponential identity. `generator_derivative.csv` checks the finite-time approximation to the derivative of the slack. `symbolic_checks.json` contains exact identities and finite-dimensional exact eigenbasis checks. `negative_controls.json` records the transpose violation and state-dependent-loss example.

The scalar comparison tolerance is

```text
atol + rtol * max(1, abs(Tr), abs(d*m) + abs((d*d-d)*s))
```

with default `atol=1e-10`, `rtol=1e-9`. An outcome within this tolerance is not an interval certificate. Matrix elements are never rounded before eigendecomposition, and Hermiticity-preserving superoperators are never silently replaced by Hermitian matrices.

## Manuscript build

```bash
python scripts/build_paper.py
```

Build outputs are written under `build/paper/`, leaving `paper/` untouched. pdfLaTeX and the TeX packages used by the source must already be installed. The script runs two passes, rejects undefined references/citations and overfull boxes, and checks that the expected cross-reference names (for example "By Lemma 4") appear in the PDF.

## Source integrity

`paper/SHA256SUMS` records the final manuscript bytes. Root `MANIFEST.sha256` covers the delivered repository files except itself. A ZIP checksum is supplied alongside the archive. The root manifest describes this release snapshot and must be regenerated after intentional edits; a new local reproduction directory naturally contains additional, unlisted files.