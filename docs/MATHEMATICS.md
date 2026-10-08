# Mathematical conventions and admissible examples

This file documents additional examples and implementation identities. They supplement, rather than replace, the uploaded manuscript. Numbered theorem references below refer to that manuscript.

## 1. Matrix and superoperator conventions

With column stacking, `vec(E_ij)` has its nonzero entry at `i+j*d`. Thus

$$\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}X,$$

and a Kraus map has matrix `sum(kron(K.conj(),K))`. The Hilbert–Schmidt adjoint is represented by the conjugate transpose of the superoperator. A Hermiticity-preserving superoperator is not generally a Hermitian matrix. Use general eigenvalue routines, not `eigvalsh` on the superoperator.

We use the unnormalized Choi matrix

$$J(\Phi)=\sum_{ij}E_{ij}\otimes\Phi(E_{ij}),\qquad
J(\operatorname{id})=|\Omega\rangle\langle\Omega|,\quad
|\Omega\rangle=\sum_i|ii\rangle.$$

For a single Kraus operator, `J=vec(K)vec(K)^dagger` under these index conventions. In particular `Delta(X)=tr(X)I` has `J(Delta)=I_(d^2)` and superoperator trace `d`, whereas the channel `Delta/d` has superoperator trace `1`. These identities are unit-tested independently against action on matrix units.

## 2. Slack and its transformations

Define

$$B_d(A)=d\,m(A)+(d^2-d)\,s(A)-\operatorname{Tr}A.$$

The theorem is `B_d(A)>=0` on its specified class. For `r>0`, `B_d(rA)=r B_d(A)`. Similarity preserves `B_d`, and for real `c`, `B_d(A+c id)=B_d(A)`. These facts justify the normalization in Lemma 5 and explain why scalar shifts preserve the relative rates.

For generators, the map-level inequality is applied to `exp(tL)`. With finitely many eigenvalues, the uniform Taylor error gives

$$B_d(e^{tL})=t B_d(L)+O(t^2).$$

No diagonalizability assumption enters. The code evaluates this relationship at finite positive time steps; those samples are not the proof of the limit.

## 3. Why the non-CP map examples are 2-positive

Let `R_alpha(X)=tr(X)I-alpha X`, with `0<=alpha<=1/2`. We give a direct proof of 2-positivity, not a numerical sampling argument.

For a rank-one input `|psi><psi|` in `C^2 tensor C^d`, use its Schmidt decomposition

$$|\psi\rangle=\sum_{i=1}^r\sqrt{p_i}|a_i b_i\rangle,\quad r\le2,
\qquad \rho_A=\sum_i p_i|a_i\rangle\langle a_i|.$$

For every vector `w`, Cauchy–Schwarz yields

$$|\langle\psi,w\rangle|^2
\le r\sum_i p_i|\langle a_i b_i,w\rangle|^2
\le r\langle w,(\rho_A\otimes I)w\rangle.$$

Consequently

$$(\operatorname{id}_2\otimes R_\alpha)(|\psi\rangle\langle\psi|)
=\rho_A\otimes I-\alpha|\psi\rangle\langle\psi|\ge0.$$

Decomposing an arbitrary positive semidefinite input into rank-one inputs proves 2-positivity. On the other hand,

$$J(R_\alpha)=I-\alpha|\Omega\rangle\langle\Omega|$$

has eigenvalue `1-alpha*d` on `Omega`. At `alpha=1/2` and `d>=3` this is negative, so the map is not completely positive. Indeed its `d`-fold amplification already fails on `|Omega><Omega|`.

`R_(1/2)/(d-1/2)` is trace-preserving. The map used in `non_cp_non_tp_map` is

$$\Phi=1.7\,\mathcal C_H^{-1}\,T\,\mathcal C_H,
\qquad H=\operatorname{diag}(1,2,\ldots,d).$$

Both congruences and their inverses are CP, so 2-positivity is preserved in both directions. Complete positivity is likewise preserved in both directions: if this Phi were CP, so would T be. Thus Phi remains non-CP. Its leading eigenvalue is 1.7; it cannot be positive and trace-preserving, because such maps have spectral radius 1. This also provides a known faithful adjoint eigenmatrix H.

## 4. A semigroup genuinely outside complete positivity

Let U be the cyclic shift on `C^d`, `d>=3`; `tr(U)=0`. Define

$$T(X)=\frac{\operatorname{tr}(X)I-\tfrac12 UXU^\dagger}{d-\tfrac12},
\qquad L=T+(c-1)\operatorname{id},\quad c=0.3.$$

T is 2-positive because it is `R_(1/2)` composed with unitary conjugation and scaled positively. It is trace-preserving. The exponential series gives

$$e^{tL}=e^{(c-1)t}\sum_{n\ge0}\frac{t^n}{n!}T^n,$$

which is 2-positive for every `t>=0`. Its trace is multiplied by `e^(ct)`, so it is not trace-preserving for `t>0`.

Put `u=vec(U)/sqrt(d)`. Because `tr(U)=0`, `u` is orthogonal to `Omega`. Therefore

$$\langle u,J(e^{tL})u\rangle
=t\,\frac{1-d/2}{d-1/2}+O(t^2)<0$$

for sufficiently small positive t. The zeroth-order identity-channel Choi term vanishes on u, as does the identity term in L. The semigroup is therefore not completely positive at small times. The code checks a negative Choi eigenvalue at `t=1e-4`; the analytic derivative, not the numerical test, establishes the assertion. This family exercises the actual gap between 2-positive and CP semigroups.

## 5. Completely positive generators without trace preservation

For any complex matrix A and Kraus family K_mu, set

$$L(X)=AX+XA^\dagger+\sum_\mu K_\mu X K_\mu^\dagger.$$

The drift part has CP exponential `X -> exp(tA) X exp(tA)^dagger`. The jump map is CP, and its exponential is CP by the nonnegative power series. The finite-dimensional Lie product formula then proves that `exp(tL)` is CP. This does not impose the TP condition

$$A+A^\dagger+\sum_\mu K_\mu^\dagger K_\mu=0.$$

The random tilted examples keep the original GKSL anticommutator and multiply each jump term by a positive real weight. Thus the random examples are admissible without inferring semigroup positivity from a few sampled times. Their Choi checks at `t=0.2` are additional diagnostics.

## 6. Sharpness and the necessary negative control

For the manuscript's amplitude-damping generator,

$$\sigma(G)=\{0^{[(d-1)^2]},\;(-1/2)^{[2(d-1)]},\;(-1)^{[1]}\}.$$

Shifting by c gives `Tr(L_c)=c*d^2-d`, `m=c-1`, `s=c`, and exactly zero slack. Replacing d by a in the coupled expression `a*m+(d^2-a)*s` changes `RHS-Tr` to `d-a`, negative for `a>d`.

The symbolic runner checks these expressions for a symbolic dimension d and scalar c, and verifies the stated eigenbasis on matrix units exactly for d=2,...,6. Exact checks in finitely many dimensions do not on their own prove an all-dimensional theorem; the manuscript supplies the general matrix-unit argument.

For the qubit transpose map, `Tr=2`, `m=-1`, `s=1`, hence `B_2=-2`. It is positive but not 2-positive. The test suite must detect this failure, preventing an implementation that simply labels every input valid.

## 7. Why general loss is not just a spectral shift

This is a correction to the final paragraph of manuscript Section 8, not a counterexample to its theorem. Let

$$K=\operatorname{diag}(0,2),\qquad L_K(X)=-\tfrac12(KX+XK).$$

Then `exp(tL_K)(X)=D_t X D_t` with `D_t=diag(1,exp(-t))`, a completely positive, trace-nonincreasing semigroup. Its spectrum is `0,-1,-1,-2`, and its relative rates are `0,1,1,2`. Starting from the zero generator, adding this state-dependent loss changes relative rates rather than uniformly shifting all eigenvalues. Only the special change `L -> L+c id` automatically leaves each relative rate unchanged. The universal inequality remains valid in the general case.

## 8. Numerical limits are not existence proofs

The paper uses Brouwer's theorem to establish a faithful eigenmatrix for every epsilon>0. The program instead selects the leading adjoint eigenpair and checks its residual and positivity. It does not approximate the fixed-point theorem itself. Eigenmatrix symmetrization removes only an already-checked small skew-Hermitian component; there is no eigenvalue clipping or positivity projection.

A small residual does not control eigenvalue error for all nonnormal matrices. The perturbation-limit examples include a nilpotent Kraus shift, making the defective limiting case explicit. The implementation is intended for modest d, uses dense matrices of size d^2 by d^2, and does not provide interval-certified bounds, uniform perturbation estimates, or exhaustive coverage of the 2-positive cone.
