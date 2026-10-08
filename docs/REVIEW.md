# Correctness and reproducibility review

**Manuscript:** *Removing Trace Preservation from a Universal Spectral Bound for 2-Positive Maps*, Petri Tanninen, October 2026.  
**Materials:** the supplied nine-page pre-revision PDF and its LaTeX source (superseded; not included in this release).  
**Outcome:** all recommended corrections were adopted, and the author approved the revised text as the final manuscript in `paper/`. The review text below is retained as the record of the assessment.
**Review date:** 7 October 2026.  
**Scope:** mathematical verification of the supplied argument, targeted checks of its essential external references, source/PDF comparison, and executed support computations. 

## Verdict

**No gap was found in the proofs of Theorems 1 and 2, assuming the established trace-preserving Theorem A.** The normalization is in the correct direction, it genuinely preserves 2-positivity and produces trace preservation, the perturbation argument supplies the necessary faithful adjoint eigenmatrix, and the generator limit is valid without diagonalizability. The sharpness example and the relative-rate reformulation are also correct.

**The manuscript should nevertheless be revised before publication.** Section 8 makes an incorrect general statement about loss/gain merely shifting a spectrum. Section 9 describes an impossible closure property for a trace-preserving subclass and overstates a generic transfer principle. Neither passage is used in the proof. There is also a reproducible LaTeX cross-reference mismatch, and reference [6] can be updated to its journal publication.

The theorem being proved agrees with the conjecture in Section 3, equation (10), of the source paper. This review does not establish that the proof is the first resolution in the literature. The targeted source checks do not amount to an exhaustive priority search.

## 1. The external dependency is correctly stated

Manuscript Theorem A, equation (2.1), is the trace-preserving inequality

$$\operatorname{Tr}T\le d\,m(T)+d^2-d$$

for every 2-positive, trace-preserving map on `M_d`. This matches Theorem 1 in the original authors' accessible paper. Its Section 3 explicitly poses the homogeneous extension for conditionally 2-positive maps. The journal metadata recorded by arXiv agrees with the manuscript's reference [13]: *Linear Algebra and its Applications* **730** (2026), 262–275, DOI `10.1016/j.laa.2025.10.022`. Sources: [arXiv record](https://arxiv.org/abs/2506.02145), [original authors' PDF](https://arxiv.org/pdf/2506.02145), printed pages 4 and 11.

Using this theorem as an established input is legitimate. The manuscript is a proof by reduction, not a self-contained re-proof of that input. No extra hypothesis such as complete positivity, unitality, irreducibility or a faithful invariant state is required by the cited Theorem 1.

## 2. Detailed proof audit

### Preliminaries and Lemma 4: correct

A positive complex-linear map preserves Hermitian matrices. Representing such a map in a real basis of the Hermitian subspace gives a real matrix for its complex-linear extension, so its characteristic polynomial and superoperator trace are real. Differentiation at zero gives the same property for a generator of positive maps.

The Hilbert–Schmidt adjoint of a positive map is positive by duality of the positive semidefinite cone. The argument in (P2) applies with the stated inner-product convention. A congruence `X -> S X S^dagger` is CP, including the inverse congruence when S is invertible. This is stronger than the false general assertion that every invertible CP map has a CP inverse; the manuscript does not make that false assertion.

For the spectral lemma, `U=T*` is positive and unital, so `U(I)=I`. The order interval bound on Hermitian matrices gives `||U^n(A)||<=||A||`. Writing a general matrix as `A+iB` gives the uniform factor-two estimate. Applied to an eigenvector, `|lambda|^n<=2` forces `|lambda|<=1`. The factor two is harmless; its nth root tends to one. The conjugate-spectrum relation with the adjoint is correctly written in the LaTeX and rendered PDF. Some text extractions lose the overline; this is not a manuscript error.

The supplementary citation [14, Proposition 4.26] was checked in the [author-hosted Watrous text](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf), printed page 230. It does state the spectral-radius-one property for positive trace-preserving maps.

### Lemma 5: correct and the central reduction

Given `Phi*(H)=rH` with `H>0` and `r>0`, let `C_H(X)=H^(1/2)XH^(1/2)`. The manuscript's normalization

$$T=r^{-1}C_H\Phi C_H^{-1}$$

uses an actual superoperator similarity, not merely a congruence of the superoperator matrix. Both outer maps are CP, so 2-positivity is preserved. The trace computation is valid for arbitrary complex X because `Phi*(H)` is Hermitian. Equivalently, direct adjunction gives `T*(I)=I`.

Similarity gives `sigma(T)=sigma(Phi)/r` and `Tr(T)=Tr(Phi)/r`. Lemma 4 gives `1 in sigma(T)` and confines all its eigenvalues to the unit disk. Therefore r is not merely some positive eigenvalue: it is the spectral radius of Phi and equals `s(Phi)`. This observation justifies replacing the scaling parameter r by the largest real part in the desired inequality. Applying Theorem A and multiplying by r proves the claimed bound.

No assumption about the eigenmatrix of the original unperturbed map has been smuggled into the theorem: the next step removes it.

### Section 5 / Theorem 2: correct, including defective limits

The perturbation uses the **unnormalized** map `Delta(X)=tr(X)I`. Its Kraus sum and self-adjointness are correct. With `Phi_epsilon=Phi+epsilon Delta`, the adjoint maps every density matrix H to a matrix at least `epsilon I`. Thus the denominator in the normalized fixed-point map is at least `epsilon*d`, and the map is a continuous self-map of a nonempty compact convex set.

Brouwer's theorem therefore supplies a fixed point `H_epsilon`. Defining `r_epsilon` as the trace of its image gives a positive eigenvalue, and

$$H_\epsilon\ge(\epsilon/r_\epsilon)I>0.$$

Lemma 5 applies for each positive epsilon. Taking limits in the resulting scalar inequality is legitimate because trace and the extrema of real parts of the finite spectrum are continuous. Neither the conditioning of `H_epsilon` nor its convergence is needed. Nilpotency, reducibility, zero limiting spectral radius and Jordan blocks do not obstruct this argument.

Remark 7's statement that `r_epsilon` need not have any particular limiting behavior is best reworded: no positive lower bound needs to be assumed, but in fact the proof implies `r_epsilon=s(Phi_epsilon) -> s(Phi)`. This is a clarification, not a flaw in the limit argument.

The Kato pinpoint citation [9, Chapter II, Theorem 5.14] was not independently checked against that book edition in this review. The eigenvalue-continuity fact used here is valid, and the original source paper uses the same pinpoint citation. No mathematical conclusion here depends on an unverified differentiability claim about eigenvalues.

### Section 6 / Theorem 1: correct

The definition of conditional 2-positivity explicitly provides 2-positivity of `exp(tL)` at every nonnegative time. Spectral mapping gives its eigenvalues with algebraic multiplicities even when L is defective. For finitely many fixed eigenvalues, the Taylor remainder is uniform. In particular one may take

$$C=\tfrac12\max_j \vert{}\lambda_j\vert{}^2e^{\vert{}\lambda_j\vert{}}$$

for `0<t<=1`. The real-part errors are bounded in absolute value by `Ct^2`; the absolute-value signs are present in both the source and the rendered manuscript. Taking the minimum or maximum preserves the uniform error. Consequently

$$m(e^{tL})=1+t m(L)+O(t^2),\quad s(e^{tL})=1+t s(L)+O(t^2).$$

The trace expansion and cancellation of the `d^2` constant term are correct. Dividing by positive t and taking the limit proves the main inequality. This does not differentiate an arbitrarily chosen eigenvalue branch, so spectral crossings introduce no gap.

Remark 8 correctly proves that a 2-positive map is itself conditionally 2-positive by its exponential series. The final manuscript clarifies that a **trace-preserving semigroup generator** satisfies `L*(I)=0`, not the trace-preserving-map condition `L*(I)=I`. In that case zero is an eigenvalue and all real parts are nonpositive, yielding `s(L)=0`. For a positive trace-preserving map Phi, instead `s(Phi)=1`.

### Proposition 9 and Corollaries 10–11: correct

For amplitude damping, the matrix-unit action gives eigenvalues `0`, `-1/2`, and `-1` with multiplicities `(d-1)^2`, `2(d-1)` and `1`. The displayed eigenvectors form a basis for every `d>=2`. The Kraus formula agrees with the action on every matrix unit; completeness is exactly `K0* K0+K1* K1=I`.

For `L_c=G+c id`, real c scales the semigroup by the positive factor `exp(ct)` and shifts its spectrum by c. When `c!=0`, the semigroup is not trace-preserving for any positive time. The trace and extremal real parts yield equality exactly. In the coupled family `a m+(d^2-a)s`, subtracting `Tr(L_c)` from the right side gives `d-a`, so any `a>d` fails. The final theorem statement spells out this coupled-family interpretation of optimality.

The relative-rate algebra is exact: `sum Gamma=d^2 s-Tr(L)` and `max Gamma=s-m`. The relative-rate inequality is therefore equivalent to Theorem 1. It is a spectral statement, not a general bound on finite-time transient amplification or Jordan-block polynomial prefactors.

## 3. Required textual corrections

### R1. Section 8: loss/gain is not generally a scalar spectral shift

The supplied final paragraph claims that loss or gain of norm shifts the whole spectrum without altering relative-rate geometry. This is false for a general non-trace-preserving perturbation.

For a concrete counterexample to that sentence, take

$$K=\operatorname{diag}(0,2),\qquad L_K(X)=-\tfrac12(KX+XK).$$

Its exponential is the CP, trace-nonincreasing map `X -> diag(1,e^-t) X diag(1,e^-t)`. The spectrum is `0,-1,-1,-2`, with relative rates `0,1,1,2`. Adding this state-dependent loss to the zero generator does not uniformly shift the four zero eigenvalues. It changes the relative rates. The theorem still holds, with equality; only the explanatory sentence is wrong.

**Replacement implemented in the final manuscript:** scalar shifts `L -> L+c id` leave each relative rate unchanged; general non-TP perturbations can change the rates, but the same universal inequality remains valid.

### R2. Section 9: use an ambient class, plus a homogeneous continuous inequality

A nonempty class consisting solely of trace-preserving maps cannot be closed under arbitrary positive rescaling: `2T` is not TP when T is. Adding `epsilon Delta` also changes traces, and general normalization congruences do not preserve TP individually. Thus the concluding subclass condition cannot stand literally as written.

**Replacement implemented:** let C be an ambient class of positive maps stable under the required congruences, positive rescaling and addition of `epsilon Delta`. Let F be a continuous, positively homogeneous spectral functional, with `F(T)>=0` known for its trace-preserving members. Then the normalization gives `F(Phi_epsilon)=r_epsilon F(T_epsilon)>=0`, and continuity gives `F(Phi)>=0`. This is a precise valid transfer principle. The closure requirements concern C, not its TP subset.

The final manuscript also avoids presenting the present status of other open questions as independently verified. It says they were raised in [13] and are not addressed here.

## 4. PDF/source reproducibility and references

The original PDF is readable, has nine pages, and shows no clipped equations in the inspected renders. The source compiles successfully in two pdfLaTeX passes, with no warnings on the second pass in the local environment. However, this fresh compilation changes five reference occurrences: Corollary 11, Lemma 4 twice, Lemma 5, and Proposition 9 are labeled as Theorem 11/4/5/9. The original PDF has the correct names. The mismatch comes from shared theorem counters and cleveref type inference in the local build; it is not a changed mathematical statement.

The final source explicitly supplies label types, for example `\label[lemma]{lem:spectrum}`. The build script checks the expected reference names. Every source change from the pre-revision draft is listed in `docs/CHANGES.md`.

Reference [6] has a journal publication: *Reports on Progress in Physics* **88** (2025), 097602, DOI `10.1088/1361-6633/ae075f`. This was checked against the [authors' arXiv record](https://arxiv.org/abs/2505.24467), which records the revised version and journal reference. The final manuscript updates the bibliography.

Reference [3] supports the Doob-transform connection: the original authors normalize tilted generators using their left/adjoint eigenmatrix. See the “Quantum Doob Transform” section of [Carollo et al.](https://arxiv.org/html/1711.10951). “Map-level analogue” is more precise than the manuscript's “dual form” because the displayed generator construction uses the same orientation of congruence and adjoint eigenmatrix.

The remaining historical references were not subjected to an exhaustive bibliographic/content audit. The most load-bearing source, its hypotheses and the conjecture statement were checked directly.

## 5. Executed code validation

The local suite passed **105 tests**. The recorded reproduction used seed `20261007` and included 200 random CP maps and 200 random tilted CP-semigroup generators over d=2,...,6. It additionally checked 20 equality cases, eight explicit cases outside complete positivity, 19 normalizations, 16 epsilon-limit cases and nine generator-limit time points. All internal checks passed.

The largest recorded normalization TP residual was approximately `7.02e-15`; the largest adjoint-eigenmatrix residual was `2.11e-15`. The maximum absolute sharpness slack was `2.78e-16`. The qubit-transpose negative control returned slack approximately `-2`, as required.

The outside-CP tests are not disguised Kraus sampling. Their admissibility is proved explicitly in `MATHEMATICS.md`, including a generator whose semigroup is 2-positive but has a negative Choi expectation at sufficiently small times. The code checks those examples as well as CP examples. Symbolic checks verify the sharpness identities and exact finite-dimensional eigenbases.

These results check implementation and examples. They do **not** numerically prove the universal theorem, decide arbitrary 2-positivity, certify spectra with interval arithmetic, or establish mathematical originality. Local execution is documented; the GitHub Actions workflow is supplied but has not been run remotely.

## Recommendation

Retain the theorem and its core proof. Correct R1 and R2, repair the reference labels, and update reference [6]. (Done: the author approved the revision as the final manuscript.)