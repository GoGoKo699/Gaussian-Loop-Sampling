# Verification and reference implementation

## Run

The recorded run uses Python **3.12.14** and NumPy **2.3.5**. From the repository root:

```sh
python -m pip install -r verification/requirements.txt
python verification/run_all.py --output verification/latest-results.json
```

The default statistical run uses 12,000 independent draws per case with fixed seeds. `--samples N` changes this count; at least 1,000 draws are required. The optional output file contains the same JSON printed to standard output. A successful run exits with code zero and marks each check group `"status": "passed"`. The saved [expected-results.json](expected-results.json) is a recorded diagnostic output, not a requirement for bit-for-bit agreement across numerical platforms.

## What is checked

| Module | Independent comparison | What it establishes |
|---|---|---|
| [check_identities.py](check_identities.py) | Explicit rooted-walk sums versus matrix traces and determinant moments | Factors of two, periodic-walk normalization, projected occupations, phase pairing, and the Fourier truncation formula on a real signed case |
| [check_sampling.py](check_sampling.py) | Seeded loop samples versus small hafnian probabilities and exact moments | Numerical agreement for one-mode, pair-squeezed, positive dense, gauge-equivalent, real signed, and complex-phase inputs |
| [check_topology.py](check_topology.py) | Vertex determinants, directed-edge determinants, and nonbacktracking traces | Magnetic Bass conventions, unequal magnitudes, and the explicit graph certificate calculations |

The identity checks enumerate **680 ordered rooted walks** through length eight on a two-mode support. They compare the omitted terms against an analytical tail bound; the finite sums are not treated as the complete infinite measure. The hafnian recursion is a separate exponential algorithm used only for tiny systems. Its probabilities through fourteen photons are not renormalized, so the omitted mass remains visible.

Statistical assertions use deliberately conservative sampling-scale tolerances and fixed seeds. The JSON includes observed errors, normalization, reference photon mean, and predicted and measured second moments. Passing them is evidence against implementation mistakes, not a proof that an arbitrary distribution is correct or a numerical verification of an asymptotic theorem.

The topology checks use small dense matrices for the determinant identities. The degree-six LPS example is constructed as sparse projective-group adjacency lists; its connectivity and girth checks do not allocate a dense 2,184-mode Gaussian matrix. The table of LPS total-variation bounds is calculated from the proved formula and the published girth lower bound, rather than estimated by sampling.

## Reference API

Import from the repository root:

```python
from gaussian_loops import GaussianLoopSampler, exact_probabilities
```

`GaussianLoopSampler(B, seed=None, stability_tol=1e-12, max_steps=None)` requires a finite nonempty complex symmetric matrix. It accepts symmetry differences up to absolute tolerance `1e-12` and averages them. It rejects `norm(abs(B)) >= 1 - stability_tol`. The guard is a numerical precaution; changing it does not provide a conditioning certificate.

| Method or property | Meaning |
|---|---|
| `sample_reference()` | Independent sample from the positive reference law of `C = abs(B)` |
| `iter_loops()` | A fresh reference soup, yielding complete Poisson loops as sparse visit lists with phases; exhaust the iterator to obtain a complete soup |
| `draw_cosine()` | Reference counts and a signed weight; multiplying a bounded observable by this weight gives its ideal unbiased target estimator |
| `negative_mass_estimate(K)` | Fourier approximation of negative loop mass using `K` determinant harmonics; returned truncation bounds exclude floating-point error |
| `draw_signed(negative_mass=r)` | Absolute-cosine thinning, with constant-magnitude signed weight; exact unbiasedness requires the true `r` |
| `sample_positive_part()` | Sample of the positive-part approximation, whose target error needs a proved negative-mass certificate |
| `sample_positive_part(maximum_length=L)` | Also deletes loops longer than `L`; the additional tail error must be included |
| `delta`, `gamma` | Numerical evaluations of the determinant gap and cosine second-moment exponent |
| `mean_C`, `mean_B` | Numerical reference and target photon means |

The loop construction uses a Poisson number of logarithmic-series clusters. Each cluster contains first-return excursions and is retained as one loop before evaluating a cosine. Replacing clusters by individually reweighted excursions changes the second-moment calculation.

### Example: bounded-observable estimation

```python
import numpy as np
from gaussian_loops import GaussianLoopSampler

B = np.array([[0.4, 0.2], [0.2, -0.4]], dtype=complex)
sampler = GaussianLoopSampler(B, seed=11)
values = []
for _ in range(12000):
    counts, weight = sampler.draw_cosine()
    values.append(weight * float(counts[0] > 0))
estimate = np.mean(values)
print(estimate)
```

This estimates a probability as a weighted average; the returned reference counts are not themselves samples of the complex target. The theorem in [Phase estimators](../science/estimators.md) gives confidence bounds and a median-of-means option. The code example is only an API demonstration.

For real signed matrices every loop cosine is `+1` or `-1`, so the exact identity `r = delta/2` is available analytically. For general complex matrices, a computed Fourier approximation `r_hat` used as the `draw_signed` normalizer scales the expectation by `exp(2*(r_hat-r))`. The [ratio estimator](../science/estimators.md) avoids that normalizer but requires a certified upper bound for its sample budget; its theorem is not a claim that the NumPy harmonic calculation is interval-certified.

## Numerical and runtime limits

- All matrix inverses, phases, logarithms, categorical probabilities, and random variates use floating-point routines. The implementation does not certify a total variation tolerance on arbitrary input.
- Tiny negative values of analytically nonnegative resource parameters are clipped within a stated roundoff tolerance; larger inconsistencies raise an error. Near instability, elimination and determinant evaluation can be ill conditioned.
- A final importance weight can exceed the floating-point range. `draw_cosine()` accumulates the product logarithmically, but conversion back to a float can still raise `OverflowError`.
- `max_steps` is an optional work guard. Exceeding it raises `WorkLimitExceeded` and returns no sample. Repeatedly discarding such failures would condition the law; a correct capped algorithm must budget the error and handle failures as specified in the [finite-bit theorem](../science/finite-precision.md).
- A length cutoff in `sample_positive_part()` discards loops after generation. It changes the output law but does not limit generation time.
- The code makes dense matrix preparation explicit. Neighbor lists reduce per-visit work on sparse supports, but this is not an optimized sparse-factorization package or a performance comparison against other simulators.

The [finite-bit analysis](../science/finite-precision.md) supplies a separate constructive algorithm using rational matrix powers, finite random choices, explicit tail cutoffs, and absolute-error budgets. That theorem, the ideal cubic arithmetic theorem, and the behavior of this numerical implementation are three different guarantees.
