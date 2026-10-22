# Physical meaning and computational boundaries

## What the classical variables represent

A Gaussian state's photon amplitudes coherently sum over pairings. The loop construction reorganizes the squared amplitudes after summing over output occupations. Its classical random variables are closed walks with repeated vertices and backtracking, not candidate perfect matchings of a selected output pattern. This is why a sample can be generated without evaluating a hafnian.

The determinant is useful at the level of the *whole generating function*. A determinant identity alone is not yet a sampler: its logarithmic coefficients must form a positive measure, or their signs must be handled with a quantitative error or variance bound. The [direct construction](direct-sampler.md) supplies this missing probabilistic step for nonnegative matrices. The [phase analysis](estimators.md) identifies what survives when the edge weights are complex.

The resulting algorithm is entirely classical. “Quantum inspired” describes the input distribution and the quantum-optical motivation; it does not assert a classical speedup over the best unrelated algorithm, and it requires no quantum oracle.

## Interference is organized by closed-walk phases

Local changes of optical phase act as $`B\mapsto DBD`$ and do not change photon counts. Accordingly, an individual entry phase is not an invariant measure of difficulty. The lifted graph turns the physically relevant alternating phase products into ordinary loop holonomies. Backtracking traverses an edge and its conjugate, so it contributes no phase. On a locally tree-like support, phase-sensitive loops must be long.

This observation connects the algorithmic results:

- Positive or gauge-equivalent matrices have an exact positive loop representation.
- General stable-magnitude matrices have signed observable estimators. Their efficiency depends on the reference output size and a variance or negative-mass certificate.
- Small negative loop mass gives a positive approximation to the full photon law.
- Large girth and controlled nonbacktracking growth can make the *total* phase-sensitive loop weight small even when every individual edge phase is arbitrary.

The last mechanism requires no small-angle promise. Conversely, a small-angle expansion gives a quartic variance improvement after pairing orientations, but it should not be described as an exclusive regime for expectation estimation: the positive short-loop construction can also provide full-distribution sampling in part of that regime. The precise bounds and examples are in [Phase estimators](estimators.md).

## An extensive-photon example

[Topology](topology.md) applies the loop certificate to a specified degree-six Lubotzky–Phillips–Sarnak family. These are connected regular expanders with logarithmic girth and linear support-treewidth. Their construction and spectral properties come from [Lubotzky, Phillips, and Sarnak, *Ramanujan graphs* (1988)](https://doi.org/10.1007/BF02126799); the photon-sampling conclusion is derived in this repository.

For $`C=(t/6)A_G`$ and arbitrary symmetric phases on its edges, $`B=C\circ e^{i\theta}`$ obeys the stability promise for any fixed $`t<1`$. The total mean photon number is extensive in the number of modes. In the useful illustrative case $`t=1/2`$,

```math
\frac{m}{24}\le\overline N_B,\overline N_C\le\frac{m}{18}.
```

The Ramanujan nonbacktracking bound proves that sampling $`P_C`$ approximates $`P_B`$ with a specified error tending to zero along this family. The finite-instance certificate, its asymptotic rate, and the larger admissible coupling interval are stated in the topology page. An accuracy that shrinks with family size is not a claim to meet every requested error at a fixed matrix without changing the algorithm or parameters.

This example avoids an explanation based solely on a bounded photon count or bounded support-treewidth. It does not prove a volume law for entanglement, a quantum hardness result for nearby instances, or optimality of the classical algorithm. Support-treewidth is a property of a graph representation; Gaussian entanglement depends on the covariance matrix and on the bipartition. The two should not be identified.

## Gaussian does not mean an easy photon distribution

Gaussian states admit compact first- and second-moment descriptions, and Gaussian quadrature measurements are particularly tractable. Photon-number measurement is a different, non-Gaussian measurement. The optical background and this measurement distinction are explained in [Weedbrook et al., Sections II.A–II.E](https://doi.org/10.1103/RevModPhys.84.621). The hafnian formula captures the additional combinatorics of photon detection.

The positive loop representation is a condition on a count-generating function. It is not a claim that the quantum state is a classical mixture of coherent states or that it lacks squeezing or entanglement. Even a single squeezed mode has a nontrivial positive count-loop representation.

## What each resource actually pays for

| Resource | Operational role | Limitation |
|---|---|---|
| Number of modes $`m`$ | Matrix preparation, categorical scans, output representation | Cubic dense preprocessing may dominate practical runtime |
| Reference mean $`\overline N_C`$ | Expected number of generated visits | May be much larger than the target mean $`\overline N_B`$ |
| Gap below $`\|C\|_2=1`$ or a mean bound | Tail cutoffs and finite precision | An exponentially unstable reference need not give an efficient algorithm |
| $`\Gamma`$ | Exact cosine-estimator second moment through $`e^\Gamma`$ | Sufficient certificate for this estimator, not a hardness criterion |
| Negative loop mass $`r`$ | Absolute-cosine sign cost and positive-part sampling error | Its efficient approximation or an upper bound must be charged |
| Observable evaluation cost | Computes the requested bounded function on each draw | A formula specifying a hard-to-evaluate function is not free access |
| Requested error and failure probability | Repetition count, rounding, truncation | Additive accuracy does not imply relative accuracy for rare events |

The [finite-bit theorem](finite-precision.md) makes these promises explicit. It supplies a constructive polynomial implementation with exact rational bridge probabilities and controlled approximate random choices. Its matrix-power preprocessing is different from the cubic ideal-arithmetic excursion algorithm. The included NumPy implementation is a reproducible numerical realization of the latter; it is not certified interval software or a large-scale performance claim.

## Detection, conditioning, and extensions

Independent detector loss and any efficiently sampled classical processing of the observed counts can be applied after the loop algorithm. Total variation cannot increase under this common stochastic map. Bounded observables of the processed counts can likewise be evaluated by including that randomness in the observable estimator. This does not automatically cover a circuit that changes later quantum operations in response to earlier measurements.

Conditioning on a rare event requires an additional cost and accuracy analysis. For example, postselecting on a photon sector or on collision-free outcomes can make an otherwise inexpensive unconditional sampler expensive. The event probability and the denominator in a conditional expectation must be included; the unconditional theorems do not hide this cost. [Phase estimators](estimators.md) gives the quantitative statements.

The established contribution concerns zero-displacement pure Gaussian inputs with a stable entrywise-magnitude reference. General mixed inputs, arbitrary displacement, or an unstable magnitude reference need a separate representation and analysis. Failure of these sufficient conditions is not evidence of classical hardness. The [related-work comparison](../references/RELATED_WORK.md) records the relevant earlier algorithms, overlapping observables, and inherited loop identities so that these boundaries remain visible.
