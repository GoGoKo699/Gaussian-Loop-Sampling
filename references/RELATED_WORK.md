# Related work and the scope of the results

The results in this repository combine a direct sampler for a nonnegative Gaussian reference matrix with an instance-dependent measure of destructive loop interference. The principal phase-sensitive statements are the negative-loop estimator and positive-part sampling bound in [estimators](../science/estimators.md), with explicit graph families in [topology](../science/topology.md).

The underlying determinant expansion, Markov loop soup, and absolute-weight importance-sampling principle are established tools. Polynomial classical simulation of nonnegative graph-based GBS also has substantial prior work. The contribution should therefore be assessed through the concrete sampling cost, the phase certificate, and its proved consequences for supplied complex matrices. This document does not make a first-polynomial-simulation claim or a claim of optimal estimator variance.

## Determinants, photon counts, and loop soups

Parthasarathy and Sengupta derive multivariate photon-count generating functions for Gaussian states in [*From particle counting to Gaussian tomography*](https://arxiv.org/abs/1504.07054), *Infinite Dimensional Analysis, Quantum Probability and Related Topics* **18**, 1550023 (2015); see Eq. (3.1). Their discussion of infinite divisibility concerns particular photon-number distributions and should not be read as an algorithm for arbitrary complex full GBS output.

Le Jan's [*Markov loops, determinants and Gaussian fields*](https://arxiv.org/abs/math/0612112), arXiv:math/0612112v4 (2007), develops determinant occupation laws for Markov loop soups. Lawler's [*Topics in loop measures and the loop-erased walk*](https://www.math.uchicago.edu/~lawler/fall16.pdf) (2017), Sections 5.3–5.4, presents the negative-binomial growth and bubble constructions. These sources supply the mathematical ancestry of the direct sampler. Maintaining compound-Poisson loop clusters is essential when attaching one cosine to each complete loop; matching the occupation law alone does not preserve a phase estimator.

The use of absolute quasiprobability weights as a proposal, followed by a signed estimator, is explicit in Pashayan, Wallman, and Bartlett, [*Estimating outcome probabilities of quantum circuits using quasiprobabilities*](https://arxiv.org/html/1503.07525v2), *Physical Review Letters* **115**, 070501 (2015), Eqs. (8)–(9). In the present construction, the application-specific ingredients are the Gaussian loop representation, an efficiently samplable dominating loop measure, and determinant certificates for its negative mass. No general invention of sign sampling is claimed.

## Nonnegative graph-based GBS

### Published all-graph sampling theorem

Anand, Chen, Cryan, Freifeld, Goldberg, Guo, and Zhang, [*Simulating Gaussian boson sampling on graphs in polynomial time*](https://arxiv.org/html/2511.16558v2), arXiv:2511.16558v2 (25 June 2026), prove polynomial-time approximate sampling from the distribution on subsets of an unweighted graph with weights $`c^{2|S|}\operatorname{PM}(G[S])^2`$. Their Theorem 1.1 has no graph-density restriction. Its proof samples perfect matchings on two graph copies joined by corresponding-vertex edges. Lemma 4.2 bounds the near-perfect/perfect matching partition-function ratio. The distribution in this theorem is a collision-free subset distribution; the paper distinguishes it from the general repeated-index Gaussian formula in Section 4.

This result supersedes density restrictions relevant to Zhang et al., [*Efficient classical sampling from Gaussian boson sampling distributions on unweighted graphs*](https://doi.org/10.1038/s41467-025-64442-7), *Nature Communications* **16**, 9335 (2025). A direct loop sampler can improve the construction and cost without establishing the first polynomial algorithm for positive GBS.

### Weighted and collision-inclusive extensions: deductions from the reduction

The following extensions explain why the distinction between unweighted subsets and weighted full counts is insufficient, by itself, to establish a new polynomial-time simulability class. They are deductions presented here, rather than statements attributed verbatim to Anand et al.'s theorem.

For a nonnegative symmetric matrix $`C`$, let $`h(S)=\operatorname{haf}(C_S)`$. Use weight $`C_{ij}`$ on corresponding edges in each of two copies and weight one on each edge joining a vertex to its copy. A perfect matching then gives subset weight $`h(S)^2`$. In the near-perfect comparison, the products appearing in the two possible unmatched-vertex cases have the forms $`h(S)h(S\cup\{u,v\})`$ and $`h(S\cup\{u\})h(S\cup\{v\})`$. Applying $`ab\le(a^2+b^2)/2`$ gives the same bound

```math
Z_{n-1}<2n^2 Z_n.
```

Thus the weighted version follows through the same matching reduction, with the numerical weight and bit costs included. In the physical reference regime $`\|C\|_2<1`$, every matrix entry is below one, so the largest original edge weight is bounded.

To recover multiplicities, replace each optical mode by $`K`$ distinct clones. Give the edge between a clone of mode $`i`$ and a distinct clone of mode $`j`$ weight $`C_{ij}/K`$. This includes pairs of different clones of the same mode when $`C_{ii}>0`$. A subset with $`n_i`$ selected clones of mode $`i`$ has hafnian $`K^{-|\mathbf n|/2}\operatorname{haf}(C[\mathbf n])`$. Summing over its $`\prod_i\binom K{n_i}`$ possible clone subsets yields a projected law

```math
p_K(\mathbf n)
=\frac{p_C(\mathbf n)a_K(\mathbf n)}{\mathbb E_{P_C}a_K},
\qquad
a_K(\mathbf n)=\prod_i\frac{(K)_{n_i}}{K^{n_i}},
```

where $`(K)_n=K(K-1)\cdots(K-n+1)`$ for $`n\le K`$ and is zero for $`n>K`$. This is rejection reweighting with acceptance probability $`a_K`$, so

```math
\operatorname{TV}(P_K,P_C)
\le1-\mathbb E_{P_C}a_K.
```

For $`|\mathbf n|\le L\le K`$, the product inequality $`1-\prod_j(1-x_j)\le\sum_jx_j`$ gives $`1-a_K(\mathbf n)\le L^2/(2K)`$. Consequently,

```math
\operatorname{TV}(P_K,P_C)
\le\Pr_{P_C}(|\mathbf n|>L)+\frac{L^2}{2K}
\le\frac{\overline N_C}{L}+\frac{L^2}{2K}.
```

Choosing $`L=\max\{1,\lceil2\overline N_C/\varepsilon\rceil\}`$ and $`K\ge\max\{L,\lceil L^2/\varepsilon\rceil\}`$ gives a polynomial-size collision-inclusive reduction. The matching sampler's own error is added separately. This is a costly route, but it makes a first-polynomial claim for the positive full distribution inappropriate. The direct loop method instead has a reusable cubic dense preprocessing step and expected sample cost proportional to the reference photon number, as specified in [the sampler](../science/direct-sampler.md).

### Related positive approximations

Oh, Fefferman, Jiang, and Quesada, [*Quantum-Inspired Classical Algorithm for Graph Problems by Gaussian Boson Sampling*](https://doi.org/10.1103/PRXQuantum.5.020341), *PRX Quantum* **5**, 020341 (2024), construct a useful classical graph sampler with a different matching-weight dependence. It should not be identified with the full squared-hafnian target merely because both favor dense subgraphs.

Quesada and Arrazola, [*Exact simulation of Gaussian boson sampling in polynomial space and exponential time*](https://doi.org/10.1103/PhysRevResearch.2.023005), *Physical Review Research* **2**, 023005 (2020), already claim in their abstract a polynomial-runtime approximate sampler for encoded nonnegative matrices. Section II.D.4 uses nonnegative hafnian estimators. The claim must be distinguished from a uniform error-certified runtime: a cost containing the number of random determinant estimates does not by itself bound that number polynomially at a specified sampling error. Uvarov and Vinichenko, [*On randomized estimators of the Hafnian of a nonnegative matrix*](https://arxiv.org/html/2312.10143v2), arXiv:2312.10143v2 (2024), give exponential-relative-variance examples for standard estimators and discuss this issue. Our analysis uses a direct probability representation rather than assuming an uncontrolled hafnian-estimation variance is small.

## Recent phase-sensitive and observable-estimation results

### Lim: lossy simulation and an explicit even-cycle expansion

Youngrong Lim, [*Classical Simulation of Lossy Gaussian Boson Sampling beyond Square-Root Regime*](https://arxiv.org/html/2610.06526v1), arXiv:2610.06526v1 (5 October 2026), analyzes full PNR simulation in loss-dependent regimes. Supplement SIV.C explicitly expands the pure-core generating function as even-cycle traces. With the paper's $`W=UU^{\mathsf T}`$, $`q=\tau^2`$, and $`Y=WZW^\dagger Z`$, Eq. (S44) is the remainder

```math
\log(P_c/Q_c)
=\sum_{j\ge2}\frac{q^j}{2j}
\left(\operatorname{Tr}Y^j-N\right).
```

The subsequent analysis controls signed coefficients through convolution. Therefore, neither the even-cycle expansion of GBS nor the use of its signed remainder is new here. The present result instead uses the absolute matrix $`C=|B|`$ as a samplable envelope and quantifies the negative cosine mass for each supplied matrix, without requiring loss as the source of the approximation. Its scope is restricted by $`\|C\|_2<1`$ and by the resulting resource bounds.

### Dou et al.: Haar-average function computation

Dou, Zhan, Wu, Zhang, and Yao, [*Classical Algorithms for Function Computation in Gaussian Boson Sampling*](https://arxiv.org/html/2609.38882v1), arXiv:2609.38882v1 (30 September 2026), state a function-estimation theorem for Haar-random passive networks with bounded squeezing. Theorem 1 includes readout moment conditions and a function fixed independently of the random network; the success guarantee is over the network and algorithmic randomness. Their assumptions should not be replaced by the broader headline that every supplied GBS instance is classically simulable.

The loop theorem is instancewise and permits a bounded callable readout selected using the supplied $`B`$. Its cost can nevertheless be exponential when the negative loop mass is extensive. These are different guarantees, and the present result does not displace the Haar-average theorem.

### Lim and Oh: expectation values and operator norms

Lim and Oh, [*Classical Algorithms for Estimating Expectation Values in Linear Optical Circuits*](https://doi.org/10.1103/bl9r-2gyb), *Physical Review Letters* **137**, 120602 (15 September 2026), analyze product inputs and operator expectation values. The full theorem is available in [arXiv:2502.12882v3](https://arxiv.org/html/2502.12882v3). Their Theorem 1 has phase-space sampling/evaluation assumptions and a cost involving the squared Hilbert–Schmidt norm of the measured subsystem operator and the subsystem purity. A polynomial operator-norm factor is a sufficient efficient regime; Corollary 1 adds a sparsity premise for full sampling.

A bounded classical function of all PNR outcomes need not have polynomial Hilbert–Schmidt norm as a diagonal observable, even after a physically controlled photon cutoff. Nor does a callable representation automatically provide the phase-space sampling and evaluation routines. The loop result supplies an alternative guarantee for this readout class when its negative-loop certificate is small. It does not imply that observables already covered by Lim and Oh require loops.

Their earlier [*Approximating outcome probabilities of linear optical circuits*](https://arxiv.org/html/2211.07184v2), *npj Quantum Information* **9**, 124 (2023), Theorem 4, also gives additive approximations to individual and marginal probabilities in a bounded-squeezing/noise regime. Additive accuracy for such probabilities does not alone imply total-variation sampling of a distribution with exponentially many relevant outcomes.

### Seron et al.: distributions of linear statistics

Seron, Thomas, Araujo, Arkhipov, Oh, and Novo, [*Efficient classical algorithm for estimating linear statistics of Boson Sampling*](https://arxiv.org/html/2609.35447), arXiv:2609.35447v1 (28 September 2026), prove in Theorem 1 an additive estimator for the cumulative distribution of $`\omega\cdot\mathbf n`$, with polynomial dependence on the logarithm of the statistic's range. Section IV.3 extends the characteristic-function construction to Gaussian inputs, including displaced squeezed vacuum, using determinant formulas and photon-number truncation. Section V discusses nonlinear statistics through interacting transition amplitudes and leaves their general efficient estimation open.

Consequently, linear occupation statistics, detector bins, or their cumulative probabilities are not persuasive examples of a unique advantage for the loop estimator. The latter permits arbitrary bounded, efficiently evaluated nonlinear readouts within its certified matrix class. This does not establish a hardness result for those readouts outside that class.

## Structural and perturbative simulation methods

Oh, Lim, Fefferman, and Jiang, [*Classical Simulation of Boson Sampling Based on Graph Structure*](https://arxiv.org/abs/2110.01564), *Physical Review Letters* **128**, 190501 (2022), obtain algorithms controlled by graph width. Kolarovszki, Kaposi, Kozsik, and Zimborás, [*Simulating sparse and shallow Gaussian Boson Sampling*](https://www.iccs-meeting.org/archive/iccs2023/papers/140770201.pdf), ICCS 2023, DOI [10.1007/978-3-031-36030-5_17](https://doi.org/10.1007/978-3-031-36030-5_17), use banded/shallow structure and sufficiently long vacuum runs to factor the problem. Fixed degree of a Bargmann-matrix support graph does not imply low treewidth or a shallow one-dimensional optical implementation. The high-girth expander family must therefore be compared using these actual structural hypotheses, rather than the word “sparse.”

Li, Solomons, Bulmer, Patel, and Walmsley, [*A complexity transition in displaced Gaussian Boson sampling*](https://doi.org/10.1038/s41534-025-01062-5), *npj Quantum Information* **11**, 119 (2025), use Taylor approximation of loop hafnians in displacement-controlled regimes. Their Theorem 3 for nonnegative displacement-normalized matrices has quasipolynomial-type dependence on matrix size, row-sum norm, and precision. The normalization involves displacement, so its useful regime cannot be transferred uniformly to zero displacement. This is relevant approximate-counting methodology, rather than an existing instancewise negative-loop certificate for full undisplaced PNR output.

Sánchez-Segovia and Alhambra, [*Efficient Quantum Monte Carlo through Cluster Expansions*](https://arxiv.org/html/2610.03893v1), arXiv:2610.03893v1 (2 October 2026), sample connected cluster expansions to approximate $`\log Z`$ within a convergence regime, including complex/signed settings. Their distinction between connected-cluster sampling and physical-configuration sampling is important. For the Gaussian model here, $`Z`$ is already a determinant. Efficient computation of this normalizer does not by itself supply an estimator for an arbitrary callable function of the complete photon-count vector.

## A conservative support-polymer corollary

Small constant squeezing on bounded-degree graphs already has a generic route to full classical sampling, independently of phase or girth. This section gives a concrete comparison, including multiplicities and finite input costs. The condition derived below is deliberately sufficient; it is not a sharp GBS threshold and does not identify the largest regime of any existing algorithm.

### Published theorem being used

Chen, Galanis, Goldberg, Perkins, Stewart, and Vigoda, [*Fast algorithms at low temperatures via Markov chains*](https://arxiv.org/pdf/1901.06653), *Random Structures & Algorithms* **58**, 294–321 (2021), DOI [10.1002/rsa.20968](https://doi.org/10.1002/rsa.20968), define a subset-polymer mixing condition in Definition 1 and prove in Theorem 2 that

```math
\sum_{\gamma'\not\sim\gamma}|\gamma'|w(\gamma')
\le\vartheta|\gamma|,\qquad \vartheta<1,
```

implies an $`O(m\log(m/\varepsilon))`$ mixing-time bound for a single-polymer-update chain. A chain step still needs an implementation. Their stronger near-linear algorithm, Theorem 5, assumes both a stronger weight-decay condition and a weight oracle polynomial in polymer size. We do not invoke that runtime for Gaussian activities, whose evaluation below is exponential in the local cutoff.

### Exact support-polymer representation

Let $`G`$ be the graph of nonzero off-diagonal entries of $`B`$, with maximum degree at most a fixed $`d\ge3`$. Diagonal entries are allowed. Assume the matrix entries have rational real and imaginary parts and that a fixed bound $`\|B\|_2\le t<1`$ is supplied.

For a connected nonempty vertex set $`S`$, define

```math
a_B(S)=\sum_{\substack{n_i\ge1\ (i\in S)}}
\frac{|\operatorname{haf}(B_S[\mathbf n])|^2}{\prod_{i\in S}n_i!}.
```

If the occupied support is disconnected, its hafnian weight factors over its connected components. Thus the support law is exactly a positive polymer model: a polymer is a connected occupied set, and distinct polymers must have graph distance at least two. Its partition function is

```math
Z_B=\det(I-BB^\dagger)^{-1/2}.
```

Inclusion–exclusion gives the identity

```math
a_B(S)=\sum_{T\subseteq S}(-1)^{|S|-|T|}
\det(I-B_TB_T^\dagger)^{-1/2}.
```

This formula explains the activity but is not, on its own, a numerically stable relative-precision oracle for a very small activity. The finite construction below uses nonnegative rational sums instead.

### Activity decay and an explicit sufficient constant

Suppose $`t<1/\sqrt2`$ and set $`x=1/(\sqrt2t)>1`$. Every summand of $`a_B(S)`$ contains at least $`|S|`$ photons. Positivity and the generating-function identity imply

```math
\begin{aligned}
a_B(S)
&\le x^{-|S|}\det(I-x^2B_SB_S^\dagger)^{-1/2}\\
&\le x^{-|S|}2^{|S|/2}
=(2t)^{|S|}.
\end{aligned}
```

The connected-set count stated as Lemma 14 in Chen et al. bounds the number of size-$`k`$ connected vertex sets containing a given vertex by $`(ed)^{k-1}`$. Every set incompatible with $`\gamma`$ intersects its closed neighborhood, which has at most $`(d+1)|\gamma|`$ vertices. Therefore, whenever $`2edt<1`$,

```math
\begin{aligned}
\sum_{S\not\sim\gamma}|S|a_B(S)
&\le(d+1)|\gamma|
\sum_{k\ge1}k(ed)^{k-1}(2t)^k\\
&=(d+1)|\gamma|\frac{2t}{(1-2edt)^2}.
\end{aligned}
```

The explicit sufficient condition

```math
t\le\frac1{4e(d+1)}
```

makes the last coefficient at most $`2/e<1`$. This verifies the published mixing condition for the Gaussian support model. The parameter is the matrix-norm bound; for a $`d`$-regular graph with $`C=tA_G/d`$ and arbitrary entry phases, $`\|B\|_2\le\|C\|_2=t`$, so the same condition applies.

### Truncating both support size and multiplicities

The mixing theorem alone does not settle the cost of computing activities or sampling photon multiplicities. Both can be made polynomial here by finite enumeration at logarithmic cutoffs.

First omit polymers larger than $`K`$. Their total activity is bounded by

```math
\sum_{|S|>K}a_B(S)
\le\frac{2mt(2edt)^K}{1-2edt}.
```

Choose $`K=O(\log(m/\varepsilon))`$ so this is at most $`\varepsilon/6`$. The probability that any one polymer is present is at most its activity: its compatible remainder partition function is at most $`Z_B`$. A union bound therefore controls the probability of an omitted large polymer by the same sum.

Next truncate the photon number **inside each retained polymer** to $`L`$. For $`|S|\le K`$, the omitted activity satisfies

```math
b_{S,L}
=\sum_{\substack{n_i\ge1\ (i\in S)\\|\mathbf n|>L}}
\frac{|\operatorname{haf}(B_S[\mathbf n])|^2}{\prod_i n_i!}
\le x^{-L}2^{|S|/2}.
```

There are at most $`mK(ed)^{K-1}`$ retained connected sets. Hence it suffices to choose $`L\ge K`$ with

```math
x^{-L}mK(ed)^{K-1}2^{K/2}\le\varepsilon/6.
```

At fixed $`d`$ and $`t`$, this again requires only $`L=O(\log(m/\varepsilon))`$. This is a bound on the **sum of omitted activities**, not a relative bound divided by a potentially tiny activity.

To make the coupling precise, regard a polymer as a support together with its photon-count decoration. The original joint law assigns it its nonnegative hafnian weight. Conditioning that law on no omitted support or decoration gives exactly the finite truncated polymer model. The union bounds above show that this conditioning changes the law by at most $`\varepsilon/3`$ in total variation.

### Implementation and bit costs

For each retained support, enumerate all positive count vectors of total number at most $`L`$ and compute their weights. There are at most $`2^L`$ such vectors for each support. A subset dynamic program computes each hafnian in $`2^{O(L)}`$ arithmetic operations. With rational complex input, the weights are nonnegative rationals; their numerator and denominator bit lengths grow polynomially in the input length and $`L`$. Exact arithmetic handles zero weights and cancellations within the hafnian without assuming a lower bound on a nonzero activity.

Summing these weights gives a truncated activity no larger than $`a_B(S)`$, so the mixing condition persists. All retained polymers and their rational weights can be listed in polynomial time in $`m`$, $`1/\varepsilon`$, and the input bit length, with exponents depending on fixed $`d,t`$. The single-polymer chain of Theorem 2 can then be implemented using these finite lists and rational categorical draws. After selecting supports, draw each count decoration from its precomputed finite conditional list. Random-bit rejection routines can be given an explicit failure cap and its probability added to the total-variation budget.

This proves a polynomial full-count sampling corollary in the displayed small-norm regime. It does **not** inherit the near-linear runtime of Theorem 5. For general real-number oracles, matrix and weight approximation need an additional stability analysis; the rational-input construction suffices for the comparison made here.

### What this comparison establishes

The statement “arbitrary phases on bounded-degree graphs are easy at sufficiently small constant squeezing” is already implied by generic positive-polymer machinery. A high-girth theorem is more informative when it gives a quantitatively broader sufficient regime, an explicit efficient sampler, or a simpler certificate. The [topology result](../science/topology.md) provides such a separate analysis. Exceeding the conservative constant above establishes only that this particular sufficient condition is not being used; it does not prove that every other classical method fails.

## Claim boundaries

The defensible content is a concrete direct sampler and a phase-sensitive extension with auditable resources. The determinant Fourier series certifies the negative mass with a stated finite truncation error. Its three-determinant upper certificate, positive-part total-variation bound, and high-girth examples have distinct roles and are proved in their topic documents.

The following distinctions remain part of every complexity claim:

- Full-count sampling and arbitrary bounded readout estimation are different tasks; a signed estimator is not automatically a sampler.
- Additive estimates of a few probabilities, cumulants, or linear statistics do not imply full total-variation sampling.
- The source cost depends on the absolute-envelope mean photon number $`\overline N_C`$, which can be larger than the physical mean $`\overline N_B`$.
- Exact-real arithmetic identities do not certify a floating-point implementation. Input representation, gap bounds, truncation, and numerical errors are treated in [finite precision](../science/finite-precision.md).
- Output loss is a classical channel when placed after the modeled circuit. Other loss placements require their own model.
- Rare photon-number or collision-free postselection incurs an acceptance or conditional-estimation cost.
- None of these sufficient conditions is a hardness boundary, a general cure for a sign problem, or proof of an optimal simulation method.
