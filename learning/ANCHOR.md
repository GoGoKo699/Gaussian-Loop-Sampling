# From Gaussian quantum information to Gaussian loop sampling

Consider an optical experiment that prepares squeezed vacuum, mixes the modes in an interferometer, and counts photons in every output mode. The classical task is to produce a complete count vector with the same distribution, or to estimate a specified bounded function of that vector. A compact description of the optical state is the starting point. The additional question is whether the photon measurement admits a classical random construction with controlled error and cost.

The project's single educational anchor is **Christian Weedbrook, Stefano Pirandola, Raúl García-Patrón, Nicolas J. Cerf, Timothy C. Ralph, Jeffrey H. Shapiro, and Seth Lloyd, [*Gaussian Quantum Information*, Reviews of Modern Physics 84, 621–669 (2012)](https://doi.org/10.1103/RevModPhys.84.621)**. The accessible [arXiv HTML](https://arxiv.org/html/1110.3234v1) and [arXiv PDF](https://arxiv.org/pdf/1110.3234v1) are versions of the same anchor. The PDF uses lettered subsections such as II.A; the HTML uses numerical labels such as II.1.

This page supplies the steps from that background to the project's input matrix, count distribution, and loop representation. The [30-milestone map](MILESTONES.md) organizes the complete scientific story. Primary research references remain in [Related work](../references/RELATED_WORK.md): they establish attribution and proof dependencies, and are not additional educational anchors.

## Relevant sections of the anchor

The useful route is concentrated in the notation discussion and Section II. The later communication, cryptography, and computation chapters are not prerequisites for the loop algorithm.

| PDF section | HTML section | Role in this project |
|---|---|---|
| I.D | I.4 | Quadrature normalization and vacuum variance |
| II.A.1–2 | II.1.1–2 | Moments, Gaussian states, and symplectic transformations |
| II.B.1–2 | II.2.1–2 | Vacuum, thermal noise, and displacement |
| II.B.3–4 | II.2.3–4 | Squeezing and optical phase |
| II.B.6 | II.2.6 | Passive mode mixing |
| II.C.1–2 | II.3.1–2 | Purity, Williamson form, and Euler/Bloch–Messiah decomposition |
| II.E.4, with II.E.1–3 for comparison | II.5.4, with II.5.1–3 | Photon counting and its distinction from Gaussian measurements |

The anchor supplies the optical and Gaussian-state foundations. It does not supply the GBS hafnian law, the Poisson-loop sampler, the phase certificates, or the nonbacktracking/Bass bounds. The derivations below and the linked science pages supply those connections, with research attribution where each ingredient enters.

## Notation dictionary

Several letters in the review have different jobs here. In particular, a covariance matrix and a pair matrix describe the same pure Gaussian state through different coordinates; they are not interchangeable matrices.

| Object | Anchor convention | Convention here |
|---|---|---|
| Number of modes | $`N`$ | $`m`$ |
| Quadrature vector | $`\hat x=(q_1,p_1,\ldots,q_m,p_m)^{\mathsf T}`$ | The same interleaved ordering, except where a grouped ordering is explicitly introduced |
| First moments | $`\bar x`$ | $`d=\langle\hat x\rangle`$; the main state model assumes $`d=0`$ |
| Covariance | Real symmetric $`2m\times2m`$ matrix $`V`$ | The same $`V`$, with vacuum covariance $`I`$ |
| Symplectic form | $`\Omega=\bigoplus_i\begin{pmatrix}0&1\\-1&0\end{pmatrix}`$ | The same form in interleaved order |
| Squeeze parameter | $`r`$ in $`S(r)`$ | $`\xi`$ in this bridge, to distinguish it from the negative loop mass $`r`$ |
| Pair matrix | Not the review's covariance blocks or beam-splitter operator | Complex symmetric $`m\times m`$ matrix $`B`$ in the creation-operator exponential |
| Magnitude reference | Not the review's covariance cross block $`C`$ | $`C_{ij}=\lvert B_{ij}\rvert`$, with the algorithmic promise $`\lVert C\rVert_2<1`$ |
| Photon output | Number-basis measurement | $`n=(n_1,\ldots,n_m)\in\mathbb N^m`$; $`\overline N_B=\mathbb E_B\sum_i n_i`$ |
| Logarithm | Section I.D allows bits or nats according to context | All loop intensities, certificates, and formulas here use the natural logarithm |

## Quadratures, covariance, and the state restriction

Use the anchor's normalization $`\hbar=2`$. The operator definitions and commutator are

```math
q_i=a_i+a_i^\dagger,\qquad
p_i=i(a_i^\dagger-a_i),\qquad
a_i=\frac{q_i+ip_i}{2},\qquad
[q_i,p_j]=2i\delta_{ij}.
```

Each vacuum quadrature has variance one. A Gaussian state is specified by its first moments and the symmetrized covariance

```math
V_{ab}=\frac12\langle
(x_a-d_a)(x_b-d_b)+(x_b-d_b)(x_a-d_a)
\rangle.
```

Physical covariance matrices satisfy $`V+i\Omega\succeq0`$. For a pure Gaussian state, all symplectic eigenvalues equal one; equivalently, $`V\Omega V=\Omega`$. These conditions refer to symplectic eigenvalues, not to the ordinary eigenvalues of $`V`$: a squeezed pure state has ordinary variances both above and below one.

The project assumes a **pure state with zero displacement**. Under these assumptions, the entire state can be written using one symmetric pair matrix:

```math
|\psi_B\rangle
=\det(I-B^\dagger B)^{1/4}
\exp\!\left(\frac12\sum_{i,j}B_{ij}a_i^\dagger a_j^\dagger\right)|0\rangle,
\qquad B=B^{\mathsf T},\quad\|B\|_2<1.
```

The exponential creates pairs. Its normalizability condition is that every singular value of $`B`$ lies below one. Mixed thermal noise requires more information, and nonzero displacement adds a linear creation term. Neither modification is silently included in this formula.

There is a useful physical check on every translation of conventions. The commutator gives

```math
a_i^\dagger a_i=\frac{q_i^2+p_i^2-2}{4},
\qquad
\overline N_B=\frac{\mathrm{Tr}V-2m}{4}
\quad\text{when }d=0.
```

Thus the excess total quadrature variance determines the mean photon number. With nonzero first moments the numerator would also contain $`\|d\|^2`$. The factor of four is specific to the chosen quadrature normalization.

## Translating covariance into the pair matrix

The canonical statement and proof are in [Covariance and pairing conventions](../science/model.md#covariance-and-pairing-conventions). The following calculation explains the translation operationally.

First permute the quadratures into grouped order, $`x_g=(q_1,\ldots,q_m,p_1,\ldots,p_m)^{\mathsf T}`$, and form the blocks

```math
V_g=\begin{pmatrix}V_{qq}&V_{qp}\\V_{pq}&V_{pp}\end{pmatrix},
\qquad V_{pq}=V_{qp}^{\mathsf T}.
```

Define two complex moment matrices with their indices fixed as

```math
K_{ij}=\langle a_i a_j^\dagger\rangle,
\qquad M_{ij}=\langle a_i a_j\rangle.
```

Here $`K`$ includes the vacuum contribution; it is not the matrix with entries $`\langle a_i^\dagger a_j\rangle`$. Expanding the annihilation operators gives

```math
K=\frac{V_{qq}+V_{pp}+2I+i(V_{pq}-V_{qp})}{4},
\qquad
M=\frac{V_{qq}-V_{pp}+i(V_{qp}+V_{pq})}{4}.
```

The $`2I`$ in $`K`$ comes from commuting a position and a momentum into symmetric order. Omitting it would already give the wrong answer for the vacuum, where $`K=I`$ and $`M=0`$.

For a pure pair state, the canonical nullifier proof gives $`M=KB`$. The operational conversion is therefore to solve the matrix equation $`KB=M`$:

```math
B=K^{-1}M.
```

The order matters for multimode matrices. The inverse relation is $`K=(I-BB^\dagger)^{-1}`$, $`M=KB`$; the canonical page gives all four real covariance blocks explicitly. This also makes $`\overline N_B=\mathrm{Tr}(K-I)`$ agree with the covariance trace formula.

Purity is essential to the argument: the annihilation identities are identities for the pure state vector. For example, a thermal mode with nonzero mean occupation $`\bar n_{\rm th}`$ has $`V=(2\bar n_{\rm th}+1)I`$, $`M=0`$, and $`K=(\bar n_{\rm th}+1)I`$. Applying $`K^{-1}M`$ would give $`B=0`$, whose pair state is vacuum. That computation would discard the thermal noise rather than represent the original state.

## Squeezing sign and passive mixing

The sign of a pair coefficient depends on the squeeze-operator convention. With the anchor's Eq. (31), written using $`\xi`$ here,

```math
S(\xi)=\exp\!\left[\frac{\xi}{2}(a^2-a^{\dagger2})\right],
\qquad
S^\dagger(\xi)aS(\xi)=a\cosh\xi-a^\dagger\sinh\xi.
```

For $`S(\xi)|0\rangle`$, the annihilation identity is $`(a\cosh\xi+a^\dagger\sinh\xi)|\psi\rangle=0`$. Comparing it with $`(a-ba^\dagger)|\psi_b\rangle=0`$ gives

```math
b=-\tanh\xi,\qquad
V=\mathrm{diag}(e^{-2\xi},e^{2\xi}).
```

The sign can also be checked at first order: $`S(\xi)|0\rangle=|0\rangle-(\xi/2)a^{\dagger2}|0\rangle+O(\xi^2)`$. A positive coefficient $`b=\tanh\xi`$ corresponds instead to $`S(-\xi)|0\rangle`$, with the two squeezed axes exchanged.

**Version convention.** The arXiv v1 Fock expansion printed in Eq. (33) uses positive powers of $`\tanh\xi`$. The operator definition in Eq. (31) and covariance transformation in Eq. (32) fix the negative sign derived above. The formulas here follow those operator and covariance conventions; using the positive coefficient entails the stated change of squeeze axis. The two states have the same photon-count probabilities.

For many modes, write the Takagi factorization

```math
B=U\mathrm{diag}(s_1,\ldots,s_m)U^{\mathsf T},
\qquad s_j=\tanh\xi_j\in[0,1).
```

To fix the passive convention, let $`\widehat U^\dagger a\widehat U=Ua`$. Applying $`\widehat U`$ to a pair state sends $`B`$ to $`UBU^{\mathsf T}`$. The Takagi form therefore describes independent positive pair coefficients, prepared by $`S(-\xi_j)`$, followed by this passive transformation. The minus sign of the opposite squeeze convention can be absorbed in column phases of $`U`$.

This is the pair-matrix expression of the anchor's Euler/Bloch–Messiah preparation. The first passive transformation in that decomposition leaves vacuum invariant; the squeezers and the final interferometer determine the state. For photon counts, the resulting mean is

```math
\overline N_B=\sum_j\frac{s_j^2}{1-s_j^2}
=\sum_j\sinh^2\xi_j.
```

## From pair creation to the photon generating function

The following steps connect the Gaussian background to the count law used in [the model](../science/model.md). The GBS hafnian connection is attributed to [Hamilton et al.](https://doi.org/10.1103/PhysRevLett.119.170501), separately from the educational anchor.

Expand the creation exponential. To obtain an amplitude with $`2k`$ photons, choose $`k`$ pair-creation terms. Their products enumerate pairings of the $`2k`$ labeled photon positions. The expansion factor $`1/(2^k k!)`$ removes the ordering of the pairs and the ordering within each pair. The remaining sum is the hafnian.

If $`B[n]`$ repeats row and column $`i`$ exactly $`n_i`$ times, the number-basis normalization gives

```math
\langle n|\psi_B\rangle
=\det(I-B^\dagger B)^{1/4}
\frac{\mathrm{haf}(B[n])}{\sqrt{\prod_i n_i!}}.
```

Odd total occupation has zero amplitude. Squaring this expression gives the probability law in the model. Repeated photons in a mode are allowed, and the factorials are required for those collisions. The normalizing determinant follows mode by mode after Takagi factorization: an unnormalized single-mode pair state has squared norm $`(1-s^2)^{-1/2}`$.

Introduce a marker $`z_i`$ for each mode and put $`Z=\mathrm{diag}(z_1,\ldots,z_m)`$. Multiplying a Fock amplitude by $`\prod_i z_i^{n_i/2}`$ multiplies each pair coefficient $`B_{ij}`$ by $`\sqrt{z_i z_j}`$. The squared norm of the resulting state is the probability generating function. Taking the ratio of squared norms for the corresponding unnormalized pair states, and using cyclic invariance of a determinant, gives

```math
F_B(z)=\mathbb E_B\prod_i z_i^{n_i}
=\sqrt{\frac{\det(I-B^\dagger B)}
{\det(I-B^\dagger ZBZ)}}.
```

Initially take $`0\le z_i\le1`$ so all operations are directly bounded. The normalized branch has $`F_B(1,\ldots,1)=1`$. This identity encodes every collision pattern and every even photon total at once. It does not require computing the hafnian of a selected output pattern.

## Why the logarithm produces a classical loop process

First take a real, entrywise nonnegative pair matrix $`C=C^{\mathsf T}`$ with $`\|C\|_2<1`$. The trace-log expansion gives

```math
\log F_C(z)
=\sum_{k\ge1}\frac{
\mathrm{Tr}(CZ)^{2k}-\mathrm{Tr}C^{2k}}{2k}.
```

A term in a trace is a closed walk of length $`2k`$ through the modes. Its weight is the product of the traversed matrix entries. Its markers record one occupation for each visited vertex, excluding the duplicated terminal vertex. All such weights are nonnegative. Walks may backtrack, revisit a mode, and repeat a shorter loop; restricting the sum to simple cycles would change the distribution.

The factor $`1/(2k)`$ keeps the multiplicity associated with the walk's root. The [model](../science/model.md#lifted-walks) makes this precise with a doubled graph and a loop measure $`\mu_C`$. Its total mass is

```math
\Lambda_C=\int d\mu_C
=-\frac12\log\det(I-C^2)<\infty.
```

Let $`n(w)`$ be the occupation vector of one complete loop. The same generating function becomes

```math
F_C(z)=\exp\!\left[\int\bigl(z^{n(w)}-1\bigr)\,d\mu_C(w)\right].
```

This is the generating function of a compound Poisson distribution: draw a Poisson number of loops of mean $`\Lambda_C`$, draw their shapes with probabilities proportional to their loop weights, and add their occupation vectors. The determinant identity has now become a sampling prescription because its logarithmic coefficients form a positive measure.

There are infinitely many possible loop shapes, but only finitely many are drawn almost surely. The [direct sampler](../science/direct-sampler.md) generates them through first-return paths and logarithmic-series clusters, without enumerating the infinite set. That decomposition comes from Markov-loop theory, as attributed in the science pages. These auxiliary loops describe a classical probability construction for the photon outputs; they are not physical paths taken by individual photons.

## A complete one-mode example

Choose the pair coefficient $`b=1/2`$. It is prepared by $`S(-\xi)|0\rangle`$ with $`\xi=\mathrm{arctanh}(1/2)`$. The covariance and complex moments are

```math
V=\begin{pmatrix}3&0\\0&1/3\end{pmatrix},\qquad
K=\frac43,\qquad M=\frac23,\qquad K^{-1}M=\frac12.
```

The mode is squeezed: one quadrature variance is $`1/3`$. Both descriptions give the same mean photon number,

```math
\overline N=\frac{3+1/3-2}{4}
=\frac{b^2}{1-b^2}
=\frac13.
```

Expanding the state in the number basis yields

```math
|\psi_b\rangle=(1-b^2)^{1/4}
\sum_{k\ge0}\frac{\sqrt{(2k)!}}{2^k k!}\,b^k|2k\rangle,
```

```math
P_b(2k)=\sqrt{1-b^2}\binom{2k}{k}\left(\frac{b^2}{4}\right)^k,
\qquad P_b(2k+1)=0.
```

For this value of $`b`$, the first probabilities are

| Photon count | Probability | Approximate value |
|---|---|---|
| $`0`$ | $`\sqrt3/2`$ | $`0.8660254`$ |
| $`2`$ | $`\sqrt3/16`$ | $`0.1082532`$ |
| $`4`$ | $`3\sqrt3/256`$ | $`0.0202975`$ |

The count generating function and its logarithm are

```math
F_b(z)=\sqrt{\frac{1-b^2}{1-b^2z^2}},\qquad
\log F_b(z)=\sum_{k\ge1}\frac{b^{2k}}{2k}(z^{2k}-1).
```

Thus an exact ideal sampler can draw independent Poisson variables $`L_k`$ with means $`\lambda_k=4^{-k}/(2k)`$ and return $`N=\sum_{k\ge1}2kL_k`$. An implementable equivalent avoids iterating over all $`k`$:

1. Draw $`J\sim\mathrm{Poisson}(\Lambda)`$, where $`\Lambda=\tfrac12\log(4/3)\simeq0.143841`$.
2. For each of the $`J`$ loops, draw $`K\ge1`$ with logarithmic-series probability $`\Pr(K=k)=(1/4)^k/[k\log(4/3)]`$.
3. Add $`2K`$ photons for that loop and return the total.

The vacuum occurs when no loop is drawn, giving $`e^{-\Lambda}=\sqrt3/2`$. Two photons require one length-two loop, giving $`e^{-\Lambda}\lambda_1=\sqrt3/16`$. Four photons have two possible loop decompositions: one length-four loop, or two length-two loops. Their combined probability is

```math
e^{-\Lambda}\left(\lambda_2+\frac{\lambda_1^2}{2}\right)
=\frac{\sqrt3}{2}\left(\frac1{64}+\frac1{128}\right)
=\frac{3\sqrt3}{256}.
```

This calculation exposes why independent pairs with a single Poisson rate would be insufficient: longer clusters contribute to the same output alongside collections of shorter clusters. The expected generated occupation is also exact,

```math
\mathbb E N=\sum_{k\ge1}2k\lambda_k
=\sum_{k\ge1}(1/4)^k=\frac13.
```

Replacing $`b`$ by $`e^{i\theta}/2`$ rotates the covariance ellipse but leaves every count probability unchanged. In one mode each closed-walk phase cancels against its conjugate, so all three phase resources below vanish. Multimode loops are needed to exhibit nontrivial phase interference. The positive count sampler in this example coexists with quadrature squeezing; it does not assert that the state is a classical mixture of coherent light.

## Complex phases and the algorithmic resource

For a general complex symmetric $`B`$, the construction uses $`C=|B|`$ entrywise. The additional promise $`\|C\|_2<1`$ is stronger than the physical requirement $`\|B\|_2<1`$. It provides an absolutely convergent positive reference ensemble whose preparation and output cost must be counted.

The doubled graph alternates edges of $`B`$ and $`B^\dagger`$. A loop therefore acquires an alternating phase. For example, a length-four loop on modes $`i,j,k,\ell`$ has phase

```math
\phi=\arg B_{ij}-\arg B_{jk}+\arg B_{k\ell}-\arg B_{\ell i}
\quad\text{modulo }2\pi.
```

Reversing the loop conjugates the phase without changing its occupation. Pairing the two orientations replaces the complex factor by $`c(w)=\cos\phi(w)`$. The complete count law obeys

```math
\log F_B(z)=\int c(w)\bigl(z^{n(w)}-1\bigr)\,d\mu_C(w).
```

Negative cosine values prevent a direct interpretation of this expression as a positive Poisson intensity. The project quantifies their effect through

```math
\Delta=\int(1-c)\,d\mu_C,\qquad
\Gamma=\int(1-c)^2\,d\mu_C,\qquad
r=\int(-c)_+\,d\mu_C.
```

Here $`r`$ is the negative loop mass, not a squeezing parameter. [Phase estimators](../science/estimators.md) gives their precise operational statements:

- Keeping the reference loops and using $`W=e^\Delta\prod_{w\text{ drawn}}c(w)`$, with one cosine per **complete Poisson loop**, gives an unbiased bounded-observable estimator with weight second moment $`e^\Gamma`$.
- Thinning the loop intensity by $`|c|`$ and using $`W_{\rm abs}=e^{2r}S`$, where $`S`$ is the product of the retained loop signs, gives second moment $`e^{4r}`$. A ratio estimator can instead use a certified upper bound on $`r`$.
- Keeping only intensity $`\max(c,0)\,d\mu_C`$ gives an ordinary positive sampler whose total variation error is at most $`\min\{1,(e^{2r}-1)/2\}`$.

The distinction between sampling and estimation remains essential. A signed weighted draw estimates an expectation; it is not automatically a sample from the target law. The bounded observable must also be efficiently evaluable on a count vector. Certificate computation, such as the determinant/Fourier upper bound on $`r`$, is part of the cost.

Individual entry phases do not determine these resources. A local phase change $`B\mapsto DBD`$, with diagonal unitary $`D`$, cancels around every alternating closed walk and leaves the photon law unchanged. Backtracking also cancels. This explains the connection to [Topology](../science/topology.md): on suitable graphs with long cycles, phase-sensitive loops are long, and a quantitative bound on their combined weight can certify classical approximation even with arbitrary edge phases. Girth and sparsity alone do not establish that bound; the spectral and nonbacktracking estimates in that page supply the required control.

## Completing the connection to the research results

The anchor explains why the inputs are Gaussian optical states and how their covariance, squeezing, and passive transformations fit together. Photon measurement then introduces the hafnian combinatorics. The loop method reorganizes the *whole count generating function* so that positivity, or a controlled negative part, becomes the relevant sufficient condition for a classical algorithm.

The remaining proof dependencies have direct locations:

| Question | Authoritative project page |
|---|---|
| What is the exact state, count law, reference promise, and loop measure? | [Model](../science/model.md) |
| How are complete loops sampled without enumerating all walks? | [Direct sampler](../science/direct-sampler.md) |
| How are phases, signs, observables, and certificates handled? | [Phase estimators](../science/estimators.md) |
| What turns exact arithmetic identities into a finite-bit guarantee? | [Finite precision](../science/finite-precision.md) |
| How do graph structure and arbitrary optical phases coexist with approximation? | [Topology](../science/topology.md) |
| Which physical and computational conclusions follow? | [Interpretation](../science/interpretation.md) |
| Which ingredients and overlapping results come from earlier research? | [Related work](../references/RELATED_WORK.md) |

For the positive reference law, the direct construction uses cubic dense preprocessing and expected draw work controlled by $`m`$ and $`\overline N_C`$. This is an ideal-arithmetic statement. The finite-bit theorem separately accounts for input precision, stability, tail cutoffs, and approximate random choices; the included floating-point implementation is a numerical realization rather than certified arithmetic. A small target mean $`\overline N_B`$ does not replace the reference-mean promise.

Output loss and efficiently sampled classical processing can be applied to the generated counts with their error and cost accounted for. Conditioning on a rare photon sector or on collision-free outcomes requires an additional acceptance-probability analysis. The unconditional full-count result does not remove that cost. Likewise, a compact covariance description, a small loop certificate, and support-graph treewidth are different statements: none alone establishes quantum hardness, entanglement scaling, or optimality of this classical algorithm.
