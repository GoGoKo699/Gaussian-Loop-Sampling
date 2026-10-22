# Direct sampling from positive Gaussian loop weights

## Result and provenance

Let $`C=C^{\mathsf T}\ge0`$ with $`\|C\|_2<1`$, and write $`\overline N_C=\operatorname{Tr}[C^2(I-C^2)^{-1}]`$. In an ideal real-arithmetic model with exact elementary random variates, the complete Gaussian photon-count law $`P_C`$ can be sampled after $`O(m^3)`$ preprocessing using $`O(m^2)`$ stored arithmetic values. Each draw uses expected

```math
O(m+m\overline N_C)
```

arithmetic work. If every row of $`C`$ has at most $`d`$ nonzeros and neighbor lists are used for path generation, the expected draw cost is $`O(m+d\overline N_C)`$. Dense preprocessing remains cubic. Explicitly storing the output costs $`O(m)`$, and temporary path storage is proportional to the longest realized loop.

The construction is the classical **bubble decomposition of a Markov loop soup**, applied to the doubled matrix in [the model](model.md). The decomposition, logarithmic-series clusters, and negative-binomial marginal are inherited from [Lawler, Sections 5.2–5.4](https://www.math.uchicago.edu/~lawler/fall16.pdf). The theorem above is an explicit arithmetic-cost application to the full Gaussian count law, not a claim to invent loop soups or the first polynomial algorithm for nonnegative Gaussian boson sampling. [Related work](../references/RELATED_WORK.md) gives the relevant earlier claims and reductions.

This arithmetic theorem is distinct from the [finite-bit theorem](finite-precision.md). Floating-point reference code in [sampler.py](../gaussian_loops/sampler.py) demonstrates the construction but is not an implementation of certified arithmetic.

## First-return data

Order the $`2m`$ lifted vertices. At vertex $`i`$, let $`R=\{i+1,\ldots,2m-1\}`$ and use $`T=T_C`$. Set

```math
h_i=1,\qquad h_R=(I-T_{RR})^{-1}T_{Ri},\qquad
q_i=T_{ii}+T_{iR}h_R.
```

The nonnegative vector $`h_R`$ sums the weights of paths from each trailing vertex to $`i`$, with no earlier visit to $`i`$. The scalar $`q_i`$ sums the weights of first-return excursions from $`i`$ whose internal vertices lie in $`R`$. A positive inverse follows from the convergent Neumann series. The Schur complement in the positive-definite matrix $`I-T_{\{i\}\cup R}`$ gives $`0\le q_i<1`$.

If $`q_i=0`$, that vertex produces no loops. Otherwise draw a first-return excursion by choosing the first step according to

```math
\Pr(i\to j)=\frac{T_{ij}h_j}{q_i},\qquad j\in\{i\}\cup R,
```

and subsequent steps, for $`j\in R`$ with $`h_j>0`$, according to

```math
\Pr(j\to k)=\frac{T_{jk}h_k}{h_j},\qquad k\in\{i\}\cup R.
```

The identities defining $`q_i`$ and $`h`$ normalize these distributions. Zero-weight transitions are absent, so no division by a zero $`h_j`$ is required. The process terminates at its first return to $`i`$. On its reachable states the continuation kernel is a diagonal similarity transform of a submatrix of $`T_{RR}`$, which has spectral radius below one; absorption occurs almost surely.

For an excursion $`e=(i,v_1,\ldots,v_{L-1},i)`$, all $`h`$ factors telescope and its probability is

```math
\Pr(e)=q_i^{-1}\prod_{j=0}^{L-1}T_{v_jv_{j+1}}.
```

Count the starting vertex and each internal visit, but exclude the terminal return to $`i`$. This convention gives one occupation per traversed edge and prevents double-counting when excursions are concatenated.

## Poisson loops, not ungrouped excursions

Independently for each $`i`$, draw

```math
J_i\sim\operatorname{Poisson}\!\left(-\frac12\log(1-q_i)\right).
```

For each of these $`J_i`$ loops, draw a positive integer $`K`$ from the logarithmic-series law

```math
\Pr(K=k)=\frac{q_i^k}{k[-\log(1-q_i)]},\qquad k\ge1,
```

then concatenate $`K`$ independent first-return excursions at $`i`$. Each resulting cluster is **one complete Poisson loop**. Project every lifted visit to its original mode and sum the occupations over all clusters and all $`i`$.

### Why the construction has the right measure

Every loop has a unique minimum vertex $`i`$. Cutting at its visits to $`i`$ gives a cyclic collection of first-return excursions. For an **ordered** tuple of $`k`$ excursions, the cluster intensity is $`1/(2k)`$ times the product of their edge weights: the Poisson rate multiplied by the logarithmic-series probability is exactly $`q_i^k/(2k)`$, and the normalized excursion probabilities cancel $`q_i^k`$. If an unrooted loop repeats a primitive loop $`s`$ times, there are $`k/s`$ distinct ordered tuples obtained by cutting at $`i`$. Their total intensity is $`1/(2s)`$ times the edge product, which is precisely the rooted-walk measure after its $`1/L`$ factor. Thus periodic loops retain the correct multiplicity.

A determinant proof also avoids any informal orbit-counting argument. Gaussian elimination gives

```math
\det(I-T)=\prod_i(1-q_i),\qquad
\Lambda_C=-\frac12\sum_i\log(1-q_i).
```

Inserting fugacities on visited vertices replaces $`q_i`$ by the corresponding excursion generating function $`q_i(z)`$. The compound-Poisson generating function is

```math
\prod_i\exp\!\left[\frac12\sum_{k\ge1}
\frac{q_i(z)^k-q_i(1)^k}{k}\right]
=\sqrt{\frac{\det(I-T)}{\det(I-T\widehat Z)}}
=F_C(z).
```

This proves equality of the complete count laws, not only agreement of their first moments.

For count sampling alone, the total number of excursions at $`i`$ is negative binomial with shape $`1/2`$ and success probability $`1-q_i`$, and the cluster boundaries can be discarded. **Phase estimators must retain the cluster boundaries:** taking a cosine once per excursion does not produce the second moment $`e^\Gamma`$ proved in [Phase estimators](estimators.md). The implementation always retains complete loops so that the same underlying representation supports all APIs.

## Preprocessing and cost

Suppose $`G=(I-T_{RR})^{-1}`$ is already known. With $`v=T_{Ri}`$, set $`g=Gv`$, $`q=T_{ii}+v^{\mathsf T}g`$, and $`s=1-q`$. The inverse on the enlarged trailing set is

```math
\begin{pmatrix}
1/s & g^{\mathsf T}/s\\
g/s & G+gg^{\mathsf T}/s
\end{pmatrix}.
```

Processing vertices in reverse order takes $`\sum_{j=1}^{2m}O(j^2)=O(m^3)`$ operations, while storing $`h`$, $`q`$, and the current inverse takes $`O(m^2)`$ arithmetic values.

The total loop count has mean $`\Lambda_C`$. Every loop has length at least two, so $`2\Lambda_C\le\overline N_C`$. The expected number of generated visits is exactly $`\overline N_C`$. A dense categorical step takes $`O(m)`$ work; scanning only actual neighbors takes $`O(d)`$. Initializing the count vector and visiting the possible minimum vertices costs $`O(m)`$. Logarithmic-series and Poisson variates are charged unit cost in this arithmetic statement. These observations give the claimed bounds.

The reference code stores each loop as a sparse list of visits, so it does not initialize an $`m`$-entry vector for every loop. Otherwise that extra cost would invalidate the stated $`O(d\overline N_C)`$ path-generation refinement.

## Tail bounds and runtime promises

Let $`\Lambda_{>2K}`$ be the reference intensity of loops longer than $`2K`$. Since each such loop has length at least $`2K+2`$,

```math
\Lambda_{>2K}\le\frac{\overline N_C}{2K+2}.
```

The sharper norm bound is

```math
\Lambda_{>2K}
=\sum_{k>K}\frac{\operatorname{Tr}C^{2k}}{2k}
\le\frac{m\rho^{2K+2}}{(2K+2)(1-\rho^2)}.
```

Deleting these loops changes the positive reference count law in total variation by at most their occurrence probability, $`1-e^{-\Lambda_{>2K}}\le\Lambda_{>2K}`$. Signed and reweighted algorithms need the stronger error accounting in [Finite precision](finite-precision.md); a small unweighted tail probability alone is insufficient for an estimator with large possible weights.

Expected work is polynomial in $`m`$ and $`\overline N_C`$, not necessarily in the binary description length alone. A gap exponentially close to instability can create an exponentially large output. A known polynomial bound on the reference mean, or a suitable inverse-polynomial gap, is therefore a genuine resource assumption.
