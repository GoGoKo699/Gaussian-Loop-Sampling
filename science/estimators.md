# Phase-sensitive estimation and sampling

The central resource in this document is the **negative mass of a real loop measure**. It gives two consequences for the same supplied Gaussian matrix: an estimator for arbitrary bounded, computable functions of the complete photon-count vector, and a positive sampling distribution with a certified total-variation error. Determinants provide computable upper bounds and additive approximations to this resource.

The target is the full, pure, undisplaced photon-number-resolving distribution. Collisions and fluctuating photon number are included. The definitions of the Gaussian state and its generating function are in [the model](model.md); the source loop sampler is in [the direct sampler](direct-sampler.md). Arithmetic-operation counts below must be distinguished from the certified implementation conditions in [finite precision](finite-precision.md).

## 1. Input, output, and cost parameters

Let $`B=B^{\mathsf T}\in\mathbb C^{m\times m}`$, set $`C=|B|`$ entrywise, and assume

```math
\rho=\|C\|_2<1.
```

This condition is stronger than physical normalizability, $`\|B\|_2<1`$. It guarantees a finite positive reference loop measure. An arbitrary physical Gaussian matrix need not satisfy it.

For any symmetric matrix $`X`$ with $`\|X\|_2<1`$, define

```math
\Lambda_X=-\frac12\log\det(I-X X^\dagger),
\qquad
p_X(\mathbf n)=e^{-\Lambda_X}
\frac{|\operatorname{haf}(X[\mathbf n])|^2}{\prod_i n_i!}.
```

The empty hafnian is one, and an odd-order hafnian is zero. Replacing $`XX^\dagger`$ by $`X^\dagger X`$ leaves the determinant unchanged. The reference photon number is

```math
\overline N_C
=\operatorname{Tr}\!\left[C^2(I-C^2)^{-1}\right].
```

The algorithms pay for $`\overline N_C`$, which can exceed the physical mean $`\overline N_B`$ substantially. In particular,

```math
\overline N_B\le\overline N_C,
\qquad
2\Lambda_C\le\overline N_C,
\qquad
\rho^2\le\frac{\overline N_C}{1+\overline N_C}.
```

The first inequality follows by expanding the traces in powers and bounding each alternating product by its absolute value. The second follows term by term from $`1/k\le1`$ in the logarithmic series. The third retains the largest eigenvalue's contribution to $`\overline N_C`$.

The readout is a callable function $`F:\mathbb N^m\to[-1,1]`$. It may depend on the supplied matrix $`B`$ and may be nonlinear or involve all modes. Its evaluation cost is part of the algorithm. The theorem does not supply an efficient procedure for a readout that is itself computationally hard. Complex-valued bounded readouts can be treated through their real and imaginary parts.

## 2. The signed loop representation

An even rooted walk is an ordered tuple

```math
w=(i_1,\ldots,i_{2k}),\qquad i_{2k+1}=i_1,\qquad k\ge1.
```

It contributes occupation vector $`v_i(w)=|\{j:i_j=i\}|`$. Define its positive intensity and alternating phase by

```math
\mu_C(w)=\frac1{2k}\prod_{j=1}^{2k}C_{i_j i_{j+1}},
\qquad
\phi(w)=\sum_{j=1}^{2k}(-1)^{j-1}\arg B_{i_j i_{j+1}},
\qquad c(w)=\cos\phi(w).
```

Walks containing a zero-weight edge have zero intensity and need no phase convention. Summation over these ordered tuples, with the displayed $`1/(2k)`$ factor, is one precise rooted-walk convention; it must not be replaced by a sum over distinct cycles with the same factor. Equivalent unrooted or doubled-graph conventions are valid when their multiplicities are included.

The logarithmic determinant expansion gives

```math
\int d\mu_C=\Lambda_C,
\qquad
\int |v(w)|\,d\mu_C(w)=\overline N_C,
\qquad
\int c(w)\,d\mu_C(w)=\Lambda_B.
```

Writing $`z^{v(w)}=\prod_i z_i^{v_i(w)}`$, the photon generating function is

```math
G_B(z)
=\sum_{\mathbf n}p_B(\mathbf n)z^{\mathbf n}
=\exp\left\{\int c(w)\bigl(z^{v(w)}-1\bigr)d\mu_C(w)\right\}.
```

To see the cosine explicitly, expand $`-\tfrac12\log\det(I-BZB^\dagger Z)`$, where $`Z=\operatorname{diag}(z)`$. Each term is an alternating closed walk. Rotating a walk by one edge conjugates its magnetic weight and preserves its occupation monomial. Pairing these terms replaces the phase factor by its real part. Absolute convergence follows from $`\|C\|_2<1`$ for $`|z_i|\le1`$.

The measure $`c\mu_C`$ is generally signed. The generating function has nonnegative coefficients because it is the physical Gaussian probability law; the signed measure itself is not a sampling distribution.

These determinant and loop identities have established antecedents in Gaussian photon-count generating functions and Markov loop soups. The GBS cycle expansion also appears explicitly in Lim's lossy-simulation analysis. See [the related-work account](../references/RELATED_WORK.md) and references [1–4] below.

### Gauge invariance

A mode-phase transformation sends $`B_{ij}`$ to $`e^{i(\alpha_i+\alpha_j)}B_{ij}`$. Every vertex phase cancels in the alternating sum around an even walk. All quantities below, including $`r`$, $`\Delta`$, and $`\Gamma`$, are therefore invariant under these transformations.

## 3. Main theorem: negative loop mass

Define

```math
r=\int(-c(w))_+\,d\mu_C(w),
\qquad
\mu_+=c_+\mu_C,
\qquad
\mu_-=(-c)_+\mu_C,
\qquad
\mu_{\mathrm{abs}}=|c|\mu_C.
```

Their total masses satisfy

```math
\mu_-(\mathcal W)=r,
\qquad
\mu_+(\mathcal W)=\Lambda_B+r,
\qquad
\mu_{\mathrm{abs}}(\mathcal W)=\Lambda_B+2r\le\Lambda_C.
```

**Theorem.** Suppose genuine independent Poisson loops with intensity $`\mu_C`$ can be sampled, retaining each loop as a complete object.

1. Retain each source loop independently with probability $`|c(w)|`$. Sum the retained occupations to obtain $`\mathbf n`$, and let $`S=(-1)^{J_-}`$, where $`J_-`$ is the number of retained loops with negative cosine. Then

   ```math
   \mathbb E\!\left[e^{2r}S F(\mathbf n)\right]
   =\mathbb E_{P_B}F,
   \qquad
   \mathbb E\!\left[(e^{2r}S)^2\right]=e^{4r}.
   ```

   Consequently, $`O(e^{4r}\varepsilon^{-2}\log(1/\delta))`$ independent draws suffice for additive error $`\varepsilon`$ with failure probability at most $`\delta`$.

2. Retain each source loop independently with probability $`c(w)_+`$. The resulting occupation law $`Q_+`$ is a positive compound-Poisson distribution and satisfies

   ```math
   \operatorname{TV}(P_B,Q_+)
   \le\min\left\{1,\frac{e^{2r}-1}{2}\right\}.
   ```

Both procedures have source-loop cost bounded by that of the $`C`$ sampler. In particular, positive long loops can be retained; the theorem does not require a length cutoff in the positive part.

### Proof of the estimator

Independent thinning changes a Poisson intensity from $`\mu_C`$ to $`\mu_{\mathrm{abs}}`$. The Poisson exponential formula yields

```math
\begin{aligned}
\mathbb E[S z^{\mathbf n}]
&=\exp\left\{\int\left[
  \operatorname{sgn}(c)z^{v(w)}-1\right]|c|\,d\mu_C\right\}\\
&=\exp\left\{\int\left[c z^{v(w)}-|c|\right]d\mu_C\right\}\\
&=e^{-2r}G_B(z).
\end{aligned}
```

Equality of the absolutely convergent coefficients proves the identity for every bounded $`F`$. The weight magnitude is exactly $`e^{2r}`$, so its second moment is $`e^{4r}`$ and Hoeffding's inequality gives the stated sample count. This is an absolute-weight importance-sampling construction of the kind used in quasiprobability simulation [5]; the representation and efficiently samplable proposal are the application-specific ingredients.

### Proof of the total-variation bound

Push the loop measures forward to occupation vectors. Let $`\delta_0`$ be unit mass at the zero vector, and use $`*`$ for convolution. The signed finite measure

```math
\sigma_-
=e^r\sum_{j=0}^{\infty}\frac{(-1)^j}{j!}\mu_-^{*j}
```

has generating function $`\exp\{-\int(z^v-1)d\mu_-\}`$. Hence the signed-convolution identity is

```math
P_B=Q_+*\sigma_-.
```

The variation norm obeys

```math
\begin{aligned}
\|\sigma_--\delta_0\|_1
&\le(e^r-1)+e^r\sum_{j\ge1}\frac{r^j}{j!}\\
&=e^{2r}-1.
\end{aligned}
```

Convolution with a probability measure contracts this norm. Dividing by two gives the total-variation bound. The physical positivity of $`P_B`$ is already supplied by its Gaussian definition; no positivity inference is made from pointwise closeness of generating functions.

### What the theorem does and does not optimize

The bound $`e^{4r}`$ is a uniform second-moment guarantee for bounded readouts. No claim is made that this estimator minimizes the variance for every readout, among all loop proposals, or among classical algorithms. Different proposals can have different readout-dependent variances.

## 4. Why the Poisson loop partition must be preserved

A sampler that preserves only total occupations is insufficient for the theorem. The negative-binomial excursion representation can forget which elementary excursions belonged to the same Poisson loop. Splitting them changes the estimator because

```math
\cos(\phi_1+\phi_2)\ne\cos\phi_1\cos\phi_2
```

in general, and it also changes absolute-cosine thinning.

For the standard doubled graph

```math
K_C=\begin{pmatrix}0&C\\C&0\end{pmatrix},
```

the loop-soup intensity is $`\alpha=1/2`$. In an elimination construction, let $`q_i`$ be the first-return excursion mass at the current vertex after earlier vertices have been removed. The correct partition is obtained from

```math
J_i\sim\operatorname{Pois}\!\left[-\alpha\log(1-q_i)\right].
```

For each of these $`J_i`$ loops, draw an independent logarithmic-series count

```math
\Pr(K=k)=\frac{q_i^k}{k[-\log(1-q_i)]},\qquad k\ge1,
```

and concatenate $`K`$ independent first-return excursions into **one** loop. Its total alternating phase supplies one cosine and one thinning decision. Vertices with $`q_i=0`$ generate no loops. The factor $`1/k`$ is the cyclic-rooting factor. This is the compound-Poisson form of the standard negative-binomial bubble construction [1,2].

The doubled-graph convention and the base-walk convention in Section 2 agree after pushforward: the two possible starting layers cancel the extra factor $`1/2`$. The [direct-sampler document](direct-sampler.md) specifies the elimination and transition probabilities.

## 5. Two inexpensive determinant bounds

### Domination by the positive reference law

Set

```math
\Delta=\Lambda_C-\Lambda_B
=\int(1-c)\,d\mu_C\ge0.
```

The hafnian triangle inequality gives the pointwise domination

```math
p_B(\mathbf n)\le e^\Delta p_C(\mathbf n).
```

It follows that

```math
\operatorname{TV}(P_B,P_C)
\le1-e^{-\Delta}\le\Delta.
```

The proof is [Proposition 1 in topology and phases](topology.md). This permits sampling directly from $`P_C`$ when $`\Delta`$ is small. It does not give an efficient exact rejection sampler for general $`B`$, because evaluating its acceptance ratio would require additional work.

### The cosine-product estimator and its exact second moment

For an integer $`q\ge0`$, define the entrywise phase harmonic

```math
B^{[q]}_{ij}=
\begin{cases}
C_{ij}e^{iq\arg B_{ij}},&C_{ij}>0,\\
0,&C_{ij}=0.
\end{cases}
```

Every harmonic satisfies $`\|B^{[q]}\|_2\le\rho`$. In particular, $`B^{[0]}=C`$ and $`B^{[1]}=B`$. Write $`D=B^{[2]}`$ and define

```math
\begin{aligned}
\Gamma
&=\frac32\Lambda_C-2\Lambda_B+\frac12\Lambda_D\\
&=\int(1-c)^2\,d\mu_C.
\end{aligned}
```

For the unthinned $`C`$ soup $`\mathcal L`$, with occupations $`\mathbf n=\sum_{w\in\mathcal L}v(w)`$, the weight

```math
W_{\cos}=e^\Delta\prod_{w\in\mathcal L}c(w)
```

satisfies

```math
\mathbb E[W_{\cos}F(\mathbf n)]=\mathbb E_{P_B}F,
\qquad
\mathbb E[W_{\cos}^2]=e^\Gamma.
```

The first identity follows from the Poisson formula with $`g(w)=c(w)z^{v(w)}`$. For the second, use $`\cos^2\phi=(1+\cos2\phi)/2`$ in the same formula. Median-of-means estimation therefore uses $`O(e^\Gamma\varepsilon^{-2}\log(1/\delta))`$ source draws.

The pointwise inequalities $`4(-c)_+\le(1-c)^2\le2(1-c)`$ imply

```math
4r\le\Gamma\le2\Delta\le2\Lambda_C\le\overline N_C.
```

Thus three log determinants give an upper bound on the negative-mass cost. The absolute-cosine method has the smaller uniform second-moment bound, $`e^{4r}\le e^\Gamma`$. This comparison does not assert an ordering of the actual variance for each readout under the two different occupation proposals.

### Quartic sensitivity to small edge phases

Suppose a mode-phase gauge makes every nonzero edge phase satisfy $`|\arg B_{ij}|\le\theta`$. A length-$`L`$ alternating loop then has a real phase representative of magnitude at most $`L\theta`$. Using $`1-\cos x\le x^2/2`$ and $`\operatorname{Tr}C^{2k}\le m\rho^{2k}`$, with $`q=\rho^2`$, gives

```math
\Delta\le\frac{m\theta^2 q}{(1-q)^2},
\qquad
\Gamma\le\frac{2m\theta^4 q(1+4q+q^2)}{(1-q)^4}.
```

For example, the second bound follows from

```math
\Gamma\le\frac{\theta^4}{4}\int |v(w)|^4\,d\mu_C(w)
\le2m\theta^4\sum_{k\ge1}k^3q^k.
```

These orders can occur on connected graphs with large treewidth. Let $`H`$ be a fixed-degree $`d`$-regular expander on $`n`$ vertices and $`G=H\mathbin{\square}K_2`$, so $`m=2n`$ and $`G`$ has degree $`d+1`$. Take $`C=tA_G/(d+1)`$ at fixed $`0<t<1`$. Put phase $`\theta`$ on every edge inside either copy of $`H`$, and phase zero on the rungs. Each of the $`md/4`$ squares associated with an edge of $`H`$ has alternating phase $`\pm2\theta`$ and total rooted-loop intensity $`2[t/(d+1)]^4`$. Consequently,

```math
\begin{aligned}
\Delta&\ge\frac{md}{2}\left(\frac{t}{d+1}\right)^4(1-\cos2\theta),\\
\Gamma&\ge\frac{md}{2}\left(\frac{t}{d+1}\right)^4(1-\cos2\theta)^2.
\end{aligned}
```

Together with the upper bounds, these yield $`\Delta=\Theta(m\theta^2)`$ and $`\Gamma=\Theta(m\theta^4)`$ as $`\theta\to0`$. The physical mean photon number is $`\Theta(m)`$, since $`\operatorname{Tr}(BB^\dagger)=mt^2/(d+1)`$ and $`\|B\|_2\le t`$. The expander is a graph minor of $`G`$, so taking the Cartesian product does not remove its large support treewidth.

At $`\theta=m^{-1/4}`$, the cosine estimator has bounded second moment even though $`\Delta=\Theta(\sqrt m)`$. This demonstrates a variance-certificate improvement over a bound exponential in $`\Delta`$. It is **not** a separation between estimation and full sampling: loops shorter than $`\pi/(2\theta)`$ have nonnegative cosine, and the exponentially small long-loop tail at fixed $`t`$ also certifies efficient positive-part sampling for this family. The corresponding strong-sampling argument is given in [topology and phases](topology.md).

## 6. A Fourier certificate with an explicit finite truncation error

The absolutely convergent Fourier expansion

```math
|\cos\phi|
=\frac2\pi+\frac4\pi\sum_{k=1}^{\infty}
\frac{(-1)^{k-1}\cos(2k\phi)}{4k^2-1}
```

can be integrated against the finite loop measure. Since

```math
\int\cos(q\phi)\,d\mu_C=\Lambda_{B^{[q]}},
```

one obtains

```math
r=\frac12\left[
\frac2\pi\Lambda_C
+\frac4\pi\sum_{k=1}^{\infty}
\frac{(-1)^{k-1}\Lambda_{B^{[2k]}}}{4k^2-1}
-\Lambda_B\right].
```

Let $`r_K`$ retain the first $`K`$ harmonics in this expression. Each harmonic log determinant lies between zero and $`\Lambda_C`$. The telescoping identity

```math
\sum_{k>K}\frac1{4k^2-1}=\frac1{2(2K+1)}
```

therefore proves the exact error bound

```math
|r-r_K|\le\frac{\Lambda_C}{\pi(2K+1)}.
```

This is a finite approximation to $`r`$, not an assertion that a fixed number of determinants determines $`r`$ exactly.

### Determinant and input errors

Suppose every log determinant entering a computed $`\widetilde r_K`$ has certified absolute error at most $`\zeta`$. This budget must include error in forming the harmonic matrix, as well as the determinant calculation itself. Let $`\Lambda_C^U`$ be a certified upper bound. Then

```math
|r-\widetilde r_K|\le E_K,
\qquad
E_K=\frac{\Lambda_C^U}{\pi(2K+1)}
+\left(\frac12+\frac2\pi\right)\zeta.
```

The numerical coefficient is the sum of the absolute Fourier coefficients, including the $`\Lambda_C`$ and $`\Lambda_B`$ terms. A valid interval is

```math
\left[\max\{0,\widetilde r_K-E_K\},\ \widetilde r_K+E_K\right],
```

and it can be intersected with any separately certified upper bound $`\Gamma^U/4`$. If the claimed certificates produce an empty intersection, that inconsistency must be reported rather than hidden by clipping. When each of the three log determinants defining $`\Gamma`$ has error at most $`\zeta`$, its absolute error is at most $`4\zeta`$.

Dense evaluation uses $`O((K+1)m^3)`$ arithmetic operations. Additive accuracy $`\eta`$ is achieved with $`K=O(1+\Lambda_C/\eta)`$ and determinant error $`O(\eta)`$. Cancellation makes **absolute** accuracy the relevant criterion, even when $`r`$ or $`\Gamma`$ vanishes. Certified bit costs also require a spectral-gap certificate and controlled input approximation; see [finite precision](finite-precision.md).

### Normalization and a ratio alternative

Replacing $`r`$ by $`\widehat r`$ changes the expectation by a known scalar factor:

```math
\mathbb E[e^{2\widehat r}S F]
=e^{2(\widehat r-r)}\mathbb E_{P_B}F.
```

If $`|\widehat r-r|\le\eta`$, the absolute bias is at most $`e^{2\eta}-1`$. Thus additive $`O(\varepsilon)`$ accuracy in $`r`$ is sufficient. An upper bound alone cannot be substituted into the exponential while retaining unbiasedness.

A ratio estimator avoids the normalizer entirely:

```math
\mathbb E_{P_B}F=\frac{\mathbb E[S F]}{\mathbb E[S]},
\qquad \mathbb E[S]=e^{-2r}.
```

Given a certified $`R\ge r`$, set $`d_0=e^{-2R}`$. Estimate both expectations to absolute error at most $`\tau=\varepsilon d_0/4`$. If their empirical values are $`\widehat A,\widehat d`$, return zero when $`\widehat d<d_0/2`$; otherwise return $`\widehat A/\widehat d`$, clipped to $`[-1,1]`$. On the concentration event the fallback is never used, and

```math
\left|\frac{\widehat A}{\widehat d}-\mathbb E_{P_B}F\right|
\le\frac{2\tau}{d-\tau}\le\varepsilon,
\qquad d=e^{-2r},\quad 0<\varepsilon\le1.
```

Hoeffding's inequality and a union bound give a sufficient count $`32e^{4R}\varepsilon^{-2}\log(4/\delta)`$. The finite-sample ratio is not unbiased. A constant-additive certified upper approximation to $`r`$ changes this cost by only a constant factor, so its Fourier certificate can use $`O(1+\Lambda_C)`$ harmonics, independent of the final $`\varepsilon`$. Direct exponential normalization instead requires the additive $`O(\varepsilon)`$ approximation discussed above.

## 7. End-to-end interpretation

The dense matrix must be supplied explicitly, with $`\Omega(m^2)`$ input reading and storage access. The elimination sampler uses $`O(m^3)`$ reusable preprocessing and $`O(m^2)`$ storage. In ideal arithmetic, a genuine source-loop draw has expected cost $`O(m(1+\overline N_C))`$ when transitions are sampled by scanning their categorical probabilities. The expected number of loop visits is $`\overline N_C`$; Poisson and logarithmic-series generation must be included in the implementation. Computing phases and thinning loops fits within this bound.

For bounded readouts, the absolute-cosine method consequently has the arithmetic scaling

```math
O\!\left(
m^3+
e^{4r}\varepsilon^{-2}\log(1/\delta)
\left[m(1+\overline N_C)+T_F\right]
\right),
```

plus the chosen certificate/normalization computation. Here $`T_F`$ is the evaluation cost under the proposal, or a uniform bound after an explicitly controlled output cutoff. Reusable preprocessing is paid once for repeated draws from the same matrix. A Fourier normalizer with $`O(\varepsilon)`$ accuracy adds $`O(m^3(1+\Lambda_C/\varepsilon))`$ arithmetic operations.

Polynomial $`\overline N_C`$, $`r=O(\log m)`$, polynomially evaluable readout, and the finite-precision hypotheses together give a polynomial algorithm. A small physical mean $`\overline N_B`$ alone does not establish this claim. Neither exact-real arithmetic nor unit-cost oracle access to the input is assumed in the certified statement. The repository's floating-point reference code illustrates the identities; its outputs alone are not interval-certified determinant or total-variation bounds.

## 8. Classical output processing, loss, and conditioning

An efficiently samplable classical output channel $`K(y\mid\mathbf n)`$ can be applied to every proposal draw, followed by a bounded readout $`F(y)`$. The same weighted identities hold by conditioning on $`\mathbf n`$. Total variation contracts under the channel, so the positive-part sampling bound also persists.

This includes threshold detection and independent loss **after the modeled Gaussian circuit**:

```math
y_i\mid n_i\sim\operatorname{Binom}(n_i,\eta_i).
```

Arbitrary losses elsewhere in an interferometer are not silently replaced by this output channel. Their equivalence requires a separate optical model or commutation argument.

For a postselection event $`A`$ of physical probability $`a=P_B(A)>0`$,

```math
\mathbb E_{P_B}[F\mid A]
=\frac{\mathbb E[S F\mathbf1_A]}{\mathbb E[S\mathbf1_A]}.
```

The denominator is $`e^{-2r}a`$. A certified lower bound $`a\ge a_0>0`$ gives the sufficient signed-estimation sample count

```math
O\!\left(e^{4r}a_0^{-2}\varepsilon^{-2}\log(1/\delta)\right),
```

or the corresponding bound with a certified upper $`R\ge r`$. For example, errors at most $`\varepsilon e^{-2r}a/4`$ in numerator and denominator give conditional error at most $`\varepsilon`$ when $`0<\varepsilon\le1`$. The normalizer cancels in this ratio.

For a positive approximate sampler, if $`\operatorname{TV}(P_B,Q)\le\tau<a`$, then $`Q(A)>0`$ and the safe bound

```math
\operatorname{TV}\bigl(P_B(\cdot\mid A),Q(\cdot\mid A)\bigr)
\le\frac{2\tau}{a}
```

holds. Exact rejection sampling, when available, has expected overhead $`1/a`$. Fixed-photon-number or collision-free conditioning can therefore be costly even when the full distribution is easy. No efficiency claim for exponentially rare conditioning follows from the main theorem.

## References

1. G. F. Lawler, [*Topics in loop measures and the loop-erased walk*](https://www.math.uchicago.edu/~lawler/fall16.pdf), notes (2017), Sections 5.3–5.4.
2. Y. Le Jan, [*Markov loops, determinants and Gaussian fields*](https://arxiv.org/abs/math/0612112), arXiv:math/0612112v4 (2007).
3. K. R. Parthasarathy and R. Sengupta, [*From particle counting to Gaussian tomography*](https://arxiv.org/abs/1504.07054), *Infinite Dimensional Analysis, Quantum Probability and Related Topics* **18**, 1550023 (2015).
4. Y. Lim, [*Classical Simulation of Lossy Gaussian Boson Sampling beyond Square-Root Regime*](https://arxiv.org/html/2610.06526v1), arXiv:2610.06526v1 (2026), Supplement SIV.C.
5. H. Pashayan, J. J. Wallman, and S. D. Bartlett, [*Estimating outcome probabilities of quantum circuits using quasiprobabilities*](https://arxiv.org/html/1503.07525v2), *Physical Review Letters* **115**, 070501 (2015), Eqs. (8)–(9).
