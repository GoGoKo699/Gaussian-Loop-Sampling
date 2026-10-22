# Finite-bit sampling and phase estimation

This document gives a bounded-time implementation of the Gaussian loop representation using rational arithmetic and finitely many random bits. It separates three statements: the unbiased ideal estimator, approximation of the positive reference distribution, and the bias-certified finite-bit estimator. The implementation below uses matrix powers. Its preprocessing bound is polynomial; it is distinct from the cubic expected-arithmetic bound of the ordered-excursion implementation.

## 1. Input, output, and charged resources

Let $`C`$ be a symmetric, entrywise nonnegative $`m\times m`$ rational matrix, and let $`\rho=\|C\|_2<1`$. A nonzero entry is represented by an integer numerator and a positive integer denominator, each of at most $`b`$ bits. Zero entries are represented exactly. Supply an integer

```math
M\ge \max\{1,N_C\},\qquad
N_C=\operatorname{Tr}\!\left[C^2(I-C^2)^{-1}\right].
```

The running-time bounds are polynomial in the *value* of $`M`$, not merely its encoding length. One may compute $`N_C`$ by exact rational linear algebra, or supply an independently certified upper bound. If $`\|C\|_2\le1-\kappa`$ is known instead, $`M=\lceil m/\kappa\rceil`$ is a valid conservative choice for $`0<\kappa\le1`$.

For phase estimation, the target is $`B_{ij}=C_{ij}u_{ij}`$, with $`u_{ij}=u_{ji}`$ and $`|u_{ij}|=1`$ whenever $`C_{ij}>0`$. A deterministic phase oracle returns a rational complex number within absolute error $`2^{-p}`$ of $`u_{ij}`$ in time $`T_u(p)`$. Its output encoding length is included in that time. Values on zero entries are immaterial. This includes rational magnitudes with phases specified by polynomial-time computable angles. Section 11 extends the result to a general complex-entry oracle; no exact polar decomposition of the original input is assumed there.

An observable $`f:\mathbb N^m\to[-1,1]`$ is supplied by an algorithm which evaluates it to absolute error $`2^{-p}`$ in time $`T_f(p,F)`$ on strings with at most $`F`$ total photons. An exact rational-valued observable is a special case. Every oracle call, integer operation, random bit, and emitted photon is charged.

The accuracy and failure parameters are rational numbers in $`(0,1)`$. The target output distribution is the full, undisplaced, pure Gaussian photon-count law, including collisions and its fluctuating total count. No photon-sector or collision-free conditioning is included in these costs.

## 2. Ideal identities used by the implementation

For a symmetric complex matrix $`X`$ with $`\|X\|_2<1`$, write

```math
\Lambda_X=-\frac12\log\det(I-XX^\dagger).
```

The positive reference process is a Poisson process of even closed walks. A rooted walk $`\ell=(i_0,i_1,\ldots,i_{2k}=i_0)`$ has intensity

```math
\mu_C(\ell)=\frac1{2k}\prod_{t=0}^{2k-1}C_{i_ti_{t+1}}.
```

The formula sums over all rooted sequences. Its redundancy under cyclic shifts is intentional: this representation of the Poisson mark space produces the usual unrooted loop measure. Equivalently, it is the intensity-$`1/2`$ loop soup on the doubled kernel $`\begin{psmallmatrix}0&C\\C&0\end{psmallmatrix}`$, with the two copies of each mode aggregated. Let $`v(\ell)`$ count occurrences of each vertex among $`i_0,\ldots,i_{2k-1}`$. The total intensity and total mean occupation are

```math
\mu_C(\text{all loops})=\Lambda_C,
\qquad
\int |v(\ell)|\,d\mu_C(\ell)=N_C.
```

For a loop, use the alternating edge product

```math
z_B(\ell)=
\prod_{\substack{0\le t<2k\\t\ \mathrm{even}}}u_{i_ti_{t+1}}
\prod_{\substack{0\le t<2k\\t\ \mathrm{odd}}}\overline{u_{i_ti_{t+1}}},
\qquad c_B(\ell)=\operatorname{Re}z_B(\ell).
```

Choosing the other initial layer conjugates $`z_B`$ and leaves $`c_B`$ unchanged. Reversing a loop has the same property. Therefore no extra layer coin is needed for the cosine estimator.

Let $`\mathcal L`$ be the complete positive loop process and $`n=\sum_{\ell\in\mathcal L}v(\ell)`$. Define

```math
\Delta=\Lambda_C-\Lambda_B,
\qquad
W=e^\Delta\prod_{\ell\in\mathcal L}c_B(\ell).
```

The ideal identities are

```math
\mathbb E[Wf(n)]=\mathbb E_{P_B}f,
\qquad
\mathbb E W^2=e^\Gamma,
```

where

```math
D_{ij}=C_{ij}u_{ij}^2,
\qquad
\Gamma=\frac32\Lambda_C-2\Lambda_B+\frac12\Lambda_D
      =\int(1-c_B(\ell))^2\,d\mu_C(\ell).
```

These follow from the Poisson exponential formula and the Gaussian photon generating function. In particular,

```math
0\le\Gamma\le2\Delta\le2\Lambda_C\le N_C\le M,
\qquad
|W|\le e^\Delta\le e^{M/2}.
```

The unbiased identity applies to genuine Poisson loops. Elementary first-return excursions must be grouped into their Poisson clusters before taking one cosine per loop. The finite-power construction below generates genuine loops directly and does not have this grouping ambiguity.

## 3. Exponential tails with explicit cutoffs

Put $`A=M+1`$. Since the largest squared eigenvalue contributes $`\rho^2/(1-\rho^2)`$ to $`N_C`$,

```math
\rho^2\le\frac{M}{M+1}=1-\frac1A.
```

### 3.1 Omitted loop lengths

Let $`\tau_K`$ be the intensity of loops longer than $`2K`$. Then

```math
\begin{aligned}
\tau_K
&=\sum_{k>K}\frac{\operatorname{Tr}(C^{2k})}{2k}\\
&\le \frac{m\rho^{2(K+1)}}{2(K+1)(1-\rho^2)}
\le\frac{mA}{2}\exp\!\left(-\frac{K+1}{A}\right).
\end{aligned}
```

Consequently the probability that the full soup contains an omitted loop is at most $`\tau_K`$. For any $`0<\tau_L<1`$, the entirely integer-valued choice

```math
K=A\left\lceil\log_2\frac{mA}{\tau_L}\right\rceil
```

ensures $`\tau_K\le\tau_L`$. The ceiling of the binary logarithm is obtained by integer comparisons on the rational input; evaluating a transcendental logarithm is unnecessary.

### 3.2 Total photons

Write $`N=\sum_i n_i`$ and set $`z=1+1/(2A)`$. If $`r_i`$ are the eigenvalues of $`C^2`$ and $`t_i=r_i/(1-r_i)`$, then $`\sum_i t_i=N_C\le M`$ and

```math
\begin{aligned}
\log\mathbb E z^{N/2}
&=-\frac12\sum_i\log\left(1-\frac{t_i}{2A}\right)\\
&\le\sum_i\frac{t_i}{2A}
\le\frac12.
\end{aligned}
```

Here $`0\le t_i/(2A)<1/2`$ and $`-\log(1-x)\le2x`$ were used. Also $`\log z\ge1/(4A)`$. Therefore

```math
\Pr(N>F)\le\exp\left(\frac12-\frac{F}{8A}\right).
```

For $`0<\tau_N<1`$, take the even integer

```math
F=8A\left(1+\left\lceil\log_2\frac1{\tau_N}\right\rceil\right).
```

Then $`\Pr(N>F)\le\tau_N`$. Deleting long loops only decreases the photon count under the natural coupling, so the same bound applies after the length cutoff. These exponential tails are essential when $`\tau_L`$ or $`\tau_N`$ is exponentially small in $`M`$; a Markov bound would not provide polynomial cutoffs.

## 4. Exact rational closed-walk bridges

Precompute $`P_r=C^r`$ for $`0\le r\le2K`$ and the rational numbers

```math
\lambda_k=\frac{\operatorname{Tr}P_{2k}}{2k},
\qquad
\lambda=\sum_{k=1}^K\lambda_k\le\Lambda_C\le M/2.
```

If $`\lambda=0`$, the truncated process is empty. Otherwise its exact ideal construction is:

1. Draw $`J\sim\operatorname{Poisson}(\lambda)`$.
2. For each of the $`J`$ marks, draw $`k`$ with probability $`\lambda_k/\lambda`$.
3. Draw a root $`i`$ with probability $`(P_{2k})_{ii}/\operatorname{Tr}P_{2k}`$.
4. Starting at $`x=i`$ with $`r=2k`$ remaining edges, draw the next vertex $`y`$ with probability

```math
\Pr(y\mid x,r,i)=\frac{C_{xy}(P_{r-1})_{yi}}{(P_r)_{xi}},
```

and repeat with $`x=y`$, $`r\leftarrow r-1`$. The endpoint after $`2k`$ steps is $`i`$ with probability one. Record the $`2k`$ departure vertices and the alternating edge phases.

The bridge row sums to one by matrix multiplication. Multiplying its conditional probabilities telescopes to the product of the edge weights divided by $`(P_{2k})_{ii}`$. Including the root and length probabilities therefore assigns each rooted walk its intensity divided by $`\lambda`$, as required.

There is no numerical division by an uncertain small number. Every denominator is rational and exactly positive on a reachable bridge state. A root of zero weight is never selected; a zero transition has exactly zero probability. If a transition of positive probability is taken, its remaining bridge weight is positive, so this invariant continues inductively.

### Encoding sizes

Let $`Q`$ be the product of the denominators of the entries of $`C`$ and write $`C=H/Q`$ with integer $`H`$. Then $`\log_2Q\le m^2b`$, and $`0\le H_{ij}<Q`$ because $`\|C\|_2<1`$. An entry of $`H^r`$ has at most

```math
O\bigl(r(m^2b+\log m)\bigr)
```

bits. A common denominator for all $`\lambda_k`$ and their sum is $`2K!Q^{2K}`$. Thus their encoding lengths are

```math
B_\lambda=O\bigl(K(m^2b+\log m+\log K)\bigr).
```

Bridge probabilities are ratios of integers of comparable polynomial length. Their small numerical values introduce no dependence on inverse transition probabilities. Powers are computed by ordinary integer matrix multiplication, so this is a constructive polynomial bit bound.

## 5. Bounded finite-randomness sampling

The following implementation approximates the exact truncated process, with the rule that it returns a distinguished failure symbol if its total photon count exceeds $`F`$. Failure is later mapped to the vacuum for distribution sampling, or to the value zero for estimation.

Let $`H_0=F/2`$. Since every loop has at least two photons, $`J>H_0`$ implies a photon-cap failure. Accordingly, only the Poisson probabilities for $`0\le j\le H_0`$ and one overflow outcome are needed. If a newly drawn loop length would put the total above $`F`$, return failure before constructing its bridge. Consequently there are at most $`F`$ bridge edges, $`F/2+1`$ length draws, and $`F/2`$ root draws. A conservative bound on all categorical calls, including the Poisson draw, is

```math
R_0=2F+4.
```

### 5.1 Categorical sampling without relative-error assumptions

For any exact rational categorical distribution on at most $`S`$ outcomes, draw a uniform integer $`U\in\{0,\ldots,2^r-1\}`$ and compare $`U/2^r`$ with its exact cumulative probabilities, assigning boundary points to the interval on their right. The law differs from the desired categorical law by at most $`(S-1)2^{-r}`$ in total variation. To see this, couple the draw to a uniform real number by rounding it down; disagreement can occur only in one grid cell adjacent to each cumulative boundary.

Zero-probability intervals have zero length and are never selected, including when several cumulative boundaries coincide. In particular, finite randomness does not create a path whose exact bridge denominator vanishes.

Let $`S=\max\{K,m,H_0+2\}`$. Choosing

```math
r=\left\lceil\log_2\frac{4R_0S}{\eta_{\rm num}}\right\rceil
```

limits the accumulated categorical error to $`\eta_{\rm num}/4`$. This adaptive union bound applies because the same per-call bound holds at every reachable state. It does not require a lower bound on a positive transition probability. At most $`R_0r`$ random bits are used.

### 5.2 A rational Poisson approximation

For $`\lambda>0`$, compute a rational lower approximation $`q_0\le e^{-\lambda}`$ with

```math
0\le e^{-\lambda}-q_0\le2^{-p},
\qquad
p\ge M+\left\lceil\log_2\frac4{\eta_{\rm num}}\right\rceil.
```

Set

```math
q_j=q_0\frac{\lambda^j}{j!}\quad(0\le j\le H_0),
\qquad
q_{\rm over}=1-\sum_{j=0}^{H_0}q_j.
```

These are nonnegative exact rationals. Relative to the exact Poisson law aggregated at overflow, all removed probability is transferred to overflow. Its total variation error is therefore at most

```math
2^{-p}\sum_{j=0}^{H_0}\frac{\lambda^j}{j!}
\le2^{-p}e^{M/2}\le\eta_{\rm num}/4.
```

For $`\lambda=0`$, output zero exactly. Sampling the rational vector by the preceding grid rule adds only its already budgeted categorical error.

For completeness, $`q_0`$ can be produced without an exact exponential primitive. Let

```math
S_L=\sum_{j=0}^{L}\frac{\lambda^j}{j!},
\qquad L\ge\lceil12\lambda\rceil+p+2.
```

The exponential-series tail is at most $`2^{-p}`$: beyond this index successive terms decrease by at least a factor two, and

```math
\sum_{j>L}\frac{\lambda^j}{j!}
\le 2\frac{\lambda^{L+1}}{(L+1)!}
\le 2\left(\frac{e\lambda}{L+1}\right)^{L+1}
\le2^{-p}.
```

Then $`q_0=1/(S_L+2^{-p})`$ is a valid lower approximation with the stated error. All quantities are rational; their bit lengths are polynomial in $`L,B_\lambda,H_0`$. In particular, forming the probabilities through $`j=H_0`$ may require $`O(H_0(B_\lambda+\log H_0))`$ bits beyond the encoding of $`q_0`$. This remains polynomial because $`H_0=F/2`$ is polynomially bounded. The construction uses $`O(M+\log(1/\eta_{\rm num}))`$ Taylor terms, not an exponential number of trials.

### 5.3 Resulting coupling guarantee

Combining the Poisson and categorical errors gives, with slack,

```math
d_{\rm TV}(\widehat{\mathcal L}_{K,F},\mathcal L_{K,F})
\le\eta_{\rm num}.
```

Here $`\mathcal L_{K,F}`$ denotes the exact short-loop process with the photon-cap failure rule. The construction has a deterministic polynomial work bound. Coupling it to the full loop process adds at most $`\tau_K+\Pr(N>F)`$.

It follows immediately that the positive photon law $`P_C`$ can be sampled within any total-variation error $`\delta`$ by using $`\tau_L=\tau_N=\eta_{\rm num}=\delta/3`$ and mapping failure to the vacuum. The bit complexity is polynomial in

```math
m,\ b,\ M,\ \log(1/\delta).
```

The output itself is bounded by $`F=O((M+1)(1+\log(1/\delta)))`$ photons. This theorem does not claim exact sampling with a fixed finite number of random bits.

## 6. Certified determinants, normalizers, and phase arithmetic

No relative accuracy of $`\Delta`$, $`\Gamma`$, a determinant, or a possibly vanishing entry is required. All certificates below use absolute error.

### 6.1 Log-determinants by a convergent positive trace series

For $`X=B,C,D`$,

```math
\Lambda_X=\frac12\sum_{j\ge1}
\frac{\operatorname{Tr}[(XX^\dagger)^j]}j.
```

Every one of these matrices has norm at most $`\rho`$. The tail beyond $`L`$ is bounded by

```math
\frac{mA}{2}\exp\left(-\frac{L+1}{A}\right).
```

For phase-oracle inputs, first ask for $`v_{ij}`$ with $`|v_{ij}-u_{ij}|\le a`$ and use

```math
\widetilde u_{ij}=\frac{v_{ij}}{1+a}.
```

Then $`|\widetilde u_{ij}|\le1`$ and $`|\widetilde u_{ij}-u_{ij}|\le2a`$. Use the same value for $`ij`$ and $`ji`$. The rational matrices $`\widetilde B=C\circ\widetilde u`$ and $`\widetilde D=C\circ\widetilde u^2`$ still satisfy the norm bound $`\rho`$, because their entrywise magnitudes are bounded by $`C`$.

If $`\|X\|,\|Y\|\le\rho`$, the trace series and the power-difference identity give

```math
|\Lambda_X-\Lambda_Y|
\le \frac{m\rho}{1-\rho^2}\,\|X-Y\|_2.
```

Moreover, $`\|B-\widetilde B\|_2\le2a\rho`$ and $`\|D-\widetilde D\|_2\le4a\rho`$. To compute each $`\Lambda_X`$ within absolute error $`\xi`$, it therefore suffices to choose

```math
a\le\frac{\xi}{8mA},
\qquad
L=A\left\lceil\log_2\frac{mA}{\xi}\right\rceil,
```

compute the rational truncated trace series exactly, and retain the explicit error interval. This has polynomial bit cost. The sum of absolute coefficients in $`\Gamma=\tfrac32\Lambda_C-2\Lambda_B+\tfrac12\Lambda_D`$ is four, so intervals of width controlled at absolute error $`\xi`$ also control cancellation in $`\Gamma`$.

Compute a rational upper bound $`G`$ satisfying

```math
\Gamma\le G\le\Gamma+1.
```

A rational upper approximation $`V`$ to $`e^G`$ can then be chosen so that

```math
e^\Gamma\le V\le2e^{\Gamma+1}.
```

The same exponential-series algorithm provides this value in polynomial time in $`M`$. Although the value of $`V`$ can be exponential, its binary encoding has $`O(M)`$ bits. The eventual number of estimator repetitions depends on its value, as stated explicitly below.

### 6.2 Deterministic error in one sample value

Let $`R=2^M`$, so $`e^\Delta\le R`$. On a record with at most $`F`$ edges, the estimator can be evaluated to any uniform absolute error $`\eta`$ with polynomial bit cost:

- Compute $`\Delta`$ to absolute error at most $`\eta/(4R)`$ and clip this approximation to $`[0,M/2]`$. Evaluate its exponential to absolute error at most $`\eta/4`$, then clip the result to $`[0,R]`$. The resulting normalizer differs from $`e^\Delta`$ by at most $`\eta/2`$.
- Evaluate each unit edge phase to error $`a\le\eta/(8RF)`$ and contract it as above. Compute products as exact complex rationals. A loop of length $`\ell`$ has phase-product error at most $`2\ell a`$. Because every computed phase product has modulus at most one, telescoping across all loop cosines gives total product error at most $`2Fa`$.
- Evaluate $`f`$ to error $`\eta/(4R)`$ and clip its value to $`[-1,1]`$.

The sum of the resulting three errors is at most $`\eta`$, and the computed value has absolute value at most $`R`$. Exact rational products across at most $`F`$ edges have polynomial encoding length. Thus this prescription does not presume fixed-cost real or complex multiplication.

The required phase and observable precision has

```math
O\bigl(M+\log F+\log(1/\eta)\bigr)
```

bits. Computing $`\Delta`$ to the needed absolute precision requires a trace-series cutoff of order

```math
O\bigl((M+1)[M+\log(m(M+1)/\eta)]\bigr).
```

These are polynomial precision costs even when $`e^\Delta`$ is large.

## 7. Preserving the second-moment advantage under truncation

Using a small total-variation error alone to control an unbounded importance estimator would be insufficient. The distinction between truncation error and numerical sampling error resolves this issue.

Let $`\mathcal L_K`$ contain only loops of length at most $`2K`$, and retain the *full* normalizer $`e^\Delta`$:

```math
W_K=e^\Delta\prod_{\ell\in\mathcal L_K}c_B(\ell),
\qquad
X_K=W_Kf(n_K).
```

The omitted loop process is independent of $`\mathcal L_K`$. By the Poisson formula,

```math
\mathbb E W_K^2
=e^\Gamma\exp\left[\int_{|v|>2K}(1-c_B(\ell)^2)\,d\mu_C(\ell)\right]
\le e^{\Gamma+\tau_K}.
```

Couple the full and truncated soups. Their sample values agree if no omitted loop occurs. Cauchy–Schwarz then gives

```math
|\mathbb E[Wf]-\mathbb E X_K|
\le \sqrt{e^\Gamma\tau_K}
   +\sqrt{e^{\Gamma+\tau_K}\tau_K}.
```

Define $`X_{K,F}=X_K\mathbf1_{\{N_K\le F\}}`$. Since $`N_K\le N`$ under coupling,

```math
|\mathbb E X_K-\mathbb E X_{K,F}|
\le \sqrt{e^{\Gamma+\tau_K}\Pr(N>F)},
\qquad
\mathbb E X_{K,F}^2\le e^{\Gamma+\tau_K}.
```

These estimates justify tail budgets of order $`\epsilon^2e^{-\Gamma}`$, even when $`\Delta`$ is much larger. They apply to the original long-loop deletion and photon cap, not indiscriminately to any numerical perturbation of their distributions.

For numerical sampling, the uniform bound $`|X_{K,F}|\le R`$ is used instead. A coupling failure of probability $`\eta_{\rm num}`$ changes the bias by at most $`2R\eta_{\rm num}`$ and can add at most $`R^2\eta_{\rm num}`$ to a second-moment bound. Hence the numerical probability budget must also protect the variance. Its precision cost remains polynomial because $`\log R=M\log2`$.

## 8. Finite-bit estimator theorem

**Theorem.** Under the input model in Section 1, for $`0<\epsilon,\beta<1`$ there is a randomized finite-bit algorithm whose output $`\widehat a`$ satisfies

```math
\Pr\left[|\widehat a-\mathbb E_{P_B}f|>\epsilon\right]\le\beta.
```

It uses

```math
O\left(e^\Gamma\epsilon^{-2}\log(1/\beta)\right)
```

independent reference draws. Its remaining preprocessing and per-draw bit costs are polynomial in

```math
m,\ b,\ M,\ \log(1/\epsilon),\ \log(1/\beta),
```

and in the specified phase/observable oracle costs. In particular, for polynomial $`M`$ and $`\Gamma=O(\log m)`$ it is a polynomial-time additive estimator at inverse-polynomial accuracy and confidence. The ideal estimator is unbiased; the finite-bit estimator has an explicitly bounded systematic bias.

### Proof with explicit budgets

Compute $`V`$ as in Section 6 and set

```math
\theta=\frac{\epsilon^2}{2^{14}V},
\qquad
\tau_L=\tau_N=\theta,
\qquad
\eta_{\rm num}=\frac{\epsilon^2}{2^{14}R^2},
\qquad
\eta=\epsilon/128.
```

Use Sections 3–5 with these values. Evaluate each successful record as in Section 6 to error $`\eta/2`$, and round to a common dyadic grid with rounding error at most $`\eta/2`$. Clip to $`[-R,R]`$, which preserves this bound. Return zero on failure. Let $`Y`$ denote the resulting dyadic random value; its total evaluation error is at most $`\eta`$.

Since $`\theta<1/2`$, Section 7 gives

```math
\mathbb E X_{K,F}^2\le2V,
\qquad
|\mathbb E X_{K,F}-\mathbb E_{P_B}f|
\le4\sqrt{V\theta}=\epsilon/32.
```

Couple the finite record to the exact capped record. On agreement the evaluation error is at most $`\eta`$; on disagreement both values are bounded by $`R`$. Therefore

```math
|\mathbb E Y-\mathbb E X_{K,F}|
\le\eta+2R\eta_{\rm num}<\epsilon/64,
```

and

```math
\mathbb E Y^2
\le2\mathbb E X_{K,F}^2+2\eta^2+R^2\eta_{\rm num}
\le5V.
```

Thus $`|\mathbb E Y-\mathbb E_{P_B}f|<\epsilon/8`$ while $`\operatorname{Var}Y\le5V`$. Divide independent values into groups of

```math
s=\left\lceil80V/\epsilon^2\right\rceil
```

draws. Chebyshev's inequality makes the probability that one group mean differs from $`\mathbb E Y`$ by more than $`\epsilon/2`$ at most $`1/4`$. The median of an odd number $`t\ge8\log(1/\beta)`$ of group means then fails with probability at most $`e^{-t/8}\le\beta`$, by Hoeffding's inequality for the independent bad-group indicators. Adding the systematic bias gives the claimed accuracy.

All sample values have a common denominator of $`O(\log(1/\eta))`$ bits. A group sum therefore needs only $`O(M+\log(1/\epsilon)+\log s)`$ bits, rather than a product of unrelated sample denominators. Integer accumulation, comparison of equal-size group sums, and one final division implement the median-of-means step in time linear in the number of draws up to polynomial bit factors.

### Sizes of the cutoffs and precisions

The selected parameters give

```math
\begin{aligned}
K&=O\bigl((M+1)[\log(m(M+1))+\Gamma+\log(1/\epsilon)]\bigr),\\
F&=O\bigl((M+1)[1+\Gamma+\log(1/\epsilon)]\bigr),\\
\log(1/\eta_{\rm num})&=O(M+\log(1/\epsilon)).
\end{aligned}
```

Thus $`\Gamma`$ controls the repetition count and the loop/photon tails, while the potentially larger $`\Delta`$ is paid through a polynomial precision cost. Replacing $`R=2^M`$ by a certified tighter upper bound on $`e^\Delta`$ is permissible but unnecessary for the theorem.

## 9. Arithmetic complexity and bit complexity are different guarantees

The ordered-excursion construction of the positive soup has $`O(m^3)`$ preprocessing arithmetic, $`O(m^2)`$ storage, and expected $`O(m+mN_C)`$ arithmetic per sample, when exact scalar random primitives and exact real operations are used. That statement belongs to the ideal arithmetic model.

The implementation proved here deliberately uses a simpler finite-power construction. It requires $`O(Km^3)`$ matrix arithmetic to precompute powers and $`O(Km^2)`$ storage before accounting for integer lengths. Length and root cumulative tables can be precomputed; their draws use binary searches. Each bridge step scans at most $`m`$ weights. Poisson and phase tables are reused across repetitions.

There are polynomially many additional rational operations to construct the finite Poisson law, determinant certificates, and phase approximations. An individual draw uses at most

```math
O\bigl(F(m+\log K)+\log F\bigr)
```

categorical/bridge arithmetic operations, followed by at most $`F`$ rational phase multiplications and one observable evaluation. The integer lengths were bounded in Sections 4–6. Matrix-trace preprocessing for the normalizer may use a cutoff of order $`(M+1)[M+\log(m(M+1)/\epsilon)]`$.

Consequently there are polynomials $`P_{\rm pre}`$ and $`P_{\rm draw}`$, furnished by the explicit procedures above, for which the complete bit cost is bounded by

```math
P_{\rm pre}
+O\left(V\epsilon^{-2}\log(1/\beta)\right)P_{\rm draw},
```

where their arguments are the parameters and oracle costs stated in the theorem. This does **not** establish cubic preprocessing in the bit model. It establishes that the ideal probabilistic identities and the $`e^\Gamma`$ estimation advantage admit a complete polynomial-overhead implementation.

## 10. Bounded absolute-cosine thinning and the ratio estimator

The sharper sign measure is

```math
r_B=\int(-c_B(\ell))_+\,d\mu_C(\ell).
```

Retain each genuine reference loop independently with probability $`|c_B(\ell)|`$. Give a retained loop the sign $`\operatorname{sign}(c_B(\ell))`$, let $`S`$ be the product of its retained signs, and let $`n_{\rm abs}`$ be the total retained occupation. The resulting Poisson intensity is $`|c_B|\mu_C`$. Its ideal identities are

```math
\mathbb E S=d=e^{-2r_B},
\qquad
\mathbb E[Sf(n_{\rm abs})]=d\,\mathbb E_{P_B}f.
```

Both observables are bounded by one. Thus their ratio estimates the target expectation without computing a multiplicative normalizer. The finite-bit implementation is especially direct, but the required accuracy still depends on the lower bound for the signed denominator.

### 10.1 A certified upper bound for the negative-loop intensity

Supply a rational $`R_r\ge r_B`$. It can be replaced by $`\min\{R_r,M/4\}`$, because

```math
4(-c)_+\le(1-c)^2\quad(-1\le c\le1),
\qquad
r_B\le\Gamma/4\le M/4.
```

Compute a positive rational $`d_*`$ such that

```math
\tfrac12e^{-2R_r}\le d_*\le e^{-2R_r}\le d.
```

An upper approximation to $`e^{2R_r}`$ within a factor two, followed by inversion, provides it in polynomial bit time by Section 5's exponential-series construction. Its encoding length is polynomial in $`M`$ and that of $`R_r`$.

Two ways to obtain $`R_r`$ have different costs and strengths. A certified upper bound $`G\ge\Gamma`$ immediately gives $`R_r=G/4`$; this retains the $`e^\Gamma`$ guarantee but need not approximate $`r_B`$ closely. Alternatively use the absolutely convergent Fourier certificate

```math
r_B=\frac12\left[
\frac2\pi\Lambda_C
+\frac4\pi\sum_{k\ge1}\frac{(-1)^{k-1}}{4k^2-1}
\Lambda_{C\circ u^{2k}}
-\Lambda_B\right].
```

The tail after $`J`$ harmonics is at most $`\Lambda_C/[\pi(2J+1)]`$. If every log-determinant is evaluated to absolute error $`\xi`$, its combined contribution to the error is at most $`(1/2+2/\pi)\xi`$. These are absolute-error certificates, including when $`r_B=0`$.

For a requested additive certificate accuracy $`\eta_r`$, one conservative explicit choice is

```math
J=\max\{1,\lceil M/\eta_r\rceil\},
\qquad
\xi\le\eta_r/8,
```

and an approximation to $`1/\pi`$ of error at most $`\eta_r/[16(M+1)]`$. The sensitivity of the displayed expression to $`1/\pi`$ is at most $`2\Lambda_C\le M`$. Combining these bounds gives an interval of radius less than $`\eta_r`$ around the computed value. Its upper endpoint, intersected with the known range, gives $`r_B\le R_r\le r_B+2\eta_r`$. If a finite algorithm for $`\pi`$ is desired, the identity $`\pi=16\arctan(1/5)-4\arctan(1/239)`$ and the alternating power series yield certified intervals with a geometrically decreasing remainder.

The phase powers in the $`k`$th harmonic require no additional oracle assumption: if the contracted unit-phase approximation has error at most $`2a`$, then its $`2k`$th power has error at most $`4ka`$. Taking $`a\le\xi/(8JmA)`$ and using Section 6's trace-series cutoff evaluates every required log-determinant within $`\xi`$.

With ideal dense log-determinant arithmetic, this certificate uses $`O((J+1)m^3)`$ operations. The explicit finite-bit trace-series implementation adds its polynomial series and integer-arithmetic costs. The dependence $`J=O(1+\Lambda_C/\eta_r)`$, or the conservative $`O(1+M/\eta_r)`$ bound above, is genuinely inverse in $`\eta_r`$; it must not be replaced by logarithmic precision dependence. The ratio estimator needs only a **constant additive** certificate, for example $`\eta_r=1/8`$. This gives $`R_r\le r_B+1/4`$ and hence only a constant factor in $`e^{4R_r}`$. In contrast, implementing the direct prefactor $`e^{2r_B}`$ to control an $`O(\epsilon)`$ bias can require $`\eta_r=O(\epsilon)`$ and therefore an additional $`O(1/\epsilon)`$ Fourier truncation cost.

### 10.2 Stable thinning near a zero cosine

For one loop with cosine $`c`$, regard thinning and its sign as a single three-outcome categorical distribution:

```math
(\Pr[\mathrm{reject}],\Pr[+],\Pr[-])
=(1-|c|,c_+,(-c)_+).
```

Replacing $`c`$ by any $`\widetilde c\in[-1,1]`$ changes this distribution in total variation by at most $`|c-\widetilde c|`$. This follows from the Lipschitz bound for absolute value and

```math
|c_+-\widetilde c_+|+|(-c)_+-(-\widetilde c)_+|
=|c-\widetilde c|.
```

Accordingly, a sign error near $`c=0`$ does not require resolving the sign to relative accuracy: the probability of retaining that loop is small. Evaluating and then deterministically choosing a sign without controlling the retention decision would miss this protection.

On a reference record of at most $`F`$ edges, contracted phase approximations of accuracy $`a`$ change the joint thinning law by at most $`2Fa`$. The approximate three-outcome probabilities are exact rational numbers and are sampled by the finite-grid categorical rule. At most $`F/2`$ additional draws are needed.

### 10.3 Complete finite-bit ratio guarantee

Let

```math
t=\epsilon d_*/64.
```

Construct the reference record with $`\tau_L=\tau_N=t/4`$ and $`\eta_{\rm num}=t/8`$. Evaluate phases with $`a\le t/(16F)`$. Sample each approximate thinning/sign distribution with

```math
r_{\rm thin}=\left\lceil\log_2\frac{8F}{t}\right\rceil
```

random bits. The thinning approximation contributes at most $`t/8`$ in total variation and its categorical rounding contributes at most $`t/8`$. Including the reference approximation, the implemented signed-count record can therefore be coupled to the ideal absolute-cosine process with disagreement probability less than $`t`$.

Return $`(0,0)`$ on a reference failure. Otherwise form the two bounded observables

```math
Y=S\widetilde f(n_{\rm abs}),\qquad Z=S.
```

Evaluate $`f`$ to error at most $`t/2`$ and round it to a common dyadic grid to an additional error at most $`t/2`$, with clipping to $`[-1,1]`$. The implemented means then obey

```math
|\mathbb E Y-d\,\mathbb E_{P_B}f|\le3t,
\qquad
|\mathbb E Z-d|\le2t.
```

For

```math
L\ge\left\lceil2t^{-2}\log(4/\beta)\right\rceil
```

independent pairs, Hoeffding's inequality and a union bound ensure that both empirical means differ from their implemented expectations by at most $`t`$, except with probability $`\beta`$. On this event each total mean error is at most $`4t`$ and the empirical denominator is at least $`d-4t\ge15d/16`$.

Output zero if the empirical denominator is at most $`d_*/2`$; otherwise output the ratio of empirical means, clipped to $`[-1,1]`$. The exceptional branch is not taken on the preceding good event. Writing $`a=\mathbb E_{P_B}f`$, its error there is at most

```math
\frac{8t}{d-4t}\le\frac{2\epsilon}{15}<\epsilon.
```

This proves a bounded finite-bit estimator with

```math
L=O\left(e^{4R_r}\epsilon^{-2}\log(1/\beta)\right)
```

draws and polynomial preprocessing/per-draw bit cost. A constant-additive Fourier upper certificate gives the sharper $`O(e^{4r_B}\epsilon^{-2}\log(1/\beta))`$ dependence. The cutoffs have sizes

```math
\begin{aligned}
K&=O\bigl((M+1)[\log(m(M+1))+R_r+\log(1/\epsilon)]\bigr),\\
F&=O\bigl((M+1)[1+R_r+\log(1/\epsilon)]\bigr).
\end{aligned}
```

Unlike the direct weighted estimator, this implementation does not evaluate $`e^\Delta`$. Its numerical sampling and phase precision depends on $`R_r+\log(1/\epsilon)+\log F`$, together with the polynomial matrix/phase-certificate costs just stated. Dyadic sample values make the empirical sums polynomial in bit length and linear in the number of draws up to bit factors.

### 10.4 Finite-bit positive-part sampling

Let $`Q_+`$ be the occupation law of the Poisson process with intensity $`(c_B)_+\mu_C`$. Retain each reference loop with probability $`(c_B(\ell))_+`$ and discard its sign. This is a positive count sampler for every input in the model, with a target-dependent approximation guarantee.

For any requested algorithmic error $`0<\delta_{\rm alg}<1`$, choose

```math
\tau_L=\tau_N=\eta_{\rm num}=\delta_{\rm alg}/4
```

in the reference construction, compute contracted edge phases with $`a\le\delta_{\rm alg}/(32F)`$, and implement each retention Bernoulli draw with $`\lceil\log_2(8F/\delta_{\rm alg})\rceil`$ random bits. The positive-part map is 1-Lipschitz, so the combined error in retention probabilities is at most $`2Fa\le\delta_{\rm alg}/16`$. There are at most $`F/2`$ retention draws; their accumulated grid error is at most $`\delta_{\rm alg}/16`$. Map a reference failure to the vacuum. The resulting law satisfies

```math
d_{\rm TV}(\widehat Q_+,Q_+)\le\delta_{\rm alg}
```

and has polynomial bit cost in $`m,b,M,\log(1/\delta_{\rm alg})`$ and the phase-oracle costs. Neither a nonzero cosine gap nor a sign oracle is required.

To state the target error, let $`\nu_-`$ be the positive occupation measure of the negative loop part, of total mass $`r_B`$. The compensated signed-loop identity gives

```math
P_B=Q_+*\exp_*\left(r_B\delta_0-\nu_-\right).
```

The compensated measure has $`\ell^1`$ norm at most $`2r_B`$. Submultiplicativity of convolution and the norm-convergent exponential series imply

```math
d_{\rm TV}(P_B,Q_+)\le
\min\left\{1,\frac{e^{2r_B}-1}{2}\right\}.
```

Consequently the finite-bit full-distribution guarantee is

```math
d_{\rm TV}(P_B,\widehat Q_+)
\le\min\left\{1,\frac{e^{2r_B}-1}{2}\right\}
+\delta_{\rm alg}.
```

One may also cap this entire right-hand side at one. Running the sampler does not require evaluating $`r_B`$. Certifying a particular small structural error may require a sufficiently accurate upper bound, whose Fourier precision cost is the one stated in Section 10.1; it is separate from the algorithmic sampling error.

## 11. General complex-entry oracles and stability at zero entries

The polar input of Section 1 is convenient, but it is not required. Suppose instead that a deterministic oracle evaluates each entry of a symmetric complex matrix $`B`$ to absolute error $`2^{-p}`$ in polynomial time, and a rational $`0<\kappa\le1`$ certifies

```math
\||B|\|_2\le1-\kappa.
```

There is no promise that a nonzero entry has a minimum magnitude, nor an oracle deciding exact zero. The following reduction gives the same qualitative bit theorem with polynomial dependence on $`m,\kappa^{-1}`$ and preserves both $`r_B`$ and $`\Gamma`$ up to additive constants.

### 11.1 A rational-magnitude approximation

Set $`M_0=\lceil2m/\kappa\rceil`$. For a requested distribution tolerance $`0<\delta<1`$, choose a dyadic $`\zeta`$ satisfying

```math
\zeta\le
\min\left\{
\frac{\kappa}{4m},
\frac{\kappa}{40m^2},
\frac{\delta\kappa}{16m^2}\,2^{-2M_0}
\right\}.
```

Evaluate one rational $`v_{ij}=v_{ji}`$ within $`\zeta`$ of each entry. Choose a nonnegative rational $`c'_{ij}=c'_{ji}`$ within $`\zeta`$ of $`|v_{ij}|`$. If $`v_{ij}=0`$, choose $`c'_{ij}=0`$ and $`u'_{ij}=1`$; otherwise let $`u'_{ij}=v_{ij}/|v_{ij}|`$. Set $`B'_{ij}=c'_{ij}u'_{ij}`$ and $`C'=(c'_{ij})`$.

The square root and division needed to evaluate $`u'`$ to arbitrary absolute accuracy have polynomial bit cost: a nonzero rational $`v_{ij}`$ is separated from zero by its own explicitly represented denominator and numerator lengths. These lengths are bounded by the original oracle's output size. No lower bound on the original entry $`B_{ij}`$ is used.

Entrywise $`|B'-B|\le2\zeta`$. Also

```math
\|C'\|_2\le\||B|\|_2+2m\zeta\le1-\kappa/2,
\qquad N_{C'}\le M_0.
```

Thus $`B'`$ fits the first input model, with ordinary rational magnitudes and efficiently computable unit phases.

### 11.2 Distribution stability without phase division

The required continuity bound can be proved directly from the signed loop expansion. For two targets $`B,B'`$, let $`R_{ij}=\max\{|B_{ij}|,|B'_{ij}|\}`$ and $`E_{ij}=|B_{ij}-B'_{ij}|`$. Assume $`\|R\|_2=r<1`$. The sum of the absolute changes of all even-loop intensities is at most

```math
\begin{aligned}
\eta_{\rm loop}
&\le\sum_{k\ge1}\operatorname{Tr}(ER^{2k-1})\\
&\le\frac{m r}{1-r^2}\,\|E\|_2.
\end{aligned}
```

Indeed, expand a difference of two length-$`2k`$ edge products into $`2k`$ terms, each containing one changed edge and otherwise entries bounded by $`R`$. Summing over rooted closed walks yields $`2k\operatorname{Tr}(ER^{2k-1})`$; the loop factor $`1/(2k)`$ cancels it. Taking real parts cannot increase the bound.

Write the compensated signed measure on occupations as $`\nu_B=\int(\delta_{v(\ell)}-\delta_0)c_B(\ell)\,d\mu_{|B|}(\ell)`$. Equivalently its weights are alternating products of the entries of $`B`$, assigning zero to a product with a zero edge; this definition never divides by an original entry. Its convolution exponential is the physical law $`P_B=\exp_*\nu_B`$. The preceding bound implies

```math
\|\nu_B-\nu_{B'}\|_1\le2\eta_{\rm loop},
\qquad
\|\nu_B\|_1,\|\nu_{B'}\|_1\le2\Lambda_R.
```

The norm-convergent convolution series, or its integral derivative, therefore gives

```math
d_{\rm TV}(P_B,P_{B'})
\le \eta_{\rm loop}\,e^{2\Lambda_R}.
```

This is a conservative absolute-stability estimate. Its exponential constant costs precision bits, not an exponential number of samples. For the constructed approximation, $`R\le |B|+2\zeta\mathbf1\mathbf1^T`$, so $`r\le1-\kappa/2`$, $`2\Lambda_R\le N_R\le M_0`$, and

```math
\eta_{\rm loop}\le\frac{4m^2\zeta}{\kappa},
\qquad
d_{\rm TV}(P_B,P_{B'})
\le\frac{4m^2\zeta}{\kappa}e^{M_0}\le\delta/4.
```

The last inequality uses $`e^{M_0}\le2^{2M_0}`$ and the chosen $`\zeta`$. Hence only

```math
O\bigl(M_0+\log(m/(\kappa\delta))\bigr)
```

entry-accuracy bits are needed. For $`|f|\le1`$, replacing the target changes its expectation by at most twice this total-variation bound.

### 11.3 Stability of both sign certificates

Since the negative-part function is 1-Lipschitz, the same absolute loop-weight bound immediately gives

```math
|r_B-r_{B'}|\le\eta_{\rm loop}\le\frac{4m^2\zeta}{\kappa}\le\frac1{10}.
```

Thus a constant-additive upper certificate for $`r_{B'}`$ preserves the $`e^{4r_B}`$ repetition count up to a universal constant. This argument remains valid when arbitrarily small entries change phase or pass through zero.

Define $`d(z)=z^2/|z|`$ for $`z\ne0`$ and $`d(0)=0`$. This map is globally Lipschitz with constant three. For example, when $`|z|\ge|w|`$, split the difference using the common denominator $`|z|`$ and use $`||z|-|w||\le|z-w|`$ to obtain $`|d(z)-d(w)|\le3|z-w|`$. Thus the phase-doubled matrix remains stable even when the phase of an individual tiny entry is ill conditioned.

Combining this inequality with the log-determinant Lipschitz bound in Section 6 gives, for the constructed matrices,

```math
|\Gamma(B')-\Gamma(B)|
\le\frac{20m^2\zeta}{\kappa}\le\frac12.
```

Accordingly the estimator for $`B'`$ uses $`O(e^{\Gamma(B)}\epsilon^{-2}\log(1/\beta))`$ repetitions up to a universal constant factor. Run it at a slightly smaller error tolerance and select $`\delta`$ proportional to $`\epsilon`$ in the reduction. The preceding distribution bound then supplies the desired estimate for the original $`B`$.

Using the ratio estimator, the resulting complete running time is polynomial in $`m,\kappa^{-1},\epsilon^{-1},\log(1/\beta)`$, the entry/observable oracle costs, and $`e^{4r_B}`$, with the repetition count linear in $`e^{4r_B}`$ up to the displayed accuracy and confidence factors. A supplied inverse-polynomial gap and $`r_B=O(\log m)`$ therefore suffice. The direct weighted estimator alternatively retains its $`e^{\Gamma(B)}`$ guarantee. The reduction is deliberately conservative about numerical conditioning; it does not require exact zero detection, relative entry accuracy, or a minimum nonzero matrix element.

## References and attribution

The Gaussian probability law and the loop representation are inherited inputs to this implementation analysis. The construction of loop soups from ordered vertices and their first-return excursions is established probability theory:

- Gregory F. Lawler, *Topics in loop measures and the loop-erased walk*, Sections 5.3–5.4, lecture notes dated 22 March 2017: <https://www.math.uchicago.edu/~lawler/fall16.pdf>.
- Yves Le Jan, *Markov loops and renormalization*, Annals of Probability **38**, 1280–1319 (2010): <https://arxiv.org/abs/0802.2478>.

The elementary rational algorithms, explicit error budgets, exponential cutoffs, and stability reduction above provide the finite-bit implementation guarantee. They should not be conflated with the ideal arithmetic cost of a different loop-sampling implementation.
