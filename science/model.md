# Gaussian photon counts and the classical loop representation

## The question

When can interference in a pure Gaussian optical state be handled by a classical random process with a controlled cost? Photon probabilities contain squared hafnians, but their joint generating function is a determinant. Expanding its logarithm turns the distribution into an ensemble of closed walks. Positive walk weights give an ordinary probability distribution; complex phases introduce a signed measure whose negative part can be quantified.

The project separates three tasks throughout:

| Task | Output | Required guarantee |
|---|---|---|
| Full-distribution sampling | A photon-count vector | Total variation distance from the entire target law |
| Bounded-observable estimation | An estimate of a specified expectation | Additive error and failure probability |
| Probability evaluation | A probability of a specified count pattern | A separately specified additive or relative error |

A signed estimator solves the second task. It does not by itself solve the first, and an additive estimate of an exponentially small probability need not be useful as a relative estimate.

## Reading route and result map

1. This page defines the state, the loop measure, and the resource parameters.
2. [Direct sampler](direct-sampler.md) proves the positive sampling construction and its arithmetic cost.
3. [Phase estimators](estimators.md) proves the signed estimators, negative-mass certificates, and positive-part sampling bound.
4. [Topology](topology.md) gives finite graph bounds and an explicit expander family with extensive photon number.
5. [Finite precision](finite-precision.md) gives a separate finite-bit construction and accounts for tails and numerical error.
6. [Interpretation](interpretation.md) explains the physical and computational scope. [Related work](../references/RELATED_WORK.md) distinguishes inherited ingredients from the results assembled here.

## State and photon law

Let $`B=B^{\mathsf T}\in\mathbb C^{m\times m}`$ and $`\|B\|_2<1`$. The zero-displacement pure Gaussian state is

```math
|\psi_B\rangle
=\det(I-B^\dagger B)^{1/4}
\exp\!\left(\frac12\sum_{i,j=1}^m B_{ij}a_i^\dagger a_j^\dagger\right)|0\rangle.
```

The Takagi factorization $`B=U\operatorname{diag}(s_j)U^{\mathsf T}`$, with $`0\le s_j<1`$, identifies independent single-mode squeezers followed by a passive interferometer. The squeeze parameters satisfy $`s_j=\tanh r_j`$. No displacement, adaptive measurement, or input non-Gaussian state is included in this model. The standard Gaussian-state and hafnian formulas are reviewed by [Weedbrook et al.](https://doi.org/10.1103/RevModPhys.84.621) and used in [Hamilton et al.](https://doi.org/10.1103/PhysRevLett.119.170501).

For $`n\in\mathbb N^m`$, let $`B[n]`$ repeat row and column $`i`$ exactly $`n_i`$ times. The hafnian sums products over perfect pairings of those labeled positions; the empty hafnian is one. Odd photon totals have probability zero, and

```math
P_B(n)=\sqrt{\det(I-B^\dagger B)}\,
\frac{|\operatorname{haf}(B[n])|^2}{\prod_i n_i!}.
```

The factorial is essential when photons collide in a mode. This project samples the complete count law, including all collisions and all even photon totals; it does not condition on a fixed total.

### Generating function

Write $`Z=\operatorname{diag}(z_1,\ldots,z_m)`$, initially with $`0\le z_i\le1`$. Multiplying each occupation amplitude by $`\prod_i z_i^{n_i/2}`$ replaces $`B`$ by $`Z^{1/2}BZ^{1/2}`$ in the unnormalized state. Its squared norm is the inverse square root of the corresponding Gaussian determinant. Cyclic invariance of the determinant therefore gives

```math
F_B(z):=\mathbb E_B\prod_i z_i^{n_i}
=\sqrt{\frac{\det(I-B^\dagger B)}
{\det(I-B^\dagger ZBZ)}}.
```

The root is the analytic branch that is positive at $`z=0`$; the identity also holds as a convergent power series around the closed unit polydisc when the reference matrix below is stable. In particular, for scalar $`z_i=s`$,

```math
\overline N_B=\mathbb E_B\sum_i n_i
=\operatorname{Tr}\!\left[B^\dagger B(I-B^\dagger B)^{-1}\right]
=\sum_j\frac{s_j^2}{1-s_j^2}.
```

## The stable magnitude reference

Set $`C_{ij}=|B_{ij}|`$. All loop algorithms in this repository assume

```math
C=C^{\mathsf T}\ge0\quad\text{entrywise},\qquad \rho:=\|C\|_2<1.
```

This is stronger than the physical condition $`\|B\|_2<1`$. A physical Gaussian state can fail this reference condition; no conclusion about its classical difficulty follows from that failure. Entrywise domination gives $`\|B\|_2\le\|C\|_2`$.

For any complex symmetric matrix $`X`$ dominated entrywise by $`C`$, define

```math
\Lambda_X=-\frac12\log\det(I-X^\dagger X),\qquad
\overline N_X=\operatorname{Tr}\!\left[X^\dagger X(I-X^\dagger X)^{-1}\right].
```

The reference preparation and path-generation cost depend on $`\overline N_C`$, which can exceed $`\overline N_B`$. A small target photon number alone is not the runtime promise.

## Lifted walks

Use two copies of the mode set and define

```math
T_C=\begin{pmatrix}0&C\\C&0\end{pmatrix},\qquad
T_B=\begin{pmatrix}0&B\\B^\dagger&0\end{pmatrix},\qquad
\widehat Z=\operatorname{diag}(Z,Z).
```

Both lifted matrices have norm below one, and $`T_B`$ is Hermitian. A lifted vertex projects to its original mode. A closed walk $`w=(v_0,\ldots,v_{L-1},v_L=v_0)`$ has reference edge product $`\prod_{j=0}^{L-1}(T_C)_{v_jv_{j+1}}`$, occupation $`n(w)`$, and length $`L=|n(w)|`$. The lifted graph is bipartite, so every nonzero loop length is even. The graph may have edges between copies of the same mode when $`B_{ii}\ne0`$.

The loop measure at intensity one half is specified, for any root-invariant function $`f`$, by

```math
\int f(w)\,d\mu_C(w)
=\frac12\sum_{L\ge1}\frac1L
\sum_{v_0,\ldots,v_{L-1}}
\left(\prod_{j=0}^{L-1}(T_C)_{v_jv_{j+1}}\right)f(w).
```

The sum is over ordered rooted walks. The $`1/L`$ factor automatically gives the correct weight to periodic walks; replacing it by an informal count of distinct geometric cycles would be incorrect. Equivalently this is a measure on unrooted oriented loops with their inherited multiplicities. It includes backtracking and repeated vertices, not only simple graph cycles.

Its total mass and expected occupation are

```math
\mu_C(\mathcal L)=\frac12\sum_{L\ge1}\frac{\operatorname{Tr}T_C^L}{L}
=\Lambda_C,
\qquad
\int L\,d\mu_C=\frac12\sum_{L\ge1}\operatorname{Tr}T_C^L
=\overline N_C.
```

The trace-log identity gives

```math
\log F_C(z)
=\sum_{k\ge1}\frac{\operatorname{Tr}(CZ)^{2k}-\operatorname{Tr}C^{2k}}{2k}
=\int\bigl(z^{n(w)}-1\bigr)\,d\mu_C(w).
```

Thus a Poisson collection of loops with intensity $`\mu_C`$, projected to the total occupation vector, has exactly law $`P_C`$. This determinant/loop correspondence is an inherited result of Markov-loop theory; see [Le Jan](https://doi.org/10.1214/09-AOP509) and [Lawler, Sections 5.2–5.4](https://www.math.uchicago.edu/~lawler/fall16.pdf). The photon interpretation and the phase-dependent resource accounting are developed in the following pages.

## Holonomy and phase resources

On each nonzero lifted edge write $`(T_B)_{uv}=(T_C)_{uv}e^{i\alpha_{uv}}`$. Hermiticity gives $`\alpha_{vu}=-\alpha_{uv}`$ modulo $`2\pi`$. A loop has holonomy $`e^{i\phi(w)}=\prod_{(u,v)\in w}e^{i\alpha_{uv}}`$. Its reverse has the conjugate holonomy and the same occupation and reference weight. Pairing orientations replaces the phase by $`c(w)=\cos\phi(w)`$:

```math
\log F_B(z)=\int c(w)\bigl(z^{n(w)}-1\bigr)\,d\mu_C(w).
```

Absolute convergence follows from domination by the finite positive measure $`\mu_C`$. Define

```math
\Delta=\int(1-c)\,d\mu_C=\Lambda_C-\Lambda_B,\qquad
\Gamma=\int(1-c)^2\,d\mu_C,\qquad
r=\int(-c)_+\,d\mu_C.
```

These are nonnegative properties of the chosen loop representation. Their operational meanings are proved in [Phase estimators](estimators.md): $`e^\Gamma`$ is a cosine-estimator second moment, $`e^{4r}`$ is the second moment after absolute-cosine thinning, and small $`r`$ certifies full-distribution approximation by a positive loop process. They are sufficient algorithmic resources, not lower bounds against unrestricted classical algorithms.

For a diagonal unitary $`D`$, replacing $`B`$ by $`DBD`$ only changes local photon phases and leaves the count law invariant. If $`B=DCD`$, all lifted loop holonomies are one and $`\Delta=\Gamma=r=0`$. Large entry phases alone therefore do not measure interference difficulty; closed-walk phases are the relevant quantities here.
