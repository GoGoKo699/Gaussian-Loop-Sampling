# Phase interference, graph cycles, and full-distribution approximation

This note proves when the photon-count distribution of a complex squeezed-state matrix can be sampled approximately by a positive loop process. There are two mechanisms. Small residual edge phases make a finite collection of short loops positive. Large graph girth makes short loops exactly insensitive to arbitrary edge phases. An Ihara–Bass transformation removes immediate backtracking from the latter argument; a Ramanujan spectral bound then gives an explicit connected family with linear support treewidth, extensive photon number, and a vanishing total-variation error at fixed squeezing scale.

For the degree-six bipartite Lubotzky–Phillips–Sarnak family, the resulting sufficient interval is

```math
0<t<t_*=
\frac{6\,5^{-1/4}}{5+5^{-1/2}}
=0.7366044601\ldots,
\qquad |B|=\frac{t}{6}A_G.
```

The edge phases of $`B`$ are arbitrary. The interval is a sufficient condition for the approximation proved here; its endpoint is not a computational-complexity transition. All distributions below are **unconditioned photon counts, including collisions**, for pure Gaussian states with zero displacement. Bounded-observable reweighting is a separate task, treated in [the estimator note](estimators.md).

## 1. Conventions and the signed loop representation

The [model note](model.md) defines the state, hafnian count law, probability-generating function, and rooted loop measure, and proves their identities. Here $`P_B`$ denotes that full count law and $`F_B(z)=\sum_nP_B(n)z^n`$ its generating function. Let $`B=B^T\in\mathbb C^{m\times m}`$, put $`C=|B|`$ entrywise, and assume $`\|C\|_2<1`$. This is stronger than physicality alone, and implies $`\|B\|_2\le\|C\|_2<1`$. Use

```math
\Lambda_B=-\frac12\log\det(I-B^\dagger B),
\qquad
\Delta=\Lambda_C-\Lambda_B.
\tag{1}
```

All logarithms in this note are natural logarithms except those with an explicit base. The Hermitian doubled matrices are

```math
T_B=\begin{pmatrix}0&B\\B^\dagger&0\end{pmatrix},
\qquad
T_C=\begin{pmatrix}0&C\\C&0\end{pmatrix}.
```

They satisfy $`|T_B|=T_C`$. For an even rooted loop $`w=(v_0,\ldots,v_{L-1},v_L=v_0)`$, the model's positive intensity and its magnetic phase are

```math
\mu_C(w)=\frac{1}{2L}
   \prod_{j=0}^{L-1}(T_C)_{v_jv_{j+1}},
\qquad
\chi_B(w)=\prod_{j=0}^{L-1}
   \frac{(T_B)_{v_jv_{j+1}}}{(T_C)_{v_jv_{j+1}}}.
\tag{2}
```

The count $`n_i(w)`$ is the number of visits to the two copies of mode $`i`$, excluding the duplicated terminal vertex. Reversal conjugates $`\chi_B(w)`$ and preserves its intensity and count vector. Writing $`\chi_B(w)=e^{i\Phi(w)}`$, the model's trace-log identities take the form needed below:

```math
\begin{aligned}
\log F_B(z)
 &=\sum_w\mu_C(w)\cos\Phi(w)\,[z^{n(w)}-1],\\
\Lambda_C
 &=\sum_w\mu_C(w)
   =\sum_{k\ge1}\frac{\operatorname{Tr}C^{2k}}{2k},\\
\Delta
 &=\sum_w\mu_C(w)\,[1-\cos\Phi(w)]\ge0.
\end{aligned}
\tag{3}
```

All these series converge absolutely. Independent Poisson loops with intensities (2), projected onto the original modes, reproduce $`P_C`$ exactly. The underlying loop-measure and signed-Poisson constructions are established mathematics; see [Lawler](#references).

### Optical gauge and magnetic gauge

A diagonal optical phase change has the form $`B\mapsto UBU`$, with $`U`$ diagonal unitary. It leaves photon-count probabilities invariant. On the doubled graph it is the ordinary unitary similarity

```math
T_{UBU}=V T_B V^\dagger,
\qquad V=\operatorname{diag}(U,\overline U).
```

Consequently, every loop holonomy is unchanged. In particular, $`B=UCU`$ implies $`P_B=P_C`$ and $`\Delta=0`$. The doubled-graph construction is necessary when interpreting phases: a symmetric optical edge has the same phase in $`B_{ij}`$ and $`B_{ji}`$, while reversed edges of the Hermitian lift have conjugate phases.

## 2. A determinant certificate for full-distribution sampling

**Proposition 1.** Under the assumptions above,

```math
P_B(n)\le e^\Delta P_C(n)
\quad\text{for every }n,
\qquad
\|P_B-P_C\|_{\mathrm{TV}}
\le1-e^{-\Delta}\le\Delta.
\tag{4}
```

**Proof.** Expand the hafnian over perfect matchings. The triangle inequality gives $`|\operatorname{haf}(B_n)|\le\operatorname{haf}(C_n)`$. The ratio of the state-normalization factors is $`e^\Delta`$, proving the pointwise bound, including configurations of zero reference probability. Since

```math
\min\{P_B(n),P_C(n)\}\ge e^{-\Delta}P_B(n),
```

summing and using $`\|P-Q\|_{\mathrm{TV}}=1-\sum_n\min(P(n),Q(n))`$ proves (4). $`\square`$

Thus two ordinary log determinants certify a sampler: draw from $`P_C`$ and use the result as an approximation to $`P_B`$. This is a statement about the entire output distribution, including global events. It does not rely on enumerating output probabilities.

If an implemented reference sampler has error $`\varepsilon_{\mathrm{ref}}`$ relative to $`P_C`$, its error relative to the target is at most

```math
\varepsilon_{\mathrm{ref}}+1-e^{-\Delta}.
\tag{5}
```

The first term includes finite precision and any deliberate output truncation. The second is the structural approximation error proved here.

## 3. Positive short-loop thinning

Small residual optical phases give a different positive sampler, which can be more accurate than replacing all phases by zero.

**Theorem 2.** Suppose an optical gauge makes every supported entry satisfy $`B_{ij}=C_{ij}e^{i\theta_{ij}}`$ with real $`|\theta_{ij}|\le\theta`$. Choose a positive integer $`K`$ with $`2K\theta\le\pi/2`$. Generate the $`C`$ loop soup restricted to lengths at most $`2K`$, retaining each complete loop $`w`$ independently with probability $`\cos\Phi(w)`$. Let $`Q_K`$ be its projected count law. Then

```math
\begin{aligned}
\|P_B-Q_K\|_{\mathrm{TV}}
 &\le\frac12\left(e^{2\Lambda_{>2K}}-1\right),\\
\Lambda_{>2K}
 &=\sum_{k>K}\frac{\operatorname{Tr}C^{2k}}{2k}\\
 &\le\frac{\overline N_C}{2K+2},\\
\Lambda_{>2K}
 &\le\frac{m r^{2K+2}}{(2K+2)(1-r^2)}
 \quad\text{if }\|C\|_2\le r<1,
\end{aligned}
\tag{6}
```

where $`\overline N_C=\operatorname{Tr}[C^2(I-C^2)^{-1}]`$.

**Proof.** For a lifted loop of length $`L`$, sum the signed real edge phases before reduction modulo $`2\pi`$. The absolute value of this sum is at most $`L\theta`$. Hence the retention probability is nonnegative for every retained length.

For the error estimate, work in the convolution Banach algebra $`\ell^1(\mathbb N_0^m)`$. Push the long-loop signed intensity forward to count space:

```math
\nu_{\mathrm{tail}}
 =\sum_{|w|>2K}\mu_C(w)\cos\Phi(w)\,\delta_{n(w)},
\qquad
\kappa=\nu_{\mathrm{tail}}
       -\nu_{\mathrm{tail}}(\mathbb N_0^m)\delta_0.
```

Its norm satisfies $`\|\kappa\|_1\le2\Lambda_{>2K}`$. Equation (3), together with uniqueness of absolutely convergent generating functions, gives the exact signed-measure identity

```math
P_B=Q_K*\exp_*(\kappa),
\qquad
\exp_*(\kappa)=\sum_{j\ge0}\frac{\kappa^{*j}}{j!}.
```

The second factor has total mass one and need not be positive. Absolute convergence and submultiplicativity give

```math
\|\exp_*(\kappa)-\delta_0\|_1
\le e^{\|\kappa\|_1}-1.
```

Convolution with the probability measure $`Q_K`$ contracts the $`\ell^1`$ norm, proving the first inequality of (6). The remaining bounds follow from $`2k\ge2K+2`$ and $`\operatorname{Tr}C^{2k}\le m r^{2k}`$. $`\square`$

To target error $`\varepsilon`$, it suffices to impose both

```math
\Lambda_{>2K}\le\frac12\log(1+2\varepsilon),
\qquad
\theta\le\frac{\pi}{4K}.
\tag{7}
```

At fixed $`r<1`$, the first condition requires only $`K=O(\log(m/\varepsilon))`$. Thus residual phases of order $`1/\log(m/\varepsilon)`$ already allow full-distribution approximation. A fixed nonzero phase allowance places an upper bound on admissible $`K`$; (7) does not promise arbitrary accuracy for every fixed input.

Thinning acts on each **complete Poisson loop**. In an excursion implementation, the logarithmic-series cluster must be assembled before computing its phase and applying the length cutoff. Independently thinning constituent excursions defines a different process.

## 4. Girth eliminates short-loop phase dependence

Let $`G`$ be a finite simple $`d`$-regular graph on $`m`$ vertices, with girth $`g`$, and put

```math
C=\frac{t}{d}A_G,
\qquad 0<t<1,
\qquad B_{ij}=C_{ij}e^{i\theta_{ij}},
\quad\theta_{ij}=\theta_{ji}.
\tag{8}
```

There is no restriction on these phases. The doubled support graph is a two-sheeted cover of $`G`$, with $`2m`$ vertices and $`md`$ undirected edges. Its girth is at least $`g`$. Every closed walk shorter than $`g`$ reduces to the empty walk by deleting immediate reversals. On a Hermitian edge each reversal contributes $`\eta_e\eta_{\bar e}=1`$, so such a loop has trivial holonomy.

Set $`h=2\lceil g/2\rceil`$. Equations (3) and (4) give the elementary estimate

```math
\Delta\le
2\sum_{k\ge h/2}\frac{\operatorname{Tr}C^{2k}}{2k}
\le\frac{2m t^h}{h(1-t^2)}.
\tag{9}
```

This bound sums walks with arbitrarily many immediate reversals. The next sections sum those reversals exactly and improve the effective parameter controlling phase-sensitive cycles.

## 5. Magnetic Ihara–Bass identity with explicit weights

The determinant identity used here is classical, including twists by unitary local systems. [Bass](#references) proves the underlying formula; [Deitmar, Section 4](#references) treats weighted graphs with local systems. The elementary proof below fixes the edge orientation, weighting, and exponent conventions needed for the application.

### 5.1 Unequal magnitudes

Let $`H`$ be the Hermitian adjacency matrix of any finite simple graph, with

```math
H_{xy}=a_{xy}\eta_{xy},
\qquad a_{xy}=a_{yx}>0,
\qquad |\eta_{xy}|=1,
\qquad\eta_{yx}=\overline{\eta_{xy}}.
```

For a directed edge $`e=(x,y)`$, write $`\bar e=(y,x)`$, $`a_e=a_{xy}`$, and $`H_e=H_{xy}`$. Index the nonbacktracking matrix by directed edges and put the weight on the **next** edge:

```math
N_H[e,f]
=\mathbf1\{\operatorname{head}(e)=\operatorname{tail}(f),\ f\ne\bar e\}\,H_f.
\tag{10}
```

For $`1-u^2a_e^2\ne0`$, define a vertex matrix $`M_H(u)`$ by

```math
\begin{aligned}
[M_H(u)]_{xx}
 &=1+\sum_{y\sim x}\frac{u^2a_{xy}^2}{1-u^2a_{xy}^2},\\
[M_H(u)]_{xy}
 &=-\frac{uH_{xy}}{1-u^2a_{xy}^2}
 \quad(x\sim y),
\end{aligned}
```

with zero off-diagonal entries otherwise. Then

```math
\boxed{
\det(I-uN_H)
=\prod_{\{x,y\}\in E}(1-u^2a_{xy}^2)\,
  \det M_H(u).
}
\tag{11}
```

**Proof.** Let $`S[e,v]=\mathbf1\{\operatorname{head}(e)=v\}`$, $`R[v,f]=\mathbf1\{\operatorname{tail}(f)=v\}H_f`$, and $`J[e,f]=\mathbf1\{f=\bar e\}H_f`$. Then $`N_H=SR-J`$. Each undirected-edge block of $`I+uJ`$ has determinant $`1-u^2a_e^2`$. The matrix determinant lemma yields

```math
\begin{aligned}
\det(I-uN_H)
 &=\det(I+uJ-uSR)\\
 &=\det(I+uJ)
   \det\!\left[I-uR(I+uJ)^{-1}S\right].
\end{aligned}
```

Inverting the two-by-two edge blocks gives exactly $`M_H(u)`$ as the second factor. This proves (11); clearing denominators extends it as a polynomial identity. $`\square`$

### 5.2 Unit phases on a regular support

If all $`a_e=1`$, the vertex matrix is

```math
M_H(u)=\frac{I-uH+u^2(D-I)}{1-u^2},
```

where $`D`$ is the degree matrix. Consequently,

```math
\boxed{
\det(I-uN_H)
=(1-u^2)^{|E|-|V|}
 \det\!\left[I-uH+u^2(D-I)\right].
}
\tag{12}
```

For a $`d`$-regular graph, write $`p=d-1`$. On the doubled graph of (8), define the unit-magnitude magnetic adjacency $`H_\phi`$ by $`T_B=(t/d)H_\phi`$; write $`H_0`$ for the adjacency with all phases removed. Both have the same $`|V|=2m`$ and $`|E|=md`$.

Choose the smaller positive solution of

```math
s=\frac td=\frac{u}{1+p u^2},
\qquad
u=\frac{2s}{1+\sqrt{1-4ps^2}},
\qquad q=pu.
\tag{13}
```

For $`d\ge3`$ and $`0<t<1`$, one has $`0<q<1`$ and

```math
t=\frac{dq}{p+q^2},
\qquad q<t.
```

The scalar $`1+p u^2`$ and the factor $`(1-u^2)^{|E|-|V|}`$ in (12) do not depend on the phases. Their cancellation in the determinant ratio gives the exact identity

```math
\boxed{
\Delta
=\frac12\sum_{L\ge1}\frac{u^L}{L}
 \left(\operatorname{Tr}N_0^L-
       \operatorname{Tr}N_\phi^L\right).
}
\tag{14}
```

Indeed, $`I-uH_\phi+pu^2I=(1+pu^2)(I-sH_\phi)`$ and $`\det(I-sH_\phi)=\det(I-B^\dagger B)`$. Also $`\rho(N_\phi)\le p`$, since its entrywise absolute value is $`N_0`$ with constant row sum $`p`$. Hence $`q<1`$ makes the trace series absolutely convergent. The logarithms in the ratio are the real branches continuous from $`u=0`$; their determinants are positive by (12) and $`I-sH_\phi\succ0`$.

The traces in (14) count cyclically nonbacktracking closed walks, including the prohibition on reversal across the closing edge. Pairing each with its reverse shows

```math
0\le\operatorname{Tr}N_0^L-\operatorname{Tr}N_\phi^L
\le2\operatorname{Tr}N_0^L.
\tag{15}
```

Both traces vanish for odd $`L`$, and there are no such walks of length below $`g`$. The elementary estimate $`\operatorname{Tr}N_0^L\le2md\,p^L`$ therefore gives

```math
\Delta\le A_{m,d,h}(t)
:=\frac{2md\,q^h}{h(1-q^2)},
\qquad h=2\lceil g/2\rceil.
\tag{16}
```

Combining (4), (9), and (16) gives usable finite-instance certificates; the smaller available upper bound on $`\Delta`$ can always be selected. In particular, when $`g\ge c\log m-O(1)`$ and $`q<e^{-1/c}`$, (16) tends to zero as an inverse power of $`m`$, up to a logarithmic factor.

The general unequal-weight identity (11) does **not** give the scalar substitution (13) for arbitrary edge magnitudes. The regular high-girth bounds (16) and those below assume the uniform magnitude $`t/d`$ in (8).

## 6. Ramanujan spectra strengthen the cycle bound

Suppose now that $`G`$ is connected, bipartite, and $`d`$-regular, with all adjacency eigenvalues other than the simple eigenvalues $`\pm d`$ in $`[-2\sqrt p,2\sqrt p]`$, where $`p=d-1>1`$. Such a graph is Ramanujan. Its doubled graph is the disjoint union of two copies of $`G`$.

The following spectral count is a direct consequence of (12); the equivalent Chebyshev trace identity is also given by [Rangarajan, Lemma 4](#references). The edge matrix $`N_0`$ has $`2md`$ eigenvalues, counted with algebraic multiplicity:

| Source | Eigenvalues | Multiplicity |
|---|---|---:|
| The two copies of adjacency eigenvalue $`d`$ | $`p,1`$ | Two of each |
| The two copies of adjacency eigenvalue $`-d`$ | $`-p,-1`$ | Two of each |
| The other $`2m-4`$ adjacency eigenvalues $`\lambda`$ | Roots of $`x^2-\lambda x+p=0`$ | $`4m-8`$ roots in total |
| The Bass prefactor | $`1,-1`$ | $`m(d-2)`$ additional copies of each |

Every root in the third row has modulus $`\sqrt p`$, including repeated roots at the endpoints. Possible Jordan blocks do not change power traces. Thus, for even positive $`L`$,

```math
0\le\operatorname{Tr}N_0^L
\le4p^L+(4m-8)p^{L/2}+2m(d-2)+4.
\tag{17}
```

The four multiplicities sum to $`2md`$, as required. The terms $`2m(d-2)+4`$ and $`4p^L`$ retain the contributions of both disconnected components of the doubled graph.

Substituting (17) into (14)–(15), summing only even lengths at least $`h`$, and replacing $`4m-8`$ by $`4m`$ gives the explicit bound

```math
\boxed{
\begin{aligned}
\Delta\le R_{m,d,h}(t):=\frac1h\Bigg[
 &\frac{4q^h}{1-q^2}
 +\frac{4m(q/\sqrt p)^h}{1-q^2/p}\\
 &+\frac{[2m(d-2)+4](q/p)^h}{1-q^2/p^2}
\Bigg].
\end{aligned}
}
\tag{18}
```

This estimate is uniform over all symmetric choices of edge phases. It bounds the unsigned nonbacktracking cycle count; no spectral assumption on the phased matrix is imposed beyond $`|B|=tA_G/d`$.

## 7. An explicit scalable family

### 7.1 Graph construction and established graph properties

Let $`\ell`$ range over primes congruent to $`13`$ modulo $`20`$, and take the bipartite LPS graph $`G_\ell=X^{5,\ell}`$. These primes satisfy $`\ell\equiv1\pmod4`$ and $`(5/\ell)=-1`$. The graph is a connected six-regular Cayley graph on $`\operatorname{PGL}_2(\mathbb F_\ell)`$, of order

```math
m=\ell(\ell^2-1),
\qquad
g\ge4\log_5\ell-\log_5 4,
\qquad
|\lambda|\le2\sqrt5\quad(\lambda\ne\pm6).
\tag{19}
```

These are the LPS construction and Theorems 3.4 and 4.1 of [Lubotzky–Phillips–Sarnak](#references). Infinitely many such primes exist by Dirichlet's theorem.

For an explicit description, choose $`\iota\in\mathbb F_\ell`$ with $`\iota^2=-1`$. Use the following six matrices, considered modulo nonzero scalar multiples, as Cayley generators:

```math
\begin{pmatrix}1\pm2\iota&0\\0&1\mp2\iota\end{pmatrix},
\qquad
\begin{pmatrix}1&\pm2\\\mp2&1\end{pmatrix},
\qquad
\begin{pmatrix}1&\pm2\iota\\\pm2\iota&1\end{pmatrix}.
\tag{20}
```

They arise from the six integer quaternions with real part one and one imaginary coordinate equal to $`\pm2`$. Each determinant is five, and the generating set is closed under inversion in the projective group. Enumerating the group and multiplying by these six generators constructs an explicit sparse adjacency list in polynomial time in $`m`$; there is no graph-access oracle assumption.

### 7.2 Fixed-parameter accuracy

Define the even integer guaranteed by the paper's girth bound,

```math
h_\ell=2\left\lceil
\frac{4\log_5\ell-\log_5 4}{2}
\right\rceil.
```

For any $`0<t<1`$ and any symmetric edge phases on $`G_\ell`$, sampling $`P_C`$ obeys the finite-size certificate

```math
\|P_B-P_C\|_{\mathrm{TV}}
\le1-\exp[-R_{m,6,h_\ell}(t)].
\tag{21}
```

An independently established larger girth can replace $`h_\ell`$. Equation (21) uses only the proved lower bound (19).

Since $`m<\ell^3`$, one has $`h_\ell\ge(4/3)\log_5m-\log_5 4`$. For fixed $`t`$, set

```math
a(t)=-\frac13-\frac{4\log q(t)}{3\log5}.
\tag{22}
```

The three terms in (18) are respectively

```math
O_t\!\left(\frac{m^{-a(t)-1/3}}{\log m}\right),
\qquad
O_t\!\left(\frac{m^{-a(t)}}{\log m}\right),
\qquad
O_t\!\left(\frac{m^{-a(t)-2/3}}{\log m}\right).
```

Consequently, for every fixed $`t<t_*`$, where $`t_*`$ is the value in the introduction, $`a(t)>0`$ and

```math
\boxed{
\sup_{\{\theta_{ij}=\theta_{ji}\}}
\|P_B-P_C\|_{\mathrm{TV}}
=O_t\!\left(\frac{m^{-a(t)}}{\log m}\right).
}
\tag{23}
```

At $`t=t_*`$, the same bound is $`O(1/\log m)`$. The inverse-power statement requires the strict inequality.

For example, $`t=1/2`$ gives $`q=0.4322356371\ldots`$ and $`a(t)=0.3615547196\ldots`$. The following values are conservative certificates from (18) and (21), with the final column rounded upward; they are not empirical distances:

| Prime $`\ell`$ | Modes $`m`$ | Guaranteed $`h_\ell`$ | TV upper bound at $`t=1/2`$ |
|---:|---:|---:|---:|
| 13 | 2,184 | 6 | 0.081927 |
| 53 | 148,824 | 10 | 0.004609 |
| 113 | 1,442,784 | 12 | 0.001377 |

For a fixed graph and fixed $`t`$, (21) supplies a particular error guarantee. More reference samples do not reduce that structural bias. Arbitrarily small prescribed error requires checking (21) for the chosen family size, changing parameters, or using an additional target-dependent method. Numerical sampling error from the reference algorithm is added as in (5).

### 7.3 Extensive photons at bounded squeezing

Let $`\sigma_j`$ be the singular values of $`B`$. The total photon mean is

```math
\overline N_B
=\sum_j\frac{\sigma_j^2}{1-\sigma_j^2}.
```

For (8), $`\|B\|_2\le t`$ and

```math
\sum_j\sigma_j^2
=\|B\|_F^2
=\frac{mt^2}{d}.
```

It follows, for both the target and the positive reference, that

```math
\boxed{
\frac{mt^2}{d}
\le\overline N_B,\overline N_C
\le\frac{mt^2}{d(1-t^2)}.
}
\tag{24}
```

Thus photon number is $`\Theta(m)`$ at every fixed positive $`t<1`$. At $`d=6,t=1/2`$, both means lie between $`m/24`$ and $`m/18`$. The Takagi squeezing parameters are at most $`\operatorname{arctanh}t`$, independent of system size.

### 7.4 Linear support treewidth

For completeness, the expansion-to-treewidth implication can be proved directly. Suppose a connected $`d`$-regular graph has second-largest adjacency eigenvalue at most $`\lambda<d`$. For any vertex set $`S`$ of size at most $`m/2`$, the Rayleigh quotient of its centered indicator gives

```math
e(S,V\setminus S)
\ge(d-\lambda)\frac{|S|(m-|S|)}m
\ge\frac{d-\lambda}{2}|S|.
\tag{25}
```

A tree decomposition of width $`w`$ has a vertex separator $`X`$ of size at most $`w+1`$ such that every remaining connected component has at most $`m/2`$ vertices. To see this, assign each graph vertex to one bag containing it and take a weighted centroid bag in the decomposition tree; the connectedness condition on bags containing any vertex ensures that this bag separates the remaining components.

If $`|X|\ge m/4`$, the bound below is immediate. Otherwise choose a union $`S`$ of remaining components with $`m/4\le|S|\le m/2`$: take one component in this size range, or greedily combine smaller components. Every edge leaving $`S`$ ends in $`X`$. Hence

```math
d|X|\ge e(S,V\setminus S)
\ge\frac{(d-\lambda)m}{8}.
```

Applying $`\lambda=2\sqrt{d-1}`$ proves

```math
\operatorname{tw}(G)
\ge\frac{d-2\sqrt{d-1}}{8d}\,m-1.
\tag{26}
```

For the degree-six LPS family, the coefficient is $`(6-2\sqrt5)/48=0.0318305009\ldots`$. The support graphs therefore have linear treewidth. This statement concerns the graph of nonzero entries of $`B`$. It does not prove a lower bound on Gaussian entanglement entropy, on an optimized tensor-network bond dimension, or on unrestricted classical simulation cost.

## 8. Algorithmic scope and comparison

The constructive approximation in (21) uses the explicit positive matrix $`C=tA_G/d`$ and the nonnegative reference loop sampler. In ideal arithmetic, dense elimination preprocessing costs $`O(m^3)`$ operations and $`O(m^2)`$ storage; direct dense categorical draws cost expected $`O(m+m\overline N_C)`$ per independent count vector. Exploiting the fixed support degree in those draws reduces the latter to $`O(m+d\overline N_C)`$ after the same preprocessing. Equation (24) makes these output costs polynomial in the graph size. The dense preprocessing bound remains substantial and is not a practical large-instance benchmark. A finite-bit implementation must additionally budget its reference error as in (5).

The uniform-phase interval $`t<t_*`$ is much larger than $`1/[4e(d+1)]`$, which is $`1/(28e)=0.01313855\ldots`$ at degree six. That smaller number is a **conservative support-polymer sufficient condition derived in this project**, not a sharp published GBS boundary or a statement that polymer methods fail above it. Its assumptions, count decorations, and relationship to existing polymer algorithms are given in [the related-work analysis](../references/RELATED_WORK.md#a-conservative-support-polymer-corollary). The gain proved here comes from summing backtracking walks exactly and using the long-cycle geometry and spectrum of this explicit graph family.

The generic signed-loop construction, Bass identity, and Ramanujan graph construction are established results. The photonic statements here are their explicit consequences for full-count laws, with the normalization factors, phase convention, and finite-size TV error retained. Priority and scope comparisons for GBS algorithms belong to [the related-work analysis](../references/RELATED_WORK.md).

Applying the same classical postprocessing channel to both target and reference contracts total variation. Thus the bounds remain valid after thresholding photon counts or independently thinning detected photons. Conditioning on a rare event is different: fixed-photon-number, collision-free, or application-specific postselection requires a separate success-probability and conditional-error analysis.

## References

1. **Gregory F. Lawler.** “Topics in loop measures and the loop-erased walk.” *Probability Surveys* **15**, 28–101 (2018). DOI: [10.1214/17-PS293](https://doi.org/10.1214/17-PS293). [Author's notes](https://www.math.uchicago.edu/~lawler/fall16.pdf); [arXiv:1709.07531](https://arxiv.org/abs/1709.07531). The complex loop and Poisson constructions provide the established measure-theoretic framework.

2. **Hyman Bass.** “The Ihara–Selberg zeta function of a tree lattice.” *International Journal of Mathematics* **3**(6), 717–797 (1992). DOI: [10.1142/S0129167X92000357](https://doi.org/10.1142/S0129167X92000357). The determinant identity (12) uses the exponent $`|E|-|V|`$ on the edge-matrix side.

3. **Anton Deitmar.** “Ihara Zeta Functions of Infinite Weighted Graphs.” *SIAM Journal on Discrete Mathematics* **29**(4), 2100–2116 (2015). DOI: [10.1137/140957925](https://doi.org/10.1137/140957925). [Author-hosted published PDF](https://www.math.uni-tuebingen.de/user/deitmar/arbeiten/SIAM.pdf); [arXiv:1402.4945](https://arxiv.org/abs/1402.4945). Section 4 and Theorem 4.1 give the local-system version of the weighted determinant formula.

4. **Bharatram Rangarajan.** “A Combinatorial Proof of Ihara–Bass's Formula for the Zeta Function of Regular Graphs.” In *37th IARCS Annual Conference on Foundations of Software Technology and Theoretical Computer Science (FSTTCS 2017)*, Leibniz International Proceedings in Informatics **93**, 46:1–46:13 (published 2018). DOI: [10.4230/LIPIcs.FSTTCS.2017.46](https://doi.org/10.4230/LIPIcs.FSTTCS.2017.46). [Open-access paper](https://drops.dagstuhl.de/storage/00lipics/lipics-vol093-fsttcs2017/LIPIcs.FSTTCS.2017.46/LIPIcs.FSTTCS.2017.46.pdf). Lemma 4 gives the nonbacktracking trace formula in terms of adjacency eigenvalues.

5. **A. Lubotzky, R. Phillips, and P. Sarnak.** “Ramanujan graphs.” *Combinatorica* **8**(3), 261–277 (1988). DOI: [10.1007/BF02126799](https://doi.org/10.1007/BF02126799). [Public primary-paper PDF](https://courses.grainger.illinois.edu/cs598sgt/sp2015/LPS.pdf). Section 2 gives the generators, Theorem 3.4 gives order and girth, and Theorem 4.1 proves the Ramanujan property.
