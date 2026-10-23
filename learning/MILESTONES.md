# Thirty milestones for Gaussian Loop Sampling

This map connects the optical model to the complete classical algorithms and their guarantees. Its stable milestones describe the project's scientific dependencies and understanding targets for any reader. Proofs, full calculations, and research attribution remain in the linked canonical pages.

The single educational anchor is **Weedbrook et al., _Gaussian Quantum Information_**. Start with [the selected anchor and notation bridge](ANCHOR.md), which identifies the relevant sections and translates the review's covariance and squeezing conventions into the pairing matrix used here. Section labels below follow the review's published PDF. The loop, probability, algorithmic, and graph arguments are supplied by this repository and its cited primary research.

## Why ten background milestones and twenty narrative milestones?

**Background/preliminaries, M01–M10**, develops three connected foundations: the optical state and measurement; the probability law and computational task; and the matrix, walk, and Poisson language. It proceeds from one mode to an integrated specification of the many-mode problem. Ten milestones suffice because the project uses a narrow Gaussian setting: pure states with zero displacement and photon counting.

**Narrative/innovation, M11–M30**, follows a longer dependency chain. A local loop weight must become an executable sampler, retain its phase information, support certified estimation or distributional approximation, survive finite precision, and connect to explicit graph families. The final milestones integrate implementation evidence, prior work, and operational scope. This phase therefore receives twenty milestones. The allocation reflects those distinct scientific steps.

## Navigation

| Route | Milestones | Destination |
|---|---|---|
| Optical objects to a computational problem | [M01](#m01)–[M10](#m10) | State, counts, reference stability, and the complete input/output contract |
| Loop mechanics to certified algorithms | [M11](#m11)–[M19](#m19) | Positive sampling, phase estimators, distributional bounds, and finite-bit implementation |
| Local interference to scalable graph families | [M20](#m20)–[M26](#m26) | Small phases, short loops, girth, Bass identities, Ramanujan spectra, and LPS examples |
| Scientific and operational synthesis | [M27](#m27)–[M30](#m30) | Detection, conditioning, evidence, comparison, and a complete claim with paid costs |

The dependency field in each milestone identifies the earlier ideas it uses. M10 integrates the background phase; M30 integrates the whole project. For a direct route to the results, use [the model's result map](../science/model.md#reading-route-and-result-map). The milestones supply conceptual connections to those same results.

## Phase I — Background/preliminaries

### M01

**One optical mode and two descriptions**  
**Phase:** Background/preliminaries.

**Central question.** How can the same optical mode have continuous quadratures and discrete photon outcomes?

**Understanding target.** Relate creation and annihilation operators, quadratures, the vacuum, and number states. Recognize that a finite number of modes still has an unbounded set of count vectors. Keep the review's convention $`[q,p]=2i`$ and vacuum covariance $`I`$ consistent with its phase-space formulas.

**Dependencies and story role.** This is the local physical starting point. It establishes the objects that later appear in the Gaussian state and in the sampled output.

**Conceptual connection.** The algorithms simulate a distribution of number measurements on bosonic modes. Their auxiliary walks are classical computational variables; identifying the physical measurement first gives those variables a precise interpretation later.

**Source pointers.** Anchor §§I.D and II.A, through [Quadratures, covariance, and the state restriction](ANCHOR.md#quadratures-covariance-and-the-state-restriction); [State and photon law](../science/model.md#state-and-photon-law); [What the classical variables represent](../science/interpretation.md#what-the-classical-variables-represent).

### M02

**Gaussian states, squeezing, and the chosen physical domain**  
**Phase:** Background/preliminaries.

**Central question.** Which features of a Gaussian state are fixed by first and second moments, and which Gaussian states does this project include?

**Understanding target.** Distinguish displacement, covariance, purity, vacuum fluctuations, and squeezing. Explain why a squeezed vacuum is Gaussian and can be nonclassical. Identify the project's restriction to pure, zero-displacement states, and distinguish it from the larger class of mixed or displaced Gaussian states.

**Dependencies and story role.** Uses M01. It builds the local state preparation that will be combined across modes in M03.

**Conceptual connection.** A compact phase-space description specifies the state, while the photon-count task still requires a probability law for a non-Gaussian measurement. Positivity of a Gaussian Wigner function and positivity of the later loop weights are different statements.

**Source pointers.** Anchor §§II.A.1 and II.B.1–II.B.3, with [Squeezing sign and passive mixing](ANCHOR.md#squeezing-sign-and-passive-mixing); [State and photon law](../science/model.md#state-and-photon-law); [Gaussian does not mean an easy photon distribution](../science/interpretation.md#gaussian-does-not-mean-an-easy-photon-distribution).

### M03

**From optical preparation to the pairing matrix**  
**Phase:** Background/preliminaries.

**Central question.** How do single-mode squeezers and a passive interferometer become one complex symmetric matrix?

**Understanding target.** Interpret the Takagi form $`B=U\,\mathrm{diag}(s_j)U^{\mathsf T}`$, with $`s_j=\tanh\xi_j`$, as the supplied pure Gaussian state. Explain the condition $`\|B\|_2<1`$, the role of mode phases, and the translation between covariance data and the pairing representation. Keep the review's squeeze-operator sign convention consistent with the bridge.

**Dependencies and story role.** Uses M01–M02. It moves from individual optical elements to the global many-mode input.

**Conceptual connection.** The classical algorithms receive an explicit description of $`B`$. No quantum state preparation or amplitude oracle is assumed as a computational resource. The matrix encodes both pair strengths and interference phases.

**Source pointers.** Anchor §§II.A.2, II.B.4, II.B.6, and II.C.2; [Translating covariance into the pair matrix](ANCHOR.md#translating-covariance-into-the-pair-matrix); the canonical [Covariance and pairing conventions](../science/model.md#covariance-and-pairing-conventions); [State and photon law](../science/model.md#state-and-photon-law).

### M04

**The complete photon-count law**  
**Phase:** Background/preliminaries.

**Central question.** What exactly is drawn when the algorithm returns a photon-count vector?

**Understanding target.** Interpret $`n\in\mathbb N^m`$, repeated rows and columns in $`B[n]`$, and the factorial factors in the hafnian probability law. Explain even total parity, the vacuum outcome, and collisions in a mode. Distinguish complete counts from threshold clicks and from a law conditioned on a fixed photon total.

**Dependencies and story role.** Uses M01–M03. It gives the many-mode input a precise experimental output, including events that a collision-free model would omit.

**Conceptual connection.** The hafnian describes pairings in an individual output probability. This does not force a sampling algorithm to evaluate that hafnian. Photon detection also explains why Gaussian-state simulation with Gaussian measurements does not already settle this task.

**Source pointers.** Anchor §II.E.4, contrasted with §§II.E.1–II.E.3 through [ANCHOR.md](ANCHOR.md); [State and photon law](../science/model.md#state-and-photon-law); [Gaussian does not mean an easy photon distribution](../science/interpretation.md#gaussian-does-not-mean-an-easy-photon-distribution).

### M05

**Generating functions connect counts to matrices**  
**Phase:** Background/preliminaries.

**Central question.** How can a determinant encode a distribution whose individual probabilities contain hafnians?

**Understanding target.** Read the multivariate probability generating function as a collection of all count probabilities. Understand why inserting mode fugacities rescales the pairing matrix, how normalization appears as a determinant ratio, and how derivatives give the mean total photon number. Distinguish an identity for a generating function from an executable sampling procedure.

**Dependencies and story role.** Uses M03–M04. It supplies the algebraic bridge through which the later loop construction reaches the physical count law.

**Conceptual connection.** Gaussian normalization makes a global probability object accessible through matrix algebra. Expanding its logarithm will organize that object into additive contributions; the probabilistic meaning of those contributions is established in the narrative phase.

**Source pointers.** [Generating function](../science/model.md#generating-function), including its derivation and photon-mean formula; [A complete one-mode example](ANCHOR.md#a-complete-one-mode-example); [Determinants, photon counts, and loop soups](../references/RELATED_WORK.md#determinants-photon-counts-and-loop-soups). The determinant photon-count identity is a project bridge beyond the anchor's Gaussian foundations.

### M06

**Three computational tasks and their error guarantees**  
**Phase:** Background/preliminaries.

**Central question.** What changes when the desired output is a sample, an expectation, or one specified probability?

**Understanding target.** Distinguish total variation accuracy for an entire count law, additive estimation of a bounded callable observable, and additive or relative evaluation of a particular probability. Interpret a failure probability and the cost of computing the observable. Explain why an additive estimate of a tiny event probability need not give useful relative accuracy.

**Dependencies and story role.** Uses M04–M05. It establishes the success conditions against which every subsequent theorem is read.

**Conceptual connection.** A signed Monte Carlo output can estimate expectations even when its unweighted count vectors follow another law. This distinction later separates the phase estimators from the positive approximation samplers.

**Source pointers.** [The question](../science/model.md#the-question); [1. Input, output, and cost parameters](../science/estimators.md#1-input-output-and-cost-parameters); [What each resource actually pays for](../science/interpretation.md#what-each-resource-actually-pays-for).

### M07

**Physical stability and reference stability**  
**Phase:** Background/preliminaries.

**Central question.** Why does the algorithm require more than a normalizable physical state?

**Understanding target.** Distinguish physical stability $`\|B\|_2<1`$ from the stronger reference assumption $`\rho=\||B|\|_2<1`$. Define $`C=|B|`$, $`\Lambda_C`$, and the reference photon mean $`\overline N_C`$. Explain why phase cancellation can make the physical mean $`\overline N_B`$ much smaller than the quantity that controls source-loop work.

**Dependencies and story role.** Uses M03, M05, and M06. It supplies the central domain assumption and the first complete resource distinction.

**Conceptual connection.** Replacing complex entries by their magnitudes creates a positive reference problem. The benefit of positivity is paid for by its possibly larger occupation and its own distance from instability. Failure of the reference condition is a limit of this method, not a hardness conclusion.

**Source pointers.** [The stable magnitude reference](../science/model.md#the-stable-magnitude-reference); [1. Input, output, and charged resources](../science/finite-precision.md#1-input-output-and-charged-resources); [What each resource actually pays for](../science/interpretation.md#what-each-resource-actually-pays-for).

### M08

**Matrix powers as alternating closed walks**  
**Phase:** Background/preliminaries.

**Central question.** What walk does a term in a matrix trace represent?

**Understanding target.** Translate products of matrix entries into weighted walks and traces into closed walks. Interpret the two copies of each mode in $`T_C`$ and $`T_B`$, their even walk lengths, and projection back to physical occupations. Distinguish ordered rooted walks, periodic repetitions, immediate backtracks, and simple graph cycles.

**Dependencies and story role.** Uses M03, M05, and M07. It supplies the local combinatorial objects needed to expand the determinant correctly.

**Conceptual connection.** The doubled graph records the alternation between amplitudes and their conjugates. Counting visits gives photon occupations, while the ordered-root convention preserves the multiplicities needed for the trace-log identity.

**Source pointers.** [Lifted walks](../science/model.md#lifted-walks); [1. Conventions and the signed loop representation](../science/topology.md#1-conventions-and-the-signed-loop-representation). These walk conventions are supplied locally rather than assumed from the Gaussian review.

### M09

**Poisson collections and additive occupations**  
**Phase:** Background/preliminaries.

**Central question.** How does an exponential of independent local contributions become a probability distribution?

**Understanding target.** Interpret a finite measure on loop types as the intensity of a Poisson collection. Connect its exponential generating functional to sums of independently contributed occupation vectors. Distinguish the mean number of loops, $`\Lambda_C`$, from the mean number of visits or photons, $`\overline N_C`$. Understand thinning and superposition as operations on the intensity.

**Dependencies and story role.** Uses M05, M07, and M08. It supplies the probabilistic language needed to turn the determinant expansion into a constructive model.

**Conceptual connection.** The infinitely many possible loop types can have finite total intensity and finite expected occupation. Independence belongs to the auxiliary loop collection; physical photon modes may remain strongly correlated.

**Source pointers.** [Lifted walks](../science/model.md#lifted-walks), particularly the intensity, occupation, and generating-functional identities; [Proof of the estimator](../science/estimators.md#proof-of-the-estimator) for thinning; [Determinants, photon counts, and loop soups](../references/RELATED_WORK.md#determinants-photon-counts-and-loop-soups) for research provenance.

### M10

**The integrated problem specification**  
**Phase:** Background/preliminaries.

**Central question.** What must be specified before a claim of efficient Gaussian count simulation has a definite meaning?

**Understanding target.** Assemble the physical input, matrix-access model, measurement, computational task, error metric, and reference-stability promise into one coherent specification. Identify the roles of mode number, squeezing, input precision, reference photon mean, observable cost, and requested confidence. Locate the gap that remains: the determinant law still needs a samplable positive representation or quantified phase cancellation.

**Dependencies and story role.** Integrates M01–M09. It closes the background phase by connecting local optical objects and probability tools to the global algorithmic question.

**Conceptual connection.** The relevant notion of efficiency depends jointly on the physical family and the requested output. This specification is the common basis for comparing all later constructions and prior algorithms.

**Source pointers.** [The question](../science/model.md#the-question), [The stable magnitude reference](../science/model.md#the-stable-magnitude-reference), and [What each resource actually pays for](../science/interpretation.md#what-each-resource-actually-pays-for); [ANCHOR.md](ANCHOR.md) for the physical-to-computational bridge.

## Phase II — Narrative/innovation

### M11

**The positive loop law is the photon law**  
**Phase:** Narrative/innovation.

**Central question.** Why does a Poisson soup of weighted closed walks reproduce the complete positive-reference photon distribution?

**Understanding target.** Connect the trace-log expansion with the rooting factor $`1/(2L)`$ on the doubled graph, the product of reference edge weights, and the Poisson generating functional. Track periodic-walk multiplicities and physical visits carefully enough to explain why the resulting occupation vector has law $`P_C`$, including collisions and all even totals. Recognize the result as an inherited determinant/loop identity applied to the optical model.

**Dependencies and story role.** Uses M05 and M08–M10. It is the narrative's first mechanism: a local positive walk weight produces the global count distribution.

**Conceptual connection.** Equality of generating functions identifies the law, while the next milestones turn that representation into a direct algorithm. The representation bypasses evaluation of individual hafnians without changing the requested output.

**Source pointers.** [Lifted walks](../science/model.md#lifted-walks); [Result and provenance](../science/direct-sampler.md#result-and-provenance); [A complete one-mode example](ANCHOR.md#a-complete-one-mode-example), including the two decompositions of a four-photon event; [Determinants, photon counts, and loop soups](../references/RELATED_WORK.md#determinants-photon-counts-and-loop-soups).

### M12

**First-return paths become executable transitions**  
**Phase:** Narrative/innovation.

**Central question.** How can an entire weighted first-return excursion be sampled using local choices?

**Understanding target.** Interpret the trailing vertex set, harmonic vector, return weight $`q_i`$, and Schur-complement condition $`0\le q_i<1`$. Explain how the transformed transition probabilities telescope to the desired path weight. Identify why zero harmonic entries are unreachable and why the final repeated root contributes no additional occupation.

**Dependencies and story role.** Uses M07–M08 and M11. It resolves the local path-generation step hidden inside the loop representation.

**Conceptual connection.** Linear-system data encode sums over possible future paths. Their harmonic transformation converts those sums into normalized local probabilities, connecting matrix elimination to a stochastic excursion.

**Source pointers.** [First-return data](../science/direct-sampler.md#first-return-data); [Why the construction has the right measure](../science/direct-sampler.md#why-the-construction-has-the-right-measure). The distinct rational-bridge construction used for finite-bit guarantees appears later in M19.

### M13

**The complete Poisson loop is the phase-bearing object**  
**Phase:** Narrative/innovation.

**Central question.** Why must the sampler remember which excursions belong to the same loop?

**Understanding target.** Explain the Poisson number of clusters, the logarithmic-series number of first-return excursions within each cluster, and concatenation into a genuine loop. Understand when a negative-binomial aggregate preserves occupations and why it loses information needed for phase reweighting. Relate this loss to the failure of a cosine of a phase sum to factor into excursion cosines.

**Dependencies and story role.** Uses M09, M11, and M12. It completes the exact sampling object and preserves the structure required by every later phase theorem.

**Conceptual connection.** Two constructions can have the same positive occupation law while supplying different auxiliary randomness. An estimator's correctness depends on the representation from which its weights are derived, including this hidden partition.

**Source pointers.** [Poisson loops, not ungrouped excursions](../science/direct-sampler.md#poisson-loops-not-ungrouped-excursions); [Why the construction has the right measure](../science/direct-sampler.md#why-the-construction-has-the-right-measure); [4. Why the Poisson loop partition must be preserved](../science/estimators.md#4-why-the-poisson-loop-partition-must-be-preserved).

### M14

**The positive sampler's complete arithmetic bill**  
**Phase:** Narrative/innovation.

**Central question.** Which work is reused, which work grows with a draw, and which randomness assumptions are being charged?

**Understanding target.** Interpret reusable $`O(m^3)`$ dense arithmetic preparation and expected draw cost $`O(m+m\overline N_C)`$. Explain the sparse per-visit improvement, the continuing dense preparation cost, and the storage requirement. Distinguish expected work from a deterministic time cap and identify unit-cost elementary variates as an explicit ideal-arithmetic convention.

**Dependencies and story role.** Uses M07 and M11–M13. It turns the constructive representation into a resource statement, completing the positive sampler before complex phases are added.

**Conceptual connection.** The expected number of visited vertices equals the reference photon mean, so optical occupation controls path work directly. Near instability, that output-sensitive cost can grow even for a modest number of modes.

**Source pointers.** [Preprocessing and cost](../science/direct-sampler.md#preprocessing-and-cost); [Tail bounds and runtime promises](../science/direct-sampler.md#tail-bounds-and-runtime-promises); [9. Arithmetic complexity and bit complexity are different guarantees](../science/finite-precision.md#9-arithmetic-complexity-and-bit-complexity-are-different-guarantees).

### M15

**Gauge-invariant interference and the cosine moment**  
**Phase:** Narrative/innovation.

**Central question.** How do complex edge amplitudes change the positive loop representation?

**Understanding target.** Interpret loop holonomy, reversed orientations, and the resulting cosine $`c(w)`$. Explain why optical gauge transformations leave occupations and loop phases unchanged. Connect the signed generating function to the normalized cosine-product weight $`e^\Delta\prod_wc(w)`$ on a reference soup, its exact second moment $`e^\Gamma`$, and the determinant gap $`\Delta`$. Understand how a phase-doubled matrix makes $`\Gamma`$ computable through three log determinants.

**Dependencies and story role.** Uses M03, M08, and M11–M14. It introduces the local phase mechanism and the first global estimation guarantee.

**Conceptual connection.** Large entry phases can be removable by a gauge, whereas closed-walk phases survive. The representation can be signed even though the physical count distribution is positive; cancellations must therefore be handled through an estimator or a distributional approximation theorem.

**Source pointers.** [Holonomy and phase resources](../science/model.md#holonomy-and-phase-resources); [2. The signed loop representation](../science/estimators.md#2-the-signed-loop-representation); [The cosine-product estimator and its exact second moment](../science/estimators.md#the-cosine-product-estimator-and-its-exact-second-moment).

### M16

**Negative loop mass controls a bounded signed estimator**  
**Phase:** Narrative/innovation.

**Central question.** Which part of the phase-sensitive loop measure produces the sign overhead?

**Understanding target.** Define $`r=\int(-c)_+\,d\mu_C`$ and interpret absolute-cosine thinning. Explain the retained-loop sign, the normalization $`e^{2r}`$, and the uniform second moment $`e^{4r}`$. Connect this moment to the number of independent samples needed for a bounded observable, while retaining the source-loop and readout costs.

**Dependencies and story role.** Uses M06, M09, M13, and M15. It sharpens the phase resource from a cosine-product bound to the negative part of the loop intensity.

**Conceptual connection.** This is absolute-weight importance sampling applied to a specific efficiently samplable loop representation. The inequality $`4r\le\Gamma`$ compares uniform moment guarantees; it does not order every observable's actual variance or prove optimality among classical methods.

**Source pointers.** [3. Main theorem: negative loop mass](../science/estimators.md#3-main-theorem-negative-loop-mass); [Proof of the estimator](../science/estimators.md#proof-of-the-estimator); [What the theorem does and does not optimize](../science/estimators.md#what-the-theorem-does-and-does-not-optimize).

### M17

**From sign control to full-distribution approximation**  
**Phase:** Narrative/innovation.

**Central question.** When does a positive loop process approximate the whole complex-target count law?

**Understanding target.** Interpret positive-part thinning and its guarantee $`\mathrm{TV}(P_B,Q_+)\le\min\{1,(e^{2r}-1)/2\}`$. Compare it with the reference-law bound $`\mathrm{TV}(P_B,P_C)\le1-e^{-\Delta}`$. Explain how the signed-convolution argument controls probability coefficients and why closeness of generating functions alone would be insufficient.

**Dependencies and story role.** Uses M06, M09, and M15–M16. It moves from weighted expectation estimation to an ordinary random output with a certified distributional bias.

**Conceptual connection.** Small negative mass has two operational consequences: a low sign overhead for estimation and a small full-law error after deleting the negative intensity. The accuracy needed for one task need not match the accuracy needed for the other.

**Source pointers.** [3. Main theorem: negative loop mass](../science/estimators.md#3-main-theorem-negative-loop-mass); [Proof of the total-variation bound](../science/estimators.md#proof-of-the-total-variation-bound); [Domination by the positive reference law](../science/estimators.md#domination-by-the-positive-reference-law).

### M18

**Certifying the resource and controlling normalization bias**  
**Phase:** Narrative/innovation.

**Central question.** How can the sign resource be bounded without silently biasing the estimator?

**Understanding target.** Compare the inexpensive certificate $`r\le\Gamma/4`$ with the Fourier expansion of $`|\cos\phi|`$ and its explicit harmonic tail. Distinguish absolute determinant error from relative accuracy after cancellation. Explain why substituting an upper bound into $`e^{2r}`$ changes the expectation, and how a ratio of signed means avoids that normalizer while retaining a confidence guarantee.

**Dependencies and story role.** Uses M15–M17. It makes the mathematical resource usable when exact negative mass is unavailable.

**Conceptual connection.** A constant-additive upper certificate changes the ratio estimator's sample budget by a constant factor. A direct prefactor instead needs normalization precision tied to final accuracy. The Fourier cost is proportional to inverse certificate tolerance; it cannot generally be charged as logarithmic precision alone.

**Source pointers.** [6. A Fourier certificate with an explicit finite truncation error](../science/estimators.md#6-a-fourier-certificate-with-an-explicit-finite-truncation-error); [Determinant and input errors](../science/estimators.md#determinant-and-input-errors); [Normalization and a ratio alternative](../science/estimators.md#normalization-and-a-ratio-alternative).

### M19

**Finite-bit algorithms preserve the stated guarantees**  
**Phase:** Narrative/innovation.

**Central question.** How do infinite loops, real probabilities, and phase arithmetic become a bounded finite computation?

**Understanding target.** Assemble the input encoding, reference-mean bound, exponential length and photon tails, rational closed-walk bridges, finite categorical choices, and explicit error budgets. Explain why tiny transition probabilities need absolute-error control rather than relative division. Distinguish the unbiased ideal estimator from a bias-certified implementation, including the bounded ratio and positive-part constructions.

**Dependencies and story role.** Uses M06–M07 and M14–M18. It closes the algorithmic chain from a formal random process to finite resources.

**Conceptual connection.** Exponential tails allow very small failure budgets with polynomial cutoff lengths. Input rounding, sampling error, determinant precision, and observable rounding are separate contributions. The general complex-entry oracle extension additionally needs a known stability margin; zeros do not require an exact phase-detection oracle.

**Source pointers.** [3. Exponential tails with explicit cutoffs](../science/finite-precision.md#3-exponential-tails-with-explicit-cutoffs); [8. Finite-bit estimator theorem](../science/finite-precision.md#8-finite-bit-estimator-theorem); [10. Bounded absolute-cosine thinning and the ratio estimator](../science/finite-precision.md#10-bounded-absolute-cosine-thinning-and-the-ratio-estimator); [11. General complex-entry oracles and stability at zero entries](../science/finite-precision.md#11-general-complex-entry-oracles-and-stability-at-zero-entries).

### M20

**Small phases and quartic moment sensitivity**  
**Phase:** Narrative/innovation.

**Central question.** How can two phase-resource bounds scale differently while the photon number remains extensive?

**Understanding target.** Relate a gauge-fixed edge-phase bound $`\theta`$ to accumulated loop phases. Understand the quadratic order of $`\Delta`$ and quartic order of $`\Gamma`$ at fixed reference gap. Use the connected expander-product example to interpret the matching orders $`\Theta(m\theta^2)`$ and $`\Theta(m\theta^4)`$, including $`\theta=m^{-1/4}`$.

**Dependencies and story role.** Uses M07, M08, and M15–M18. It turns abstract resource parameters into a physical scaling regime with extensive photon number.

**Conceptual connection.** Squaring the local departure of a cosine from one changes the small-phase order of the moment bound. This improves a variance certificate relative to one exponential in $`\Delta`$; the comparison alone says nothing exclusive about the possibility of full sampling.

**Source pointers.** [Quartic sensitivity to small edge phases](../science/estimators.md#quartic-sensitivity-to-small-edge-phases), including the $`H\mathbin{\square}K_2`$ example; [An extensive-photon example](../science/interpretation.md#an-extensive-photon-example).

### M21

**Positive short loops also support sampling**  
**Phase:** Narrative/innovation.

**Central question.** What full-distribution conclusion follows when all sufficiently short loop cosines are nonnegative?

**Understanding target.** Connect the condition $`2K\theta\le\pi/2`$ to positive thinning through length $`2K`$. Interpret the omitted loop-intensity tail and the resulting total variation bound. Explain why the small-phase example in M20 admits a positive sampling approximation as well as efficient estimation.

**Dependencies and story role.** Uses M08, M17, M19, and M20. It connects local phase control to global sampling accuracy and qualifies the interpretation of the quartic regime.

**Conceptual connection.** Loop-length decay can make unfavorable long-loop interference negligible even when a crude determinant-gap bound is large. At a fixed nonzero phase allowance, however, positivity restricts the maximum admissible cutoff; arbitrary requested accuracy requires its own valid certificate.

**Source pointers.** [3. Positive short-loop thinning](../science/topology.md#3-positive-short-loop-thinning); [Quartic sensitivity to small edge phases](../science/estimators.md#quartic-sensitivity-to-small-edge-phases); [10.4 Finite-bit positive-part sampling](../science/finite-precision.md#104-finite-bit-positive-part-sampling).

### M22

**Girth removes short-loop phase dependence**  
**Phase:** Narrative/innovation.

**Central question.** How can graph structure suppress interference without requiring small edge phases?

**Understanding target.** Explain why a closed walk shorter than the support girth reduces through immediate backtracks, whose conjugate edge phases cancel. Track this statement on the doubled support and the even threshold $`h=2\lceil g/2\rceil`$. Interpret the resulting determinant-gap tail for $`C=tA/d`$ with arbitrary symmetric edge phases.

**Dependencies and story role.** Uses M08, M15, M17, and M21. It replaces perturbative phase smallness with a topological reason for short-loop positivity.

**Conceptual connection.** The relevant phase is accumulated around a cycle. High girth delays the first cycle capable of carrying nontrivial holonomy, while spectral stability suppresses sufficiently long walks. These mechanisms can coexist with large photon number and large support connectivity.

**Source pointers.** [4. Girth eliminates short-loop phase dependence](../science/topology.md#4-girth-eliminates-short-loop-phase-dependence); [Optical gauge and magnetic gauge](../science/topology.md#optical-gauge-and-magnetic-gauge); [2. A determinant certificate for full-distribution sampling](../science/topology.md#2-a-determinant-certificate-for-full-distribution-sampling).

### M23

**The Bass identity sums the backtracking contribution**  
**Phase:** Narrative/innovation.

**Central question.** How can a vertex determinant be reorganized into nonbacktracking cycle information?

**Understanding target.** Interpret the weighted directed-edge operator and the magnetic Ihara–Bass determinant identity. Distinguish the unequal-magnitude formula from its regular, unit-phase specialization. Explain the coupling substitution and why the effective nonbacktracking parameter $`q`$ is smaller than $`t`$, while phase-independent factors cancel in the determinant ratio.

**Dependencies and story role.** Uses M08, M15, and M22. It strengthens a bound based on all closed walks by resumming their backtracking structure.

**Conceptual connection.** The sampler still uses complete loops, including backtracks. The nonbacktracking representation is an analytical certificate for their net phase effect, not a replacement sample space whose occupations can be substituted without proof.

**Source pointers.** [5.1 Unequal magnitudes](../science/topology.md#51-unequal-magnitudes); [5.2 Unit phases on a regular support](../science/topology.md#52-unit-phases-on-a-regular-support), including the determinant-ratio expansion; [What is checked](../verification/README.md#what-is-checked) for the independent finite-matrix comparisons.

### M24

**Ramanujan spectra turn cycle traces into useful bounds**  
**Phase:** Narrative/innovation.

**Central question.** Which spectral information makes the nonbacktracking tail small on a large regular graph?

**Understanding target.** Connect adjacency eigenvalues to the nonbacktracking roots and separate exceptional eigenvalues from the spectral bulk. Track the doubled cover's multiplicities, including its two components for a bipartite support. Interpret the three terms in the explicit bound $`R_{m,d,h}(t)`$ and why replacing every eigenvalue by the bulk radius would be invalid.

**Dependencies and story role.** Uses M22–M23. It converts the determinant identity into a quantitative finite-size certificate.

**Conceptual connection.** Girth controls which short traces vanish, while a Ramanujan bound controls the surviving long traces. Their combination is stronger than treating every vertex as an independent worst-case source of long walks.

**Source pointers.** [6. Ramanujan spectra strengthen the cycle bound](../science/topology.md#6-ramanujan-spectra-strengthen-the-cycle-bound); [5.2 Unit phases on a regular support](../science/topology.md#52-unit-phases-on-a-regular-support) for the quadratic spectral relation.

### M25

**An explicit LPS family and its accuracy regime**  
**Phase:** Narrative/innovation.

**Central question.** Can the graph assumptions be met by a reproducibly specified sequence with a stated coupling range?

**Understanding target.** Identify the degree-six LPS construction, the allowed primes $`\ell`$, and the size $`m=\ell(\ell^2-1)`$. Connect the published girth and spectral bounds to the sufficient threshold $`t_*\simeq0.736604`$ and the size-dependent total variation guarantee. Distinguish its explicit finite-size table, fixed-coupling asymptotic rate, and boundary behavior.

**Dependencies and story role.** Uses M22–M24. It turns an abstract graph criterion into a concrete, constructible family with arbitrary edge phases.

**Conceptual connection.** The graph construction, coupling conversion, and approximation theorem jointly define the result. A guarantee that improves along this family as its size grows does not imply arbitrary accuracy for every fixed matrix at the same coupling.

**Source pointers.** [7.1 Graph construction and established graph properties](../science/topology.md#71-graph-construction-and-established-graph-properties); [7.2 Fixed-parameter accuracy](../science/topology.md#72-fixed-parameter-accuracy), including the exponent and finite-size table; [8. Algorithmic scope and comparison](../science/topology.md#8-algorithmic-scope-and-comparison).

### M26

**What the graph family means physically and structurally**  
**Phase:** Narrative/innovation.

**Central question.** Is the explicit sampling regime merely a vacuum limit or a small-treewidth construction?

**Understanding target.** Relate the fixed coupling to bounded squeezing and the proved bounds giving both physical and reference photon means of order $`m`$. Interpret the expansion-to-separator argument behind the linear support-treewidth lower bound. Distinguish support treewidth from a Gaussian entanglement measure or a lower bound against every classical algorithm.

**Dependencies and story role.** Uses M02–M03, M07, M14, and M25. It joins the graph accuracy theorem to the physical occupation scale and structural complexity of the input.

**Conceptual connection.** Extensive photons and large support treewidth make the family informative, while loop cancellations still admit a classical approximation. The sampling cost continues to include preparation and the reference occupation; these physical features do not erase that bill.

**Source pointers.** [7.3 Extensive photons at bounded squeezing](../science/topology.md#73-extensive-photons-at-bounded-squeezing); [7.4 Linear support treewidth](../science/topology.md#74-linear-support-treewidth); [What each resource actually pays for](../science/interpretation.md#what-each-resource-actually-pays-for).

### M27

**Detection channels and rare conditioning**  
**Phase:** Narrative/innovation.

**Central question.** Which experimental output changes preserve the guarantee, and which create a new cost?

**Understanding target.** Explain contraction of total variation under common classical postprocessing, including thresholding and independent output loss. Interpret a processed observable by composition with that channel. Distinguish such processing from conditioning on a rare event, where rejection and signed ratios acquire explicit dependence on the event probability.

**Dependencies and story role.** Uses M04, M06, M16–M19, and M26. It translates the main results into realistic readout choices while keeping their probability and resource contracts intact.

**Conceptual connection.** Coarse-graining cannot increase distinguishability, but postselection renormalizes a restricted distribution and can amplify error. An output-loss statement therefore has a precise location in the optical experiment; it does not cover arbitrary loss inserted throughout a circuit.

**Source pointers.** Anchor §§II.E.1–II.E.4 for measurement context, via [ANCHOR.md](ANCHOR.md); [8. Classical output processing, loss, and conditioning](../science/estimators.md#8-classical-output-processing-loss-and-conditioning); [Detection, conditioning, and extensions](../science/interpretation.md#detection-conditioning-and-extensions).

### M28

**Proof guarantees and numerical evidence have distinct roles**  
**Phase:** Narrative/innovation.

**Central question.** What does the reference implementation establish about a calculation, and which guarantees come from the proofs?

**Understanding target.** Relate exact small hafnian enumeration, explicit rooted-walk sums, matrix identities, and seeded sampling comparisons to the mechanisms they check. Distinguish the floating-point API, the ideal arithmetic theorem, and the constructive finite-bit theorem. Explain why discarded work-limit failures condition the output law, and why deleting a loop after generation does not cap generation time.

**Dependencies and story role.** Uses M11–M19 and M23–M27. It integrates the algorithmic and analytical claims with the reproducible numerical artifact.

**Conceptual connection.** Independent finite examples can expose normalization, phase, and multiplicity mistakes. Their agreement supports implementation correctness on those examples; analytical error budgets and proofs supply the general guarantees and asymptotic statements.

**Source pointers.** [What is checked](../verification/README.md#what-is-checked); [Reference API](../verification/README.md#reference-api); [Numerical and runtime limits](../verification/README.md#numerical-and-runtime-limits); [9. Arithmetic complexity and bit complexity are different guarantees](../science/finite-precision.md#9-arithmetic-complexity-and-bit-complexity-are-different-guarantees).

### M29

**The contribution in relation to existing algorithms**  
**Phase:** Narrative/innovation.

**Central question.** Which ingredients are established, and which comparison supports this project's scientific contribution?

**Understanding target.** Separate the inherited Gaussian determinant identity, Markov-loop construction, absolute-weight estimation, and Bass identity from their explicit optical application and resource accounting. Compare nonnegative full-count sampling, phase-sensitive observable estimation, structural simulation, and perturbative approximations under matched input and output promises. Interpret the conservative support-polymer threshold as a sufficient comparison bound.

**Dependencies and story role.** Uses M06–M07 and M11–M28. It places the complete mechanism and examples within the relevant research landscape.

**Conceptual connection.** A new representation or sharper certificate need not be the first polynomial algorithm, a universal speedup, or a hardness boundary. Weighted and collision-inclusive deductions, ensemble-average results, readout restrictions, and fixed-instance guarantees each require their own explicit comparison.

**Source pointers.** [Determinants, photon counts, and loop soups](../references/RELATED_WORK.md#determinants-photon-counts-and-loop-soups); [Nonnegative graph-based GBS](../references/RELATED_WORK.md#nonnegative-graph-based-gbs); [Recent phase-sensitive and observable-estimation results](../references/RELATED_WORK.md#recent-phase-sensitive-and-observable-estimation-results); [What this comparison establishes](../references/RELATED_WORK.md#what-this-comparison-establishes); [Claim boundaries](../references/RELATED_WORK.md#claim-boundaries).

### M30

**The complete operational account**  
**Phase:** Narrative/innovation.

**Central question.** Given a state family and a requested readout, what precise classical claim does the project justify?

**Understanding target.** Connect the physical model, access assumptions, stable reference, chosen output task, phase or topology certificate, preparation cost, draw count, observable cost, and finite error budget into one end-to-end statement. Select the relevant positive sampler, signed or ratio estimator, positive-part approximation, or reference-law approximation by its proved conditions. State the accuracy attainable for a fixed input or along a specified family.

**Dependencies and story role.** Integrates M01–M29. It closes the narrative by joining the local loop mechanism to its physical meaning, executable algorithms, evidence, significance, and limits.

**Conceptual connection.** The resulting scientific account explains when classical closed walks reproduce or approximate Gaussian photon statistics and what the computation pays for. Interference enters through representation-dependent loop resources; the project establishes sufficient simulation regimes with explicit guarantees.

**Source pointers.** [Reading route and result map](../science/model.md#reading-route-and-result-map); [7. End-to-end interpretation](../science/estimators.md#7-end-to-end-interpretation); [8. Algorithmic scope and comparison](../science/topology.md#8-algorithmic-scope-and-comparison); [Physical meaning and computational boundaries](../science/interpretation.md); [Claim boundaries](../references/RELATED_WORK.md#claim-boundaries).

