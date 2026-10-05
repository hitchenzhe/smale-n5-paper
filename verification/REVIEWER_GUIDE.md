# Reviewer guide: reconstructed computational proof obligations

This package reconstructs the computations from the manuscript. It is not the original author's code. Run verify_all.py and inspect certificates/verification_run.json for the exact revision, executed checks, failures and unfinished obligations. The seven previously unfinished groups are covered below. The analytic reductions from these calculations to the theorem remain mathematical arguments for reviewer inspection.

## 1. D2 mass classification and rays

The eight shape/cancellation resultants are reconstructed by d2_mass_resultants.py, using exact Sylvester determinants and bounded integer Kronecker encoding. Checked degree/coefficient bounds make decoding injective. L factors, P even/odd parts, K polynomials, constants and body exchanges are verified.

In d2_mass_rays.py, six ordered positive-mass sectors are covered by rational Bernstein leaves outside isolated holes. Uniform rational Rouche disks exhaust all shape roots in each hole; cancellation excludes all but the required conjugate pair. A four-real-coordinate rational complex Krawczyk enclosure gives existence and uniqueness of the mass/shape root. Conjugation and mass permutations supply the full ray list. d2_ray_response.py independently checks all 72 responses through rebuilt Jacobians, interval inverses/Neumann bounds, and exhaustive cubic-scale branches.

## 2. D2 whole regular mass sheets

These exact global arguments replace the original large-ball covers and numerical boundary margins.

- d2_regular_X.py uses reciprocal-polynomial Bezout identities to exclude X=0 for all unit phase ratios on S1/S3/S4/S6/S8; body exchange covers S2/S7.
- d2_regular_euler.py splits the linear source into (I+2L)U and (I-L)V. An additional singularity requires eigenvalue -1/2; zero response requires eigenvalue 0. Exact resultants restrict positive exceptions to permutations of (5,1,1) and (2+2sqrt(2),1,1). Bernstein sector bounds and a positive-definite rational interval Hessian isolate the latter. Valid shapes have unit modulus and heavy-amplitude dominance, incompatible with signed cancellation for all unit phases.
- d2_regular_lagrange.py reduces zero response on S6/S8 to a weighted Laplacian cofactor. Its real-linear mass equations have a nonzero interval-certified determinant. Two rational Rouche disks exhaust each quadratic shape polynomial; all possible mass zeros lie in interval boxes with strict amplitude dominance, hence cannot cancel.
- d2_S5_blocks.py treats the entire S5 sheet through exact Q(omega)(a,b) block inverses, positive denominator norms and sum-of-squares factors. The one-different fold has an exact rank-eight minor, centred kernel, nonzero augmentation and simple torque zero for all positive masses.

These are calculations over entire positive mass domains, not floating-grid searches.

## 3. D2 exceptional branches and all K5 labels

### K1/K2

In d2_K1_branches.py, the first linear subresultant is a^4(C1*x+C0), and the last is a^10*K1. Joint Bernstein covers of C1,C0,K1 prove that this linear subresultant is never identically zero at a positive K1 point. Complementary blow-up charts cover critical corners, including the genuine positive-mass chart face. The nonzero leading coefficient of the second polynomial excludes infinity. A finite common root is therefore unique and real.

An exact algebraic parameter replaces the unspecified original rational parameter. Write W=(-r^2,-1,c^2). The positive square-coordinate shape root is unique and simple. The interval is 0<r<sqrt(t0), where t0 is the unique root in (2/5,9/20) of

```
1 - 2*t^2 - 8*t^3 + t^4 - 8*t^5.
```

Then r<c<1. The two positive mass kernels are

```
((1-c^3)*r, c^3-r^3, (1-r^3)*c),
((1+c^3)*r, c^3-r^3, (1+r^3)*c).
```

The code checks centering, signed cancellation and nonzero mass-solve determinants. Projective immersion is checked independently by the two tangential derivatives of the mass ratios: their shape resultants have a gcd with no zero for 0<r<1. This avoids assuming that every algebraic shape admits a finite normalized scale. Uniqueness of the common root and a single vanishing signed cancellation factor give the smooth positive branches. K2 follows by body exchange. The global response proof applies at every finite compatible label.

### K4

In d2_K4_arcs.py, the rational parameter is

```
A=4*k*(1+k^4),
B=(1-k^2)*(1+6*k^2+k^4),
C=(1+k^2)*(6*k^2-1-k^4).
```

On 5^(-1/4)<k<1, exact identities give A>C>B>0. A positive sum-of-squares derivative and endpoint limits exhaust every valid shape ratio. Two positive mass kernels and body exchange give all four arcs. Exact gcds exclude interior stationary mass projections.

### K5

In d2_K5_labels.py, exact square-root relations and the isolated positive degree-eight mass root enumerate all 27 phase classes, triangle roots, 64 cross signs and pair cube roots: 9216 labels, 144 cancellation labels, 9072 noncancellation labels. The cancellation labels reduce to the verified ordinary S5 response or simple fold torque zero. Torque and the independently checked finite boundary fibre exclude the rest.

The original exceptional h^3/h^4 jet assignments are unnecessary for this route and are not reproduced. Equal-triangle-mass ordinary/fold jets that are used elsewhere are checked separately in d2_seed_jets.py, with all selected rows, denominator units and order-three/order-four minors replayed in the specified algebraic coefficient field.

## 4. D5 coefficient audit and order model

The 25 variables are explicitly named in d5_valuation.py. Rows encode necessary frame, force, C/R/A cancellation, T/B/V, shared-source and dominant central/deep relations. Each row has a descriptive label. At possibly cancelling coefficients, conservative lower bounds enlarge the feasible set. Infeasibility of this enlarged set suffices for exclusion. A feasible LP point is an order assignment, not an analytic germ.

Tropical three-term minimum alternatives partition ties: each tie belongs to one branch, with strict inequalities in later branches. Five slice-sensitive monopole corrections are treated using valid lower bounds. Pre-dipole aggregate orders are used in the raw far-forcing bounds.

The independent d5_source_audit.py checks centered real-offset formulas and nonvanishing; the near Phi/Psi, B and virial identities; emptiness of coincident-both by exact eliminants/positive-root counts; rank-two R planes and their common kernel; and all inner signs and the central deep leading-coefficient contradiction. Far-source coefficients and source uniformity are checked by the corresponding Table 5 modules. Compare these identities and necessary constraints with the manuscript's regime reductions.

## 5. D5 all LP census decisions

The ordered index in d5_census.py has 111 regime orbits and 27 patterns: 2997 cases. It reconstructs every row/disjunction and replays d5_census_inputs.json. An infeasible tree covers every alternative at each node. Every leaf is an exact rational Motzkin contradiction. Each feasible record has a complete branch path and rational point satisfying all chosen rows, strictly when required.

There are 32286 Motzkin leaves. Exactly the 26 source survivor entries remain, with no extra or missing entry. Two missing OCR table rows were restored from the D5 screenshot: inner/inner RRR and inner/coincident variants RCR. Canonical pair orientation maps the latter to coincident/inner RRC.

rational_simplex.py can regenerate certificates with exact two-phase rational simplex, free-variable splitting and common strict slack. Equality elimination preserves nonnegative provenance in original rows. Replay does not trust the solver's conclusions.

## 6. D5 all four-product decisions

The elementary subset excludes 194 cases; d5_dominant.py supplies 60 further exclusions through 4476 Motzkin leaves. d5_four_complete.py verifies the duplicate-free exhaustive 260-case index and all six full necessary-model inner witnesses. Their four-product order is -1/4. All selected V constraints, including exact inner second virial order, are checked. The independent leading-coefficient contradiction excludes the inner RRR analytic germs. LP feasibility is not used as a germ-existence claim.

## 7. D6 all leading jets and orbit labels

d6_states.py exhausts 576 states and 29 symmetry orbits: 21 regular and eight corank representatives. d6_jets.py decodes safe jet coefficients in Q(i,lambda0), represented by theta=i+lambda0 with

```
theta^8+14*theta^4-48*theta^2+49=0.
```

Mass and cubic-scale powers are explicit. Every jet is substituted into all 22 equations through its leading order 4, 3 or 2. Every inverse denominator and leading coefficient is proved to be a nonzero product of finite-atlas units. Verified symmetries transport the result to every corank state; cubic relations cover every scale branch. Atlas/rank checks cover every finite label and the relevant exceptional faces.

## Supplementary D8 Remark 11.10

Both Table-5 supplementary programs are included in the full entry point. They compute conditional orders on hypothetical solution rays with finite chart coordinates and a finite nondegenerate lower wing; they do not assert the existence of a Diagram-8 germ.

`ak_diagram8_ray_diagrams.py` reconstructs all ten squared-distance ties from the chart. It checks the source C and D forcing supports, positive leading centres and gauges, every separation/component order, and every complementary-triangle four-product. Orders are affine rational functions of p. Exact endpoint values certify every strict comparison over the entire open interval -4<p<2, including possible cancellations handled by position lower bounds. The three squared-distance orders are 4-2p (upper wing), -6 (crossed pairs), and 0 (lower wing). Eight four-products have the minimum -8-2p; the two others have orders 4-2p and 12-6p. The minimum is strictly negative throughout the interval.

`ak_diagram8_pattern_boundaries.py` checks p=(2r-4)/(r+1), its strictly positive derivative, endpoints, and corner-face orders. On the full open interval, the z-strokes are 12/15/25 without z-circles and the w-strokes are 34/35/45 with w-circles at 3/4. At p=2 the z-circles at 1/2 join the upper z-strokes, giving exactly the Diagram-6 signature. Above p=2 the z-circles have no z-stroke, contradicting Rule 1a. At p=-4 the upper w-strokes join the order of W34; below it they dominate. Those patterns are not Diagram 8. All comparisons and all ten products are checked exactly, without parameter sampling.

Certificates: `d8_ray_diagrams.json` and `d8_pattern_boundaries.json`. The original four D8 computational lemmas remain separate from these supplementary claims.

## Scope of a successful run

A successful run establishes the declared computational identities and exclusions for the recorded revision, including documented stronger replacements. The diagram classification, analytic implications, and correspondence of manuscript hypotheses with reconstructed models remain mathematical arguments for review. This is not a theorem-prover formalization of the full main theorem.
