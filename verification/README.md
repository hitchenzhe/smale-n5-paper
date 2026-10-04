# Independent verification code

This directory contains reviewer-oriented checkers rebuilt from the manuscript. They are not the author's original `ak_*.py` programs; those programs and their stored witness files were not present in the working folders examined.

## Run the current exact checks

First confirm that the current LaTeX tree matches the source revision recorded in `claim_status.json`:

```text
python verification/verify_source_revision.py --source-root C:/Users/hitch/Desktop/Smale
```

Requires Python 3.10 or later. The algebra checker uses only the standard library and exact `Fraction` arithmetic:

```text
python verification/verify_exact_algebra.py
```

It reconstructs the normalized five-body polynomial equations and checks the listed D1 and D6 identities by exact expansion. Each check raises an error and exits unsuccessfully if a polynomial difference is nonzero. D1 exponent-cone decompositions are exact rational certificates. The D6 primitive-wing norm identities are determinants over a quotient polynomial ring.

The LP certificate checker verifies rational feasible points and Farkas infeasibility witnesses:

```text
python verification/verify_lp_certificates.py path/to/certificate.json
```

Encode every rational as an integer or a string such as `"-7/12"`; decimal JSON numbers are rejected. Rows use `A x <= b`. For strict inequalities, the eventual D5 model must add a shared variable `s`, encode `0 < s <= 1`, and provide either a feasible witness or an infeasibility certificate for every branch.

The interval certificate checker uses exact rational rectangular complex intervals, including exact polynomial evaluation, Krawczyk images, a supplied exact inverse for the preconditioner, a contraction bound, and binary box-split trees:

```text
python verification/verify_interval_certificates.py path/to/boxes.json
```

It does not use floating point or claim to replay D2 yet. The D2 manuscript's actual boxes, excluded boundary strata, polynomial systems, and stored enclosures still need to be reconstructed and entered as certificates.

## Current verification coverage

| Manuscript location | Status |
|---|---|
| Common normalized equations | Exact mass-centre, torque, and virial identities pass |
| D1 Lemma 5.1 | The 20 weighted rows, boundary coefficients, absence of centre variables at exponent zero, and exact cone certificates pass |
| D1 Lemma 5.2 | All three polynomial syzygies pass |
| D1 Lemma 5.11 | Triangle virial identity and nine differential coefficients pass by exact expansion |
| D3 Lemma 7.12 | The explicit twisted family satisfies all 10 coordinate rows, 10 distance rows, and the gauge after exact ideal reduction |
| D5 Lemmas 8.20–8.21 | Reusable exact certificate checker exists; the 111 regime pairs, 27 patterns, 260 four-product cases, and their certificates have not yet been rebuilt |
| D2 interval calculations | Reusable exact rational complex-box/Krawczyk checker exists; the manuscript's actual cover and certificate data are not present |
| D6 Lemma 9.2 | Polynomial degree bounds, boundary rows, weighted centre identities, and wing separation pass; the chart substitution for item (ii) remains |
| D6 Lemma 9.4 | The torque identity passes with an explicit polynomial remainder `Q` |
| D6 Lemma 9.6 | The z-wing formula satisfies its nine boundary rows modulo the exact shape polynomial; completeness and exceptional-scale counts remain |
| D6 primitive norm identities | Four norm products appearing in the Lemma 9.9 table pass |
| Remaining D1–D8 computer-assisted claims | Not yet independently replayed |

In particular, a `PASS` in this directory is local to the named algebraic identity. It does not certify the germ coverage arguments or the main finiteness theorem. D2's interval covers need their exact boxes and an outward-rounded certificate checker; D3 needs its jet data and modular unit minors; D4 needs root-isolation witnesses; D5 needs a complete rational branch model and one checked LP witness per branch; D6 needs the full Jacobian/jet replay; D8 needs its face, rank, and leading-order certificates.

## Reviewer package still to build

1. Tie every checker input to a manuscript label and a source revision hash.
2. Rebuild the finite symbolic checks chapter by chapter and serialize compact certificates.
3. Reconstruct the full D5 order constraints and exhaustive branch enumeration, then replay all LP witnesses and Farkas certificates.
4. Reconstruct D2's actual box cover and boundary strata; extend the interval reader with Pellet/Rouché witnesses where the manuscript needs root counts.
5. Recompute D3/D4/D6/D8 algebraic tables and jets, checking denominators and minors at every specialization claimed in the proof.
6. Add a single command that refuses to report overall success while any main-proof claim or required certificate is missing.
