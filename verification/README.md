# Smale manuscript: independent computational verification

本程序包根据正文重建，采用论文第 56 页 Table 5 的 ak_*.py 名称。原作者代码和原始证书没有提供。D2、D5、D6 的剩余计算义务已补入；最终状态由当前代码的完整运行记录决定。

## Reviewer quick start / 审稿人运行方法

Requires Python 3.10 or later. SymPy 1.14.0 and mpmath 1.3.0 are bundled as unmodified wheels in vendor/. No installation or network is required; their SHA-256 hashes are checked before loading. From the LaTeX project directory containing main.tex, sections/ and verification/:

```
python verification/verify_all.py
```

The run reconstructs equations and LP constraints, checks certificates, and verifies all nine LaTeX source hashes. It takes several minutes. For another project or output directory:

```
python verification/verify_all.py --source-root PATH_TO_LATEX_PROJECT --output-dir PATH_TO_RECORDS
```

| Exit code | Meaning |
|---|---|
| 0 | All declared computational proof obligations are covered and pass |
| 1 | A check fails, sources differ/are missing, or a PASS claim lacks executed checks |
| 2 | Checks pass but computational proof obligations remain unfinished |

Exit 0 concerns computational obligations. The analytic proof has not been formalized in a theorem prover. Review the reductions in [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) and the code.

## Coverage

| Chapter | Computational coverage |
|---|---|
| Common / D1 | Normalized equations, three syzygies, chart/cone certificates, torque, shape/sign elimination, triangle rank/kernel/inertia, virial differential |
| D2 | Eight mass resultants; positivity exclusions; full L1/L8 rays and 72 responses; general-M ordinary/fold jets; K1/K2 branches, K4 arcs, K5 all 9216 labels; entire regular and one-different fold mass sheets |
| D3 | All 5184 labels; Jacobian and jets through G5; mass families/sign counts; degree-66 norm/Sturm count; complete modular G6 replay; exact twisted family |
| D4 | Torque and 64-sign census with second obstruction at every candidate |
| D5 | Far/near coefficients, uniformity, empty regimes/R planes; all 2997 census decisions, exactly 26 survivor entries; all 260 four-product decisions and separate inner leading-coefficient exclusion |
| D6 | Chart/atlas; rank at every finite label; 576 states/29 orbits; all 22 jet rows at eight corank representatives, leading orders 4/3/2 and nonzero units; transport to all labels |
| D7 | No computer-assisted claim identified |
| D8 | Four computational lemmas plus both supplementary Remark 11.10 programs: corner rows/faces, beta identity, Jacobians, signs, cofactors/pairings, second-order constants, lower atlas determinant, entire ray-order window, sharp boundary patterns and all ten four-products |

Some original D2 interval covers and exceptional higher-jet assignments are replaced by stronger exact proofs. Their unused original numerical margins and jet orders are not asserted. Both supplementary D8 ray/boundary programs in Remark 11.10 are implemented and included in the full run, although the main Diagram-8 exclusion does not need them. See claim_status.json under superseded_computations.

## Paper program names and individual runs

program_index.json maps all 54 Table 5 names to modules, source labels, functions, certificates or documented replacements. There are no dummy PASS implementations.

```
python verification/verify_exact_algebra.py --list
python verification/ak_exceptional.py
python verification/ak_p_sheet_second_order.py
python verification/ak_p5_completion.py
python verification/ak_diagram5_first_limit_relevance.py
python verification/ak_diagram5_outer_census.py
python verification/ak_diagram6_blowup.py
python verification/ak_diagram8_ray_diagrams.py
python verification/ak_diagram8_pattern_boundaries.py
```

Dependencies in smale_verify/stages.json are executed automatically. Repeated --program options select programs in verify_exact_algebra.py. D4 has a separate entry point:

```
python verification/ak_diagram4_sign_census.py --summary --output verification/certificates/d4_lemma_8_10.json
```

D5 normally replays stored certificates. Independent regeneration takes longer and resumes only after replaying completed orbits:

```
python -c "import sys; sys.path.insert(0, 'verification'); from smale_verify.d5_census import generate; generate()"
```

## Records and trust boundary

- certificates/verification_run.json: authoritative full-run status, execution bindings, source/code/certificate hashes, failures and unfinished obligations.
- certificates/exact_algebra_record.json: detailed exact/modular checks with paper program/function names.
- certificates/d4_lemma_8_10.json: sign cases, rational brackets, eliminants and obstruction gcds.
- certificates/programs/: per-program records; older subsets do not replace the full-run record.
- claim_status.json: declared obligations and required check identifiers. Each required check must actually run and pass.

Certificates are untrusted inputs. LP rows and branch choices are rebuilt from code; every multiplier and point is checked against those rows, including strict inequalities. D6 coefficients are decoded from safe polynomial data and substituted in all 22 equations. Domain coverage uses exhaustive sectors/split trees and algebraic elimination, not numerical sampling.

Arithmetic uses exact rational Laurent polynomials, quotient rings, finite fields, algebraic coefficients, Sturm chains and rational Bernstein/interval bounds. D2 resultants use injective integer encoding with checked degree/coefficient bounds. D3 modular exclusions check unit inverses and nonzero residuals; surviving labels have an exact family. Hashes identify a revision, not mathematical correctness.

Reusable generic readers:

```
python verification/verify_lp_certificates.py path/to/certificate.json
python verification/verify_interval_certificates.py path/to/boxes.json
```

Rational inputs are integers or fraction strings. Floating JSON numbers are rejected. A generic reader PASS does not establish a manuscript claim without the reconstructed model and complete coverage.
