# Smale N=5 verification package (v19)

The current source is `../Smale_5-CT-v19.tex`. `claim_status.json` binds its exact SHA-256, while `program_index.json` maps source labels to implementations and documented replacements. Hashes establish version identity, not mathematical correctness.

Run from the manuscript directory with Python 3.10 or later:

```powershell
python verification/verify_source_revision.py --source-root .
python verification/verify_all.py --source-root . --output-dir verification/run-records/v19-current
python verification/verify_v19_changes.py --source-root . --output verification/run-records/v19-current/v19_changes.json
```

SymPy 1.14.0 and mpmath 1.3.0 are bundled as hash-checked wheels in `vendor/`; verification needs no network. The suite contains the original 77 exact/modular checks, five reconstructed supplementary v19 checks, and the separate exact rational D4 census: 83 checks, covering 47 applicable claim groups and one not-applicable group. Source/program preflight rejects changed source, dangling active labels, and undocumented replacements before expensive computations.

Exit 0 means the declared computational checks pass for the bound source; exit 1 means failure or revision mismatch; exit 2 means declared computational obligations remain incomplete. These codes do not certify the analytic proof. The manuscript explicitly states its main theorem conditional on its listed computational statements.

## Implementation scope

The package is reconstructed verification code. Nine historical program names in the earlier appendix have documented replacements rather than same-named Python files. In particular:

- D2: the 1,080 rank-20 fold points are excluded by the exact first torque obstruction and finite boundary fibre. The old rank-20/h^3 calculation is not claimed as reproduced.
- D2: global algebraic and rational interval certificates establish regular/fold exclusion consequences. Not every historical box cover, intermediate jet order, or printed numerical boundary margin is replayed.
- D3: all 180 equal-mass labels are replayed through order six over both F13 and F37; 144 have a nonzero raw residual at each prime. Nonzero modular residuals prove characteristic-zero exclusion; zero residuals do not prove existence.
- D6: representative charts, exact leading jets, norm factors, and orbit coverage are checked.
- `verify_v19_changes.py` independently reconstructs the Petersen matrix, D9 order arithmetic, normalized virial/D7 inertia identities, D5 transverse Hessian and contact arithmetic, D8 scalar elimination and quadratic pairing, and selected printed D2/D3/D6 nonvanishing identities. It does not import the `ak_*` model builders or their certificates.

Far/near case coverage, chart coverage, analytic continuation, local analytic lemmas, and the use of external AK results require mathematical reading. Supplementary computations do not replace those arguments. There is no assertion here of complete Lean verification.

## Execution records

`v19-records/verification_run.json` is the preserved 83-check execution on the v19 manuscript before its verification appendix was shortened. The mathematical body and checker implementations did not change during that editorial edit. `v19-records/editorial_update.json` binds the current source and its function table to that earlier execution; it explicitly records that a full computation was not repeated on the shortened source. Running the commands above creates fresh records for the current source under the ignored `run-records/` directory, preserving the published record.

The five supplementary identities in `verify_v19_changes.py` were rerun on the current source; their result is `v19-records/v19_changes.json`.

The prebuilt PDF is at the repository root. `../build.ps1` rebuilds it using local MiKTeX with automatic package installation disabled. The local build cache is not included in this repository.

See `REVIEWER_GUIDE.md` for the trust boundary and the treatment of replacements.
