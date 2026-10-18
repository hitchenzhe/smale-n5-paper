# Verification package

The scripts and certificate data in this directory support the computational statements listed in Appendix A of `Smale_5-CT-v9.tex`. The package is bound to that exact TeX source by `claim_status.json`.

## Reproduce the PDF

From the repository root, run `./build.ps1` in PowerShell. It uses three `pdflatex` passes with shell escape disabled.

## Reproduce the computational checks

Requires Python 3.10 or later. Install the bundled, pinned dependencies and then run:

```powershell
python -m pip install --no-index --find-links verification/vendor -r verification/requirements.txt
python verification/verify_source_revision.py --source-root .
python verification/verify_all.py --source-root . --output-dir verification/run-records
```

The replay writes its run records under `verification/run-records/`; these generated records are not included in the repository. Exit code 0 means all declared computational checks pass; exit code 1 indicates a failed check or source binding; exit code 2 means implemented checks pass while proof obligations remain.

## Scope

The Appendix identifies computational inputs in Diagrams 2, 3, and 6. `program_index.json` maps source labels to implemented programs, while `claim_status.json` records the check-to-claim bindings. The code uses exact rational, algebraic, modular, interval, and linear-programming computations as described by the individual modules and certificates.

These checks do not prove the analytic deductions in the manuscript, validate every logical use of a computational statement, or formalize cited external theorems. The package is not a complete Lean proof of the theorem. Generated build output and execution records are not included.

