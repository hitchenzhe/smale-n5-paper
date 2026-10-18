# Smale N=5: current manuscript and verification

The current manuscript is version 9:

- [`Smale_5-CT-v9.tex`](Smale_5-CT-v9.tex)
- [`Smale_5-CT-v9.pdf`](Smale_5-CT-v9.pdf)

The source hash is `37599CF4B81DF2DEAC4FD8273D760404433033C7CED05DD465CF0A0D240E963B`. The PDF was compiled from this source with MiKTeX pdfLaTeX in three passes. It has 68 A4 pages; the final pass has no unresolved citation or reference warnings.

## Verification material

The scripts and certificates under [`verification/`](verification/README.md) are the version 8 verification package. Its source-revision manifest is bound to [`Smale 5-CT-v8.tex`](Smale%205-CT-v8.tex), whose SHA-256 is `4484BFD28C7B9B9FC29A5580FE2CCDF372E45C37B2C7D89A774842F38092653B`. The recorded run reports 78 computational checks passed for version 8. It is not a run bound to version 9. The v9 manuscript changes its proof exposition and the Diagram 5 argument; its complete mathematics has not been formally verified by these scripts.

The checks are reconstructed computational verifiers. They do not formalize the analytic proof or the cited Albouy–Kaloshin results. See [`verification/README.md`](verification/README.md) and [`verification/V8_AUDIT.md`](verification/V8_AUDIT.md) for scope and reproduction details.

The earlier split-source draft (`main.tex` and `sections/`) is superseded by the versioned monolithic sources above.