# Smale’s sixth problem for five bodies

This repository contains the version 9 manuscript by Zhe Chen and Zhicheng Tong, its PDF, and the computational verification package cited in the paper.

- [PDF](Smale_5-CT-v9.pdf)
- [LaTeX source](Smale_5-CT-v9.tex)
- [Verification code and reproduction guide](verification/README.md)

The v9 manuscript changes the global argument for Diagram 5 and clarifies the specialization used in the Diagram 3 coefficient identities. The Appendix lists the remaining computational inputs for Diagrams 2, 3, and 6. The verification package covers those computations and records its scope; it does not formalize the analytic proof or the cited Albouy–Kaloshin results.

This is a preprint. The manuscript says that manual review by the second author is underway and that a complete Lean proof is planned for a forthcoming version. The PDF in this repository was built from the included TeX source.

## Build

With `pdflatex` available, run from the repository root:

```powershell
.\build.ps1
```

The script performs three no-shell-escape passes, writes auxiliary files under `verification/build-v9/`, and refreshes `Smale_5-CT-v9.pdf`.

