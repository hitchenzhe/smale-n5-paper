# Smale’s sixth problem for five bodies

This repository contains version 10 of the manuscript by Zhe Chen and Zhicheng Tong, its compiled PDF, and the computational verification package.

- [PDF](Smale_5-CT-v10.pdf)
- [LaTeX source](Smale_5-CT-v10.tex)
- [Verification code and reproduction guide](verification/README.md)

Version 10 adds an analytic argument handling Diagram 9 in response to the gap reported in the proof of Albouy–Kaloshin Proposition 8. The verification programs cover the computational claims in the Appendix for Diagrams 2, 3, and 6. They do not verify the new Diagram 9 argument, the analytic proof, or the cited external theorems. The computational Appendix and its inputs are unchanged from version 9.

This is a preprint. The PDF in this repository was built from the included TeX source.

## Build

With `pdflatex` available, run from the repository root:

```powershell
.\build.ps1
```

The script performs three no-shell-escape passes, writes auxiliary files under `verification/build-v10/`, and refreshes `Smale_5-CT-v10.pdf`.
