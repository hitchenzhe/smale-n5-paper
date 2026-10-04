# Finiteness of Planar Central Configurations of Five Bodies

This repository contains the LaTeX source and compiled PDF for Zhe Chen's draft paper, “Finiteness of planar central configurations of five bodies for all positive masses.”

## Files

- `main.tex` — paper entry point and shared definitions
- `sections/d1.tex`–`sections/d8.tex` — the eight diagram sections included by `main.tex`
- `main.pdf` — compiled 57-page draft

## Build

Compile `main.tex` with pdfLaTeX three times so the table of contents and cross-references settle. The source uses standard AMS LaTeX packages and the section files in `sections/`.

The computer-assisted proof programs described in Appendix A are not included in this repository.

## Partial verification (in progress)

The `verification/` directory contains reviewer-oriented exact checkers that are still being updated. Their current coverage and limitations are documented in `verification/README.md`; they do not constitute a verification of the full theorem.

From the repository root, the current commands are:

```bash
python verification/verify_source_revision.py --source-root C:/Users/hitch/Desktop/Smale
python verification/verify_exact_algebra.py
```
