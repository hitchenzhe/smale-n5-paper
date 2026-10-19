# Reviewer guide for manuscript version 10

`claim_status.json` binds the included code to `Smale_5-CT-v10.tex` by SHA-256 and lists the computational claims with their required check identifiers. `program_index.json` maps the paper’s source labels to implemented modules and explains the scope of each program.

The computational inputs in the manuscript are listed for Diagrams 2, 3, and 6. Other registered checks are supplementary cross-checks. Version 10 also adds an analytic argument for Diagram 9; that argument is outside the scope of this code package, as are the remaining analytic arguments and the cited Albouy–Kaloshin results.

Run the commands in `verification/README.md` to check the source binding and regenerate computational run records. A passing run covers only the stated computations and certificates. It is not a formal proof in Lean.
