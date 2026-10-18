# Reviewer guide for version 9

`claim_status.json` binds the included code to `Smale_5-CT-v9.tex` by SHA-256 and lists the computational claims with their required check identifiers. `program_index.json` maps the paper’s source labels to implemented modules and explains the scope of each program.

The main theorem’s computational inputs are the entries for Diagrams 2, 3, and 6 in the Appendix. Other registered checks are additional cross-checks. They should not be read as replacing the analytic arguments in the paper.

Run the commands in `README.md` to check the source binding and regenerate the computational run records. The package checks its stated exact computations and certificates. It does not verify the analytic proof or the cited Albouy–Kaloshin results, and it is not a formal proof in Lean.

