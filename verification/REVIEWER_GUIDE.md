# Reviewer guide for manuscript v19

The bound source is `../Smale_5-CT-v19.tex`. Reproduction commands and the scope of replacements are in `README.md`. Current records live in `v19-records/`.

`program_index.json` distinguishes implemented modules, documented replacements and supplementary calculations. It maps current TeX labels, retaining removed labels only in `legacy_source_labels`. Preflight requires each active label to exist and each absent appendix program to have a documented replacement.

The original 42 applicable groups comprise 24 primary theorem-consequence groups for D2/D3/D6 and 18 supplementary groups. Five supplementary reconstructions were added for v19. PASS for a replacement means the stated exclusion consequence is covered; it does not claim a literal replay of every intermediate calculation or numerical margin of the historical program.

The suite reconstructs polynomial systems and replays rational/algebraic certificates. D3 covers both primes, every label, selected minors, all jet coefficients and raw residuals. D4 has a separate 64-sign rational Sturm census. The new supplementary checker builds its formulas directly from the manuscript without the old `ak_*` builders.

The execution record does not certify chart/case completeness, analytic deductions, external AK theorems, or an unconditional main theorem. Those interfaces require mathematical review rather than a machine check.

The latest editorial source binding is recorded in v19-records/editorial_update.json. The preserved verification_run.json is the preceding full 83-check execution; it has not been relabelled as a replay of the editorial change.
