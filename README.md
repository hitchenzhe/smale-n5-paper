# Smale N=5 paper (v19)

The current manuscript is [Smale_5-CT-v19.tex](Smale_5-CT-v19.tex), with a [compiled PDF](Smale_5-CT-v19.pdf). This repository contains the verification functions, exact certificates, and execution records used with this version.

The appendix lists the functions used at each computer-assisted step. [verification/README.md](verification/README.md) gives the reproduction commands and scope of the checks; [verification/REVIEWER_GUIDE.md](verification/REVIEWER_GUIDE.md) explains the relationship between the historical program names and the reconstructed certificates.

The recorded 83-check execution was run on the v19 manuscript before the verification appendix was shortened. The mathematical body and checker functions were unchanged by that editorial edit. The current source hash and function table were checked separately; [editorial_update.json](verification/v19-records/editorial_update.json) records this relationship. The full execution record has not been relabeled as a fresh run against the shortened source.

Run `./build.ps1` to rebuild the PDF with a local MiKTeX installation. The script disables automatic package installation. The prebuilt PDF is included here.

The computational checks establish their stated certificate results; analytic arguments and cited external theorems require mathematical review. The main theorem has the conditional status stated in the manuscript.
