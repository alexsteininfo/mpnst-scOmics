# mpnst-scOmics
In this repository, we provide code for bioinformatic analysis of single-cell sequencing data of malignant peripheral nerve sheath tumours.

Reference papers for the pilot cohort, sequencing technologies, and bioinformatic tools used (or planned) here are summarized in [literature/README.md](literature/README.md).

## Related repositories

- [`CNA-PhyloAnalysis`](https://github.com/alexsteininfo/CNA-PhyloAnalysis) (`/srv/home/aste0033/GitHub/CNA-PhyloAnalysis`) — downstream repository that takes the copy-number profiles / MEDICC2 phylogenies produced here and analyses tree topology (CNA burden, clonal dynamics, mutation rate, site-frequency spectra).
- [`MutationLoadDynamics.jl`](https://github.com/alexsteininfo/MutationLoadDynamics.jl) — Julia package simulating somatic evolution, used by `CNA-PhyloAnalysis` to generate synthetic ground-truth trees; not used directly by this repo.
