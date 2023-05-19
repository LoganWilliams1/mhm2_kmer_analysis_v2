# MetaHipMer2 - MHM2 #

[MetaHipMer (*MHM*)](https://sites.google.com/lbl.gov/exabiome/downloads?authuser=0) is a *de novo* metagenome short-read assembler. This is version 2 (MHM2), which is written entirely in
[UPC++](https://upcxx.lbl.gov), CUDA and HIP, and runs efficiently on both single servers and on multinode supercomputers, where it can scale up to
coassemble terabase-sized metagenomes. More information about MetaHipMer can be found under the [ExaBiome Project](https://sites.google.com/lbl.gov/exabiome) 
of the [Exascale Computing Project](https://www.exascaleproject.org) and in several publications:

- [E. Georganas et al., "Extreme Scale De Novo Metagenome Assembly," SC18: International Conference for High Performance Computing, Networking, Storage and Analysis, Dallas, TX, USA, 2018, pp. 122-13.](https://ieeexplore.ieee.org/document/8665813)
- [Hofmeyr, S., Egan, R., Georganas, E. et al. Terabase-scale metagenome coassembly with MetaHipMer. Sci Rep 10, 10689 (2020).](https://www.nature.com/articles/s41598-020-67416-5#citeas)
- [Awan, M.G., Deslippe, J., Buluc, A. et al. ADEPT: a domain independent sequence alignment strategy for gpu architectures. BMC Bioinformatics 21, 406 (2020).](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03720-1#citeas)
- [Muaaz Awan, Steven Hofmeyr, Rob Egan et al. "Accelerating large scale de novo metagenome assembly using GPUs.", SC 2021](https://dl.acm.org/doi/pdf/10.1145/3458817.3476212)

The quality of MetaHipMer assemblies is comparable with other leading metagenome assemblers, as documented in the results of the CAMI2 competition, where MetaHipMer was rated first for quality in two out of three datasets, and second in the third dataset:

- [F. Meyer et al., "Critical Assessment of Metagenome Interpretation: the second round of challenges", Nature Methods volume 19, pages429–440 (2022)](https://www.nature.com/articles/s41592-022-01431-4)

Information about building, installing and running MHM2 can be found in the [user guide](docs/mhm_guide.md)

MHM2 is developed and maintained by the [Exabiome Project](https://sites.google.com/lbl.gov/exabiome/) at Lawrence Berkeley National Laboratory, and is supported by the [Exascale Computing Project](https://www.exascaleproject.org) (17-SC-20-SC), a collaborative effort of the U.S. Department of Energy Office of Science and the National Nuclear Security Administration.
