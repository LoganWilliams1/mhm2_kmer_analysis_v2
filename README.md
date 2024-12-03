# MetaHipMer2 Proxy Apps

Performance-portable proxy applications derived from [MetaHipMer2 (de novo metagenome assembler)](https://bitbucket.org/berkeleylab/mhm2/src/master/)

## [mhm2-kmer-analysis](github project here)

```
git clone  tbd

```
## Build

```
env tbd

```


```
TBD
```

# Required Third-Party Library 

## [UPC++:  a PGAS library for C++](https://bitbucket.org/berkeleylab/upcxx/src/master/)

  - [UPC++ basic build recipe](https://bitbucket.org/berkeleylab/upcxx/src/master/INSTALL.md)
  - [UPC++ Developers Guide](https://bytebucket.org/berkeleylab/upcxx/wiki/docs/guide.pdf?rev=767d43b34dd00f6b2765a49d6ebf55bee91f4579)
  - [UPC++ Specification](https://bytebucket.org/berkeleylab/upcxx/wiki/docs/spec.pdf?rev=767d43b34dd00f6b2765a49d6ebf55bee91f4579)

**Nota bene:  [UPC++ uses GASNet for PGAS communication primitives](https://bitbucket.org/berkeleylab/gasnet/src/stable/)**


### UPC++ for Supercomputers

#### Architectures

  - [Power9/NVIDIA V100](https://hpc.llnl.gov/hardware/compute-platforms/lassen)

```
../configure \
--with-cxx=`which mpicxx` \
--with-cc=`which mpicc` \
--enable-udp \
--enable-mpi \
--enable-cuda \
--enable-ibv \
--with-cxxflags=-std=c++17 \
--prefix=${INSTALL_DIR}


``

  -  [MI300A](https://hpc.llnl.gov/hardware/compute-platforms/rzadams)

```
# for rocm / cray env 
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

```


```
../configure \
--with-cxx=mpicxx \
--with-cc=mpicc \
--with-pmi-rpath \
--with-ofi-spawner=pmi \
--enable-udp \
--enable-hip \
--enable-mpi \
--enable-ofi \
--with-cxxflags=-std=c++17 \
--with-default-network=ofi \
--with-ofi-provider=cxi \
--prefix=${INSTALL_DIR}

```


## Kokkos

  - [Kokkos Quick Start](https://kokkos.org/kokkos-core-wiki/quick_start.html) 


# Running

# Inputs

## Synthetic Reads to Benchmark

  - [ArcticSynth Simulated Sequences](https://portal.nersc.gov/project/hipmer/MetaHipMer_datasets_12_2019/ArcticSynth/)
  - [Detailed in, "Terabase-scale metagenome coassembly with MetaHipMer"](https://www.nature.com/articles/s41598-020-67416-5) and in [Supplemental Info](https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-67416-5/MediaObjects/41598_2020_67416_MOESM1_ESM.pdf)

## NCBI Short-Read Archive

  - [Gut Microbiota in Long-Lived Humans, 2020](https://pubmed.ncbi.nlm.nih.gov/32209716/)
    - Project:  PRJNA553191 
    - Accession:  [SRR9654970](https://www.ncbi.nlm.nih.gov/sra/?term=SRR9654970)
    - Dataset size:  1.8 Gb

  - [Ovine Gut Microbes](https://pubmed.ncbi.nlm.nih.gov/24170804/)
    - Project:  PRJNA214227
    - Accession:  [SRX331704](https://www.ncbi.nlm.nih.gov/sra/?term=SRX331704)

  - [Terabase-scale metagenome coassembly with MetaHipMer](https://pubmed.ncbi.nlm.nih.gov/32612216/)
    - Project:  [PRJNA367436](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA367436)
    - Accession:  [SRR5819383](https://www.ncbi.nlm.nih.gov/sra/SRR5819383)
    - Dataset size:  15.9 Gb 

## NCBI Tools

  - Install [SRA  Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to interact with ["Short Read Archive" (SRA)](https://www.ncbi.nlm.nih.gov/sra)



# Further Reading 

## MetaHipMer2 - MHM2 (Parent App)
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
