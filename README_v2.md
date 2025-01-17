# MHM2 KMER ANALYSIS PROXY APP
MetaHipmer2 (MHM2) is an exascale metagenome assembler software.

This application encapsulates the kmer analysis stage of the [MHM2 pipeline](https://www.nature.com/articles/s41598-020-67416-5/figures/7). 

Kmers are short DNA strings of a fixed length *k*, typically derived from DNA sequencing data; they are used in many bioinformatics applications, including genome assembly. 

For a given set of FASTQ input files and desired kmer length *k*, MKAPA (?) builds a distributed hash table of kmers, where keys are unique kmers and values include counts of that kmer as well as its most probable single base extensions.

(small example table?)
(output?)
(overall goals?)

For more MHM2 information and source code, see the [JGI home page](https://jgi.doe.gov/data-and-tools/software-tools/metahipmer/) and the [MHM2 Guide.](https://bitbucket.org/berkeleylab/mhm2/src/660ae5ab13f6a2c576b17d8c9925d44cffb98c6d/docs/mhm_guide.md)

## Install dependencies

**1. Install UPC++ [Required]** <br>
MHM2 is built on UPC++, a Partitioned Global Address (PGAS) programming library. 
- Find the latest release at the [UPC++ home page.](https://bitbucket.org/berkeleylab/upcxx/wiki/Home)
- To install, see [UPC++ INSTALL.md](https://bitbucket.org/berkeleylab/upcxx/wiki/INSTALL)
- For more info, see [UPC++ README.md](https://bitbucket.org/berkeleylab/upcxx/wiki/README.md)
- For additional build recipes and information, see (wiki?) 
<br><br>

**2. Install Kokkos [Recommended]** <br>
Kokkos is a parallel programming library for performance portability. This application offers a Kokkos implementation of GPU accelerated routines, enabling interoperability with UPC++ host code. 
- Find the latest release at the [Kokkos release page.](https://github.com/kokkos/kokkos/releases)
- To install, see [Kokkos Quick Start Guide](https://kokkos.org/kokkos-core-wiki/quick_start.html)
- For more info, see the [Kokkos Website](https://kokkos.org/)
<br><br>

## Install the Kmer Analysis App
**1. Clone the repository:** <br>
`git clone https://github.com/LoganWilliams1/mhm2_kmer_analysis_v2.git`

**2. Navigate into project directory:** <br>
`cd mhm2_kmer_analysis_v2`

**3. Run the build script or cmake?**

### Run the Kmer Analysis App
mhm2 compare script? python script?
which options to keep?

A basic run command *requires* two inputs: FASTQ input and kmer length.

`upcxx-run ./mhm2 -p seqs_1.fastq seqs_2.fastq -k 21`

Use the -k flag to denote the desired kmer length.

Use one of the following flags for FASTQ input: -p, -r, -u.
- -p denotes files containing separate paired reads. The first of paired files must be immediately followed by the second. For example, to count kmers in two libraries of separate paired reads, the option should be specified as: <br>
`-p lib1_1.fastq lib1_2.fastq lib2_1.fastq lib2_2.fastq`

- -r denotes files containing interleaved paired reads

- -u denotes files containing unpaired reads

Multiple files must be comma or space separated.







