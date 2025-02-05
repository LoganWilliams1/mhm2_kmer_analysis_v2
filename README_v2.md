# MHM2 KMER ANALYSIS PROXY APP
MetaHipmer2 (MHM2) is an exascale metagenome assembler software.

This application encapsulates the kmer analysis stage of the [MHM2 pipeline](https://www.nature.com/articles/s41598-020-67416-5/figures/7). 

Kmers are short DNA strings of a fixed length *k*, typically derived from DNA sequencing data; they are used in many bioinformatics applications, including genome assembly. 

For a given set of FASTQ input files and desired kmer length *k*, MKAPA (?) builds a distributed hash table of kmers, where keys are unique kmers and values include counts of that kmer as well as its most probable single base extensions.

(small example table?)
(output?)
(overall goals?)

For more MHM2 information and source code, see the [JGI home page](https://jgi.doe.gov/data-and-tools/software-tools/metahipmer/) and the [MHM2 Guide.](https://bitbucket.org/berkeleylab/mhm2/src/660ae5ab13f6a2c576b17d8c9925d44cffb98c6d/docs/mhm_guide.md)

<br>

## Install Dependencies

**1. Install UPC++ [Required]** <br>
MHM2 is built on UPC++, a Partitioned Global Address (PGAS) programming library. 
- Find the latest release at the [UPC++ home page.](https://bitbucket.org/berkeleylab/upcxx/wiki/Home)
- To install, see [UPC++ INSTALL.md](https://bitbucket.org/berkeleylab/upcxx/wiki/INSTALL)
- For more info, see [UPC++ README.md](https://bitbucket.org/berkeleylab/upcxx/wiki/README.md)
- For additional build recipes and information, see (wiki?) 

<br>

**2. Install Kokkos [Recommended]** <br>
Kokkos is a parallel programming library for performance portability. This application offers a Kokkos implementation of GPU accelerated routines, enabling interoperability with UPC++ host code. 
- Find the latest release at the [Kokkos release page.](https://github.com/kokkos/kokkos/releases)
- To install, see [Kokkos Quick Start Guide](https://kokkos.org/kokkos-core-wiki/quick_start.html)
- For more info, see the [Kokkos Website](https://kokkos.org/)

<br>

## Install the Kmer Analysis App
**1. Clone the repository:** <br>
`git clone https://github.com/LoganWilliams1/mhm2_kmer_analysis_v2.git`

**2. Navigate into project directory:** <br>
`cd mhm2_kmer_analysis_v2`

**3. Build and install:** <br>
Currently, the recommended method is to run the provided script `build.sh`. This is the same script included with the parent MHM2, except for a small change to allow extra cmake options. This script installs the binary `mhm2` (change?) into the `install/bin` subdirectory.

A basic build command is: <br>
`./build.sh Release`

To build the Kokkos version, make sure Kokkos is installed and 
run: <br>
`./build.sh Release -DENABLE_KOKKOS=ON`

Note: [UPC++ installed on Infiniband systems will use MPI for spawning](https://bitbucket.org/berkeleylab/upcxx/wiki/INSTALL#markdown-header-configuration-linux), and it may be necessary to specify an MPI compiler wrapper to build correctly. Ensure that this app (and Kokkos, if installed) are built with the same compiler. For example, you may need to include this option to build: <br>
```./build.sh Release -DENABLE_KOKKOS=ON -DCMAKE_CXX_COMPILER=`which mpicxx` ```

<br>

## Run the Kmer Analysis App
A simple python wrapper script is provided for convenience. Running the app *requires* two options to be specified: the number of UPC++ processes, and a FASTQ input. Therefore a minimal run command is: <br>
`./run_app.py -n 1 -p seqs_1.fastq seqs_2.fastq` 

<br>

### Run Options

**`-n INT`** <br>
Number of UPC++ processes. This option is required. 

**`-N INT`** <br>
Number of nodes. The default is 1.

**`-p STRING STRING ...`** <br>
Files containing separate paired reads. The first of paired files must be immediately followed by the second. For example, to count kmers in two libraries of separate paired reads, the option should be specified as: <br>
`-p lib1_1.fastq lib1_2.fastq lib2_1.fastq lib2_2.fastq`

**`-r STRING STRING ...`** <br>
Files containing interleaved paired reads.

**`-u STRING STRING...`** <br>
Files containing unpaired reads.

**`-k INT`** <br>
Length of kmers to analyze. The default is 21.

**`-o STRING`** <br>
Name of the output directory. The default is: <br>
`mhm2-run-<READS_FNAME1>-n<PROCS>-N<NODES>-YYMMDDhhmmss-<JOBID>`

**`--use-qf=BOOL`** <br>
This option is only applicable for **non-Kokkos GPU builds**. Setting it to true optimizes memory usage on GPUs, although it does result in a negligible amount of variation in kmer count results. See the [implementation notes]() for more details.



## Sample Run
Two small FASTQ files and a bash script are provided as a basic test case. The script runs the app with 4 UPC++ processes, counts 21-mers, and reports if the run was successful. Additional options aren't necessary for this script:

`./run_sample.sh`