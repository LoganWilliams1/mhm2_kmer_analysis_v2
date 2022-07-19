# ChangeLog.md for MetaHipMer version 2 (aka MHM2 or mhm2)


This is the ChangeLog for MetaHipMer with development at [bitbucket](https://bitbucket.org/berkeleylab/mhm2)

### 2.1.0 2021-01-24
   * Major improvements in GPU performance in kcount and adept
   * Added quotient filtering to GPU kcount
   * Improved Continuous Integration testsing
   * Optionally trim read adapters to improve load balance and assembly quality
   * Fixed issues #108, #112, #109, #107, #83, #97, #91, #89, #39, #10, #84

### 2.0.2 2022-01-24
   * Major updates to GPU offloading in kcount - PR #48
      * Removed bloom filters from CPU & GPU kcount
   * Counting quotient filter in kcount (GPU only) - PR #54
   * Fixed potential deadlock in kcount atomics PR #55, Issue #108
   * Read shuffling based on minimizer kmers, not alignments 
   * Improved GPU autodetection and initialization
   * Trimming of read adapters
   * Other PR #40, #42, #44, #45, #46, #49, #51, #53, #55
   * Minor bugfixes Issue #97, #112, #86, #63, #109, #107, #83, #91, #92, #25, #37, #89
   * Fixed mhm2.py script for PBS - Issue #32

### 2.0.1.2 2021-04-22
   * A very minor bugfix release to facilitate building against UPC++ 2021.3.0

### 2.0.1 2020-12-17
   * Use minimizer hashes to improve locality of kmers and greatly speed up deBruijn graph traversal - pull request #18
   * Aggregate messages in LocalAssem reducing overall time other optimizations in preparation for GPU acceleration - pull request #17, #16
   * Fixed scaling Issues - #52, #54, #61, #72
   * Fixed cmake/build - Issues #20, #51, #55, #68, #67
   * Fixed paired read usage - Issue #58  
   * Fixed minor UI Issues - #53, #33, #50, #46, #45, #80, #78, #77, #79
   * Fixed edge cases of improper lustre striping - Issue #62, #70
   * Added framework for unit tests
   * Added a ability to build a docker image
   * Added better support for paired reads - Issue #65, #64
   * Added better support for unpaired read - Issue #59
   * Better post-assembly sam and depths file - Issue #71
   * Corrected CIGAR alignments in some cases

### 2.0.0 2020-09-30
   * Updated documentation and added a user guide.
   * Incorporated a 3 tier aggregating store to allow scaling on 1000s of nodes Issue #4
   * Fixed stall on summit with over 200 node - Issue #34
   * Improved reading of multiple input files - Issue #38
   * Greatly improved CI testing, including quality tests - Issue #22
   * Fixed error when writing to a lagging filesystem - Issues #26 #15
   * Various bug fixes - Issues #27 #29 #31 #36 #44
   * Included build support for MacOSX (not tested)
   * Various cosmetic - Issues #38 #40 #42 #47 #43
      * Renamed from MHMXX to MHM2 - Issue #18
      * Reformatted the entire code base

### 0.1.3 2020-08-20
   * Added support for GPU in alignments with CUDA-enabled detection and build
      * major refactor of adept-sw for better performance and ease of building
   * Fixed fastq reader edge cases in record boundary detection
   * Support in options for 1 contig and 1 scaffolding round
   * Optimal execution environment on Summit with custom spawner
   * Better checkpointing support
   * Various fixes from upcxx-utils including ofstream bugs

### 0.1.2 2020-07-06
   * Semi-stable alpha release
   * Added skeleton for CI building on cori and hulk development server
   * Fixed re-build issues with upcxx-utils
   * Added warning on build and install when submodule is out of date
   * Various fixes for Summit build including vector support for SSW alignment
   * Added support for more cluster schedulers
   * Added --paired-reads support and modified Fastq class to support two-file libraries
   * Fixed LUSTRE striping on per_thread rank logs


### 0.1.1 Before 2020-07-01
   * Complete rewrite of MetaHipMer 1.0 using only UPC++
