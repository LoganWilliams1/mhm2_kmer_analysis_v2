# Motivation

Metahipmer2 (MHM2) is designed with the potential for multiple layers of parallelism. UPC++ enables PGAS parallelism, in which multiple processes can remotely access a shared global address space to execute a program in parallel. 

MHM2 also offers GPU-based SIMT parallelism if GPUs are available to the user. MHM2's GPU builds allow UPC++ processes to run certain routines on GPUs, resulting in significant performance gains (link to mhm2 gpu paper). The kmer analysis stage captured by this proxy app contains a number of these GPU-accelerated routines. Like many modern GPU codes, they are written in CUDA for Nvidia hardware (and HIP for AMD hardware).

[From the Kokkos abstract](https://kokkos.org/about/abstract/), "The Kokkos C++ Performance Portability Ecosystem is a production level solution for writing modern C++ applications in a hardware agnostic way... Kokkos Core is a programming model for parallel algorithms that use many-core chips and share memory among those cores." 

With the goal of a portable proxy app in mind, reimplementing vendor-specific GPU code was the first target of Kokkos integration. 

# Translating CUDA/HIP to Kokkos

- explain kcount-gpu vs kcount-kokkos
- list specific functions, source files
- show real code example of CUDA vs Kokkos kernel function
- explain how kokkos build does not need gpu-utils (does it..?)

### Project Structure Changes
Much of MetaHipmer's device code for the kmer analysis stage is located in the **kcount-gpu** directory. Kokkos is not a required dependency for this proxy app, so the pre-existing device code is still present in this directory for non-Kokkos GPU builds. 

(To keep new code distinct/isolated?) Kokkos was integrated through a corresponding? "mirror" directory, **kcount-kokkos**. The overall design and function of each compilation unit is preserved in kcount-kokkos; for example, parse_and_pack.cpp and parse_and_pack.hpp were rewritten as kokkos_pnp.cpp and kokkos_pnp.hpp. 

### GPU routines
Translating CUDA/HIP code to the Kokkos programming model was relatively straightforward for much of this implementation. The most fundamental concept was writing GPU kernels as Kokkos parallel dispatch operations- these routines served as the starting point around which the rest of the Kokkos code was written. 

The full list of functions containing parallelized code is below.

| source | function |
|---|---|
| parse_and_pack.cpp | parse_and_pack() |
| | build_supermers() |
| | pack_seqs() |
| gpu_hash_table.cpp | gpu_compact_ht() |
| | gpu_purge_invalid() |
| | gpu_unpack_supermer_block() |
| | gpu_insert_supermer_block() |

- note about other device functions
- note about data management and kokkos views

### Translation example

- pick function, walk through translation, screenshots?
- show example of data structure translated to kokkos view

# Future Kokkos Integration

- propose reimplementing other assembly stages' GPU code with kokkos (alignment, local assembly)

- propose replacing UPC++ PGAS layer with Kokkos Remote Spaces to achieve "maximum portability"