# Motivation

Metahipmer2 (MHM2) is designed with the potential for multiple layers of parallelism. UPC++ enables PGAS parallelism, in which multiple processes can remotely access a shared global address space to execute a program in parallel. 

MHM2 also offers GPU-based SIMT parallelism if GPUs are available to the user. MHM2's GPU builds allow UPC++ processes to run certain routines on GPUs, resulting in significant performance gains (link to mhm2 gpu paper). The kmer analysis stage captured by this proxy app contains a number of these GPU-accelerated routines. Like many modern GPU codes, they are written in CUDA for Nvidia hardware (and HIP for AMD hardware).

[From the Kokkos abstract](https://kokkos.org/about/abstract/), "The Kokkos C++ Performance Portability Ecosystem is a production level solution for writing modern C++ applications in a hardware agnostic way... Kokkos Core is a programming model for parallel algorithms that use many-core chips and share memory among those cores." 

With the goal of a portable proxy app in mind, reimplementing vendor-specific GPU code was the first target of Kokkos integration. 

# Translating CUDA/HIP to Kokkos

## Project Structure Changes
Much of MetaHipmer's device code for the kmer analysis stage is located in the **kcount-gpu** directory. Kokkos is not a required dependency for this proxy app, so the pre-existing device code is still present in this directory for non-Kokkos GPU builds. 

(To keep new code distinct/isolated?) Kokkos was integrated through a corresponding? "mirror" directory, **kcount-kokkos**. The overall design and function of each compilation unit is preserved in kcount-kokkos; for example, parse_and_pack.cpp and parse_and_pack.hpp were rewritten as kokkos_pnp.cpp and kokkos_pnp.hpp. 

## GPU routines
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
<br>

Many `__device__` annotated functions were also "translated"; in most cases, this required little more than changing the function annotation to `KOKKOS_FUNCTION` or `KOKKOS_INLINE_FUNCTION`. 

The other key translation target was the data arrays undergoing parallelized computation. Kokkos uses a special array abstraction, called a View, to handle data management in parallel programs. Data arrays used in device code were translated to Kokkos Views. 

## Translation example

To demonstrate Kokkos translation on one of the simplest parallelized functions **gpu_unpack_supermer_block()**, see the code below. <br>

<div style="display: flex; ">

  <pre style="width: 50%; font-size: 12px; ">
<code>
<span style="background-color: #fff9c4;">__global__ void</span> gpu_unpack_supermer_block(SupermerBuff unpacked_supermer_buff, SupermerBuff packed_supermer_buff, int buff_len) {

  <span style="background-color: #fff9c4;">unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadid >= buff_len) return;</span>
  uint8_t packed = packed_supermer_buff.<span style="background-color: #fff9c4;">seqs[threadid]</span>;
  if (packed == '_') return;
  uint8_t left_side = (packed & 240) >> 4;
  unpacked_supermer_buff.seqs[threadid * 2] = to_base_func(left_side, packed);
  if (packed_supermer_buff.counts) unpacked_supermer_buff.counts[threadid * 2] = packed_supermer_buff.counts[threadid];
  uint8_t right_side = packed & 15;
  unpacked_supermer_buff.seqs[threadid * 2 + 1] = to_base_func(right_side, packed);
  if (packed_supermer_buff.counts) unpacked_supermer_buff.counts[threadid * 2 + 1] = packed_supermer_buff.counts[threadid];
}
</code>
  </pre>

  <pre style="width: 50%; font-size: 12px; padding: 10px;">
<code>
<span style="background-color: #fff9c4;">void</span> gpu_unpack_supermer_block(SupermerBuff unpacked_supermer_buff, SupermerBuff packed_supermer_buff, int buff_len) {

  <span style="background-color: #fff9c4;">Kokkos::parallel_for("gpu_unpack_supermer_block", buff_len - 1, KOKKOS_LAMBDA (int i) {</span>
    uint8_t packed = packed_supermer_buff.<span style="background-color: #fff9c4;">seqs_v(i)</span>;
    if (packed == '_') return;
    uint8_t left_side = (packed & 240) >> 4;
    unpacked_supermer_buff.seqs_v(i * 2) = to_base_func(left_side, packed);
    if (packed_supermer_buff.counts_v.size() > 1) unpacked_supermer_buff.counts_v(i * 2) = packed_supermer_buff.counts_v(i);
    uint8_t right_side = packed & 15;
    unpacked_supermer_buff.seqs_v(i * 2 + 1) = to_base_func(right_side, packed);
    if (packed_supermer_buff.counts_v.size() > 1) unpacked_supermer_buff.counts_v(i * 2 + 1) = packed_supermer_buff.counts_v(i);
  });  
}
</code>
  </pre>

</div>

The first change is the annotation from device code to host code. All Kokkos code can be thought of as host code which dispatches parallel work to an Execution Space- in this case, the GPU. 

### Parallel Dispatch

Parallel dispatch is launched by the Kokkos::parallel_for() function. The arguments passed to this function in the example above are straightforward:
- an optional name, useful for debugging and profiling with Kokkos tools
- an Execution Policy- this can effectively be thought of as an iterative range over which the parallel work occurs
- a functor which contains the body of parallel work. although there are other ways to implement the functor, a simple method is the with Kokkos's lambda macro KOKKOS_LAMBDA

So, the device code's "threadid" translates to the index parameter of the lambda function, which also takes the same place as threadid in the parallel code body. Instead of setting an upper limit on the threadid, the iteration range is moved to the corresponding parallel_for() param. 

For more information on parallel dispatch, see the (link to kokkos docs)

### Views

The other change highlighted in this example is the use of a Kokkos View for data management. Views are array abstractions which can be allocated in different Execution Spaces, such as a GPU. In the original device code, data from the host must be copied to the device with a function like cudaMemcpy(). In Kokkos the data still moves from host to device, but this is done via Kokkos Views and the function deep_copy().

In the example above, the SupermerBuff is a struct with a member seqs, a pointer to a char array representing a DNA sequence. In the original code, this array is first copied to a char array on the device. <br>
```
struct SupermerBuff {
  char *seqs;
  count_t *counts;
};
```
```
ERROR_CHECK(Memcpy(packed_elem_buff_dev.seqs, elem_buff_host.seqs, buff_len, MemcpyHostToDevice));
```

In the Kokkos translation, those arrays are contained in Views. Views are allocated in the default execution space, in this case CudaSpace, unless specified otherwise. Also, in order to copy between Views, they must have the same memory layout, easily accomplished with the HostMirror typedef. 
```
struct SupermerBuff {
  Kokkos::View<char*> seqs_v;
  Kokkos::View<char*>::HostMirror h_seqs_v;

  Kokkos::View<count_t*> counts_v;
  Kokkos::View<count_t*>::HostMirror h_counts_v;
};
```
```
Kokkos::deep_copy(packed_elem_buff_dev.seqs_v, elem_buff_host.h_seqs_v);
```
Views are indexed like Fortran arrays, hence the data_v(i) syntax, seen in the example above. 

For more information on Views, see the (link to kokkos docs)


## Challenges

- device constants
- note on TCF

# Future Kokkos Integration

- propose reimplementing other assembly stages' GPU code with kokkos (alignment, local assembly)

- propose replacing UPC++ PGAS layer with Kokkos Remote Spaces to achieve "maximum portability"