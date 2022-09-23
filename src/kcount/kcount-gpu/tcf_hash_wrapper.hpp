#ifndef POGGERS_HASH_WRAPPER 
#define POGGERS_HASH_WRAPPER


#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <assert.h>

#include <cooperative_groups.h>

#include "gpu_hash_funcs.cpp"

namespace cg = cooperative_groups;


namespace poggers {

namespace hashers {




//hashers should always hash an arbitrary type to uint64_t
template <typename Key, std::size_t Partition_size = 1>
struct __attribute__ ((__packed__)) mhm_hasher {


	//tag bits change based on the #of bytes allocated per block
private:

public:

	using key_type = Key;
	//using internal_paritition_size = Partition_size;
	
	//typedef key_val_pair<Key> Key;

	//init happens by a single thread on CPU/GPU
	//no cg needed

	__host__ __device__ void init(uint64_t ext_seed){	
	}


	__host__ __device__ uint64_t hash(Key key_to_hash){

		Key copy_key = key_to_hash;
		
		return gpu_murmurhash3_64(&copy_key, sizeof(Key));


	}

	//all participate
	__device__ uint64_t hash(Key key_to_hash, cg::thread_block_tile<Partition_size> group){

		Key copy_key = key_to_hash;

		return gpu_murmurhash3_64(&copy_key, sizeof(Key));

	}


	//no need for explicit destructor struct hash no memory components

};

}

}


#endif //GPU_BLOCK_