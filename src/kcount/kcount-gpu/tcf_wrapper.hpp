#pragma once

/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/

// inclusions to build the TCF
#include <poggers/metadata.cuh>
#include <poggers/hash_schemes/murmurhash.cuh>
#include <poggers/probing_schemes/double_hashing.cuh>
#include <poggers/probing_schemes/power_of_two.cuh>

// new container for 2-byte key val pairs
#include <poggers/representations/grouped_key_val_pair.cuh>

#include <poggers/representations/key_val_pair.cuh>
#include <poggers/representations/dynamic_container.cuh>

#include <poggers/sizing/default_sizing.cuh>

#include <poggers/insert_schemes/power_of_n_shortcut.cuh>

#include <poggers/insert_schemes/power_of_n_shortcut_buckets.cuh>

#include <poggers/representations/packed_bucket.cuh>

#include <poggers/insert_schemes/linear_insert_buckets.cuh>

#include <poggers/tables/bucketed_table.cuh>

#include <poggers/representations/grouped_storage_sub_bits.cuh>

#include "tcf_hash_wrapper.hpp"

namespace two_choice_filter {

__constant__ char kmer_ext[6] = {'F', 'A', 'C', 'T', 'G', '0'};

#define TCF_SMALL 1

#define TCF_DELETE 0

#if TCF_SMALL

using backing_table = poggers::tables::bucketed_table<
    uint64_t, uint8_t,
    poggers::representations::dynamic_bucket_container<poggers::representations::dynamic_container<
        poggers::representations::bit_grouped_container<10, 6>::representation, uint16_t>::representation>::representation,
    1, 8, poggers::insert_schemes::linear_insert_bucket_scheme, 20, poggers::probing_schemes::doubleHasher,
    poggers::hashers::murmurHasher>;
using TCF = poggers::tables::bucketed_table<
    uint64_t, uint8_t,
    poggers::representations::dynamic_bucket_container<poggers::representations::dynamic_container<
        poggers::representations::bit_grouped_container<10, 6>::representation, uint16_t>::representation>::representation,
    1, 8, poggers::insert_schemes::power_of_n_insert_shortcut_bucket_scheme, 2, poggers::probing_schemes::doubleHasher,
    poggers::hashers::murmurHasher, true, backing_table>;
#else

using backing_table =
    poggers::tables::bucketed_table<uint64_t, uint16_t,
                                    poggers::representations::dynamic_bucket_container<poggers::representations::dynamic_container<
                                        poggers::representations::key_val_pair, uint16_t>::representation>::representation,
                                    1, 8, poggers::insert_schemes::linear_insert_bucket_scheme, 20,
                                    poggers::probing_schemes::doubleHasher, poggers::hashers::mhm_hasher>;
using TCF =
    poggers::tables::bucketed_table<uint64_t, uint16_t,
                                    poggers::representations::dynamic_bucket_container<poggers::representations::dynamic_container<
                                        poggers::representations::key_val_pair, uint16_t>::representation>::representation,
                                    1, 8, poggers::insert_schemes::power_of_n_insert_shortcut_bucket_scheme, 2,
                                    poggers::probing_schemes::doubleHasher, poggers::hashers::mhm_hasher, true, backing_table>;

#endif

__device__ uint16_t pack_extensions(char left, char right) {
  uint16_t ret_val = 0;

  for (uint i = 0; i < 6; i++) {
    if (left == kmer_ext[i]) {
      ret_val += i << 8;
    }

    if (right == kmer_ext[i]) {
      ret_val += i;
    }
  }

  return ret_val;
}

__device__ bool unpack_extensions(uint16_t storage, char& left, char& right) {
  uint16_t left_val = ((storage & 0xff00) >> 8);

  uint16_t right_val = (storage & 0x00ff);

  if ((left_val < 6) && (right_val < 6)) {
    left = kmer_ext[left_val];
    right = kmer_ext[right_val];

    return true;
  } else {
    return false;
  }
}

__device__ uint8_t pack_extensions_small(char left, char right) {
  uint8_t ret_val = 0;

  for (uint i = 0; i < 6; i++) {
    if (left == kmer_ext[i]) {
      ret_val += i << 3;
    }

    if (right == kmer_ext[i]) {
      ret_val += i;
    }
  }

  return ret_val;
}

__device__ bool unpack_extensions_small(uint8_t storage, char& left, char& right) {
  // grab bits 3-5
  uint8_t left_val = ((storage & 0x38) >> 3);

  // grab bits 0-2
  uint8_t right_val = (storage & 0x07);

  if ((left_val < 6) && (right_val < 6)) {
    left = kmer_ext[left_val];
    right = kmer_ext[right_val];

    return true;
  } else {
    return false;
  }
}

// returns the usage of the TCF
// this should be accurate to within a few bytes
// only deviation from the formula is based off the difference in block size.
__host__ __device__ uint64_t estimate_memory(uint64_t max_num_kmers) {
// estimate to 120% to be safe
// and 4 bytes per item pair
#if TCF_SMALL

  return (max_num_kmers * 2 * 1.2);

#else

  return (max_num_kmers * 4 * 1.2);

#endif
}

__host__ poggers::sizing::size_in_num_slots<2> get_tcf_sizing(uint64_t max_num_kmers) {
  uint64_t max_slots = max_num_kmers * 1.2;

  // 90/11 split over size for forward and backing tables.

  poggers::sizing::size_in_num_slots<2> my_size(max_slots * .9, max_slots * .1);
  return my_size;
}

__host__ poggers::sizing::size_in_num_slots<2> get_tcf_sizing_from_mem(uint64_t available_bytes) {
#if TCF_SMALL

  uint64_t max_slots = available_bytes / 2;

#else

  uint64_t max_slots = available_bytes / 4;

#endif

  // 90/11 split over size for forward and backing tables.
  poggers::sizing::size_in_num_slots<2> my_size((max_slots * 90ULL / 100ULL), (max_slots * 10ULL / 100ULL));
  // printf("%llu bytes gives %llu slots, given: %llu + %llu = %llu\n", available_bytes, max_slots, my_size.next(), my_size.next(),
  // my_size.total());
  my_size.reset();
  return my_size;
}

}  // namespace two_choice_filter
