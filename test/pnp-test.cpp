#include "kokkos_pnp.hpp"
#include <Kokkos_Core.hpp>

#include "gtest/gtest.h"


using namespace kcount_gpu;
using namespace std;


void parse_and_pack(Kokkos::View<char *> seqs, int minimizer_len, int kmer_len, int num_longs, int seqs_len,
                    Kokkos::View<int *> kmer_targets_v, int num_ranks, Kokkos::View<uint64_t[256]> twins_v,
                    Kokkos::View<uint64_t[32]> mask_v);

void build_supermers(Kokkos::View<char *> seqs_v, Kokkos::View<int *> kmer_targets_v, unsigned int num_kmers, int kmer_len,
                     int seqs_len, Kokkos::View<kcount_gpu::SupermerInfo *> supermers_v, Kokkos::View<unsigned int> num_supermers_v,
                     unsigned int &num_valid_kmers, int rank_me);

void pack_seqs(Kokkos::View<char *> dev_seqs_v, Kokkos::View<char *> dev_packed_seqs_v, int seqs_len);


class Kokkos_PnP_Test : public testing::Test {
  protected:

    ParseAndPackGPUDriver *driver;	
    int upcxx_rank_me = 0;
    int upcxx_rank_n = 4;
    int qual_offset = 33;
    int kmer_len = 21;
    int num_kmer_longs = 1;
    int minimizer_len = 15;
    double init_time = 0.0;

    // two random reads, each with length = 50 bases
    // joined with '_' as a sequence block
    string seq_block = "TCTTGCCGAGACCCAGGTTCATCCCATAGCAGGGTTTGATAGTAGGCTCC_GTGATACTACTAGATCCACTATAGTACCCGTCCTTCGACCGCCGACAAAC";

    // number of actual kmers = length of read - length of kmer + 1, so 30 kmers for each read
    // number of kmer targets = length of entire sequence block - length of kmer + 1
    int num_kmers = 60;
    int num_kmer_targets = 81;

    // rank of the upc++ process (0-3) to which each kmer will be sent
    // when the kmer window includes the read separator '_', that kmer is invalid and gets a targe of -1
    vector<int> expected_kmer_targets = {2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 0, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0};

    // the driver does not have a host kmer targets View- this is only useful for testing 
    Kokkos::View<int *>::HostMirror h_kmer_targets_v;



    /**
     * instantiate new gpu driver for each kernel test
     * copy sequence block to device
     * instantiate host kmer targets View
     */
    void SetUp() override {

      driver = new ParseAndPackGPUDriver(upcxx_rank_me, upcxx_rank_n, qual_offset, kmer_len,  
                                         num_kmer_longs, minimizer_len, init_time);

      for (int i = 0; i < seq_block.length(); i++) {
        driver->h_seqs_v(i) = seq_block[i];
      }
      Kokkos::deep_copy(driver->dev_seqs_v, driver->h_seqs_v);


      h_kmer_targets_v = create_mirror_view(driver->dev_kmer_targets_v);


    }

    void TearDown() override {
      delete driver;
    }
};
  

/**
 * parse_and_pack computes a target rank for each kmer in the sequence block
 */
TEST_F(Kokkos_PnP_Test, parse_and_pack) {

  parse_and_pack(driver->dev_seqs_v, minimizer_len, kmer_len, num_kmer_longs, seq_block.length(),
                 driver->dev_kmer_targets_v, upcxx_rank_n, driver->twins_v, driver->mask_v);	

  deep_copy(h_kmer_targets_v, driver->dev_kmer_targets_v);
  vector<int> test_kmer_targets;
  for (int i = 0; i < num_kmer_targets; i++) {
    test_kmer_targets.push_back(h_kmer_targets_v(i));
  }
   
  EXPECT_EQ(expected_kmer_targets, test_kmer_targets) << "parse_and_pack result: kmer_targets";

}



/**
 * build_supermers combines adjacent kmers with shared targets into "supermers" as a message passing optimization
 * supermers are stored as SupermerInfo structs that contain the target rank, offset (location in the block), and length of the combined sequence
 * this kernel also computes the number of valid kmers in the block
 */
TEST_F(Kokkos_PnP_Test, build_supermers) {
  
  
  for (int i = 0; i < num_kmer_targets; i++) {
    h_kmer_targets_v(i) = expected_kmer_targets[i];
  }
  deep_copy(driver->dev_kmer_targets_v, h_kmer_targets_v);
  driver->h_num_supermers_v() = 0;
  Kokkos::deep_copy(driver->dev_num_supermers_v, driver->h_num_supermers_v);
  unsigned int num_valid_kmers = 0;
  build_supermers(driver->dev_seqs_v, driver->dev_kmer_targets_v, num_kmer_targets, kmer_len, seq_block.length(), driver->dev_supermers_v, driver->dev_num_supermers_v,
                  num_valid_kmers, upcxx_rank_me);

  Kokkos::deep_copy(driver->h_num_supermers_v, driver->dev_num_supermers_v);
  Kokkos::deep_copy(driver->h_supermers_v, driver->dev_supermers_v);
  Kokkos::fence();

  vector<int> test_supermer_targets;
  vector<int> test_supermer_offsets;
  vector<int> test_supermer_lens;

  vector<int> expected_supermer_targets = {1, 3, 1, 0, 3, 0, 0, 2, 0, 2, 3, 0};
  vector<int> expected_supermer_offsets = {51, 55, 58, 65, 67, 74, 0, 1, 9, 16, 17, 19};
  vector<int> expected_supermer_lens = {26, 25, 29, 24, 29, 27, 23, 30, 29, 23, 24, 31};
  int num_supermers = 12;
  
  for (int i = 0; i < driver->h_num_supermers_v(); i++) {
    test_supermer_targets.push_back(driver->h_supermers_v(i).target);
    test_supermer_offsets.push_back(driver->h_supermers_v(i).offset);
    test_supermer_lens.push_back(driver->h_supermers_v(i).len);
  }
  
  EXPECT_EQ(num_kmers, num_valid_kmers) << "build_supermers result: num valid kmers";
  EXPECT_EQ(num_supermers, driver->h_num_supermers_v()) << "build_supermers result: num supermers";
  EXPECT_EQ(expected_supermer_targets, test_supermer_targets) << "build_supermers result: supermer targets";
  EXPECT_EQ(expected_supermer_offsets, test_supermer_offsets) << "build_supermers result: supermer offsets";
  EXPECT_EQ(expected_supermer_lens, test_supermer_lens) << "build_supermers result: supermer lens";
}


/**
 * pack_seqs compresses the seq_block to half the original size
 * some bases are converted to control characters such as '\x86' so this test converts the packed seqs to their raw byte values for straightforward comparison
 */
TEST_F(Kokkos_PnP_Test, pack_seqs) {
  
  pack_seqs(driver->dev_seqs_v, driver->dev_packed_seqs_v, seq_block.length());
  Kokkos::deep_copy(driver->h_packed_seqs_v, driver->dev_packed_seqs_v);


  int packed_seqs_len = (seq_block.length() + 1) / 2;
  vector<unsigned int> test_packed_seqs;
  vector<unsigned int> expected_packed_seqs = {134, 136, 118, 103, 87, 86, 102, 87, 120, 134, 88, 102, 101, 133, 118, 87, 119, 136, 135, 88, 87, 133, 119, 104, 102, 7, 135, 88, 86, 133, 104, 87, 88, 102, 86, 133, 133, 120, 86, 102, 120, 102, 136, 103, 86, 103, 102, 117, 101, 85, 96};

  for (int i = 0; i < packed_seqs_len; i++) {
    test_packed_seqs.push_back(static_cast<unsigned int>(static_cast<unsigned char>(driver->h_packed_seqs_v(i))));
  }

  EXPECT_EQ(expected_packed_seqs, test_packed_seqs) << "pack_seqs result: packed_seqs";
}

