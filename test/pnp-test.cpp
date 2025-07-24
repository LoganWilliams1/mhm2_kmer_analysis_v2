#include "kokkos_pnp.hpp"
#include <Kokkos_Core.hpp>

#include "gtest/gtest.h"


using namespace kcount_gpu;

class Kokkos_PnP_Test : public testing::Test {
  protected:

    ParseAndPackGPUDriver *driver;	  
    void SetUp() override {
      double init_time = 0.0;

      driver = new ParseAndPackGPUDriver(0, 1, 33, 21, 1, 15, init_time);
      // make seqs
      // add to h_seqs_v
      // deep copy to driver->dev_seqs_v
      //
      //    
    }

    void TearDown() override {
      delete driver;
    }
};
  
TEST_F(Kokkos_PnP_Test, parse_and_pack) {


  //parse_and_pack(driver->dev_seqs_v, int, int, int, int, driver->dev_kmer_targets_v, int, driver->v, driver->v	

  EXPECT_EQ(0,1) << "temp parse_and_pack";  

  
}




TEST_F(Kokkos_PnP_Test, build_supermers) {
  
  // make h_kmer_targets_v
  // deep copy to driver->dev_kmer_targets_v
  // make h_num_supermers_v = 0
  // deep copy to dev_num_supermers_v
  //
  // build_supermers(driver->dev_seqs_v, driver->v, int, int, int, driver->v, driver->v, int, int)

  EXPECT_EQ(0,1) << "temp build_sup";
}


TEST_F(Kokkos_PnP_Test, pack_seqs) {
  
  // pack_seqs(driver->v, driver->v, int)

  EXPECT_EQ(0,1) << "temp pack_seqs";
}
