// kcount - kmer counting
// Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

#include "utils.hpp"
#include "kcount.hpp"
#include "progressbar.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;

static const int QUAL_CUTOFF = 20;

extern ofstream _dbgstream;

uint64_t estimate_cardinality(unsigned kmer_len, vector<string> reads_fname_list)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_kmers = 0;
  string read_record[4];
  int64_t estimated_total_records = 0;
  int64_t total_records_processed = 0;
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    int64_t records_processed = 0;
    bool isCompressed = hasEnding(merged_reads_fname, ".gz");
    int64_t fileSize = get_file_size(merged_reads_fname);
    zstr::ifstream reads_file(merged_reads_fname);
    int64_t bytes_read = 0;
    ProgressBar progbar(merged_reads_fname, &reads_file, "Scanning reads file to estimate cardinality");
    while (!reads_file.eof()) {
      bool done = false;
      for (int i = 0; i < 4; i++) {
        getline(reads_file, read_record[i]);
        if (read_record[i] == "") {
          done = true;
          break;
        }
        bytes_read += read_record[i].length();
        num_lines++;
      }
      progbar.update(bytes_read);
      if (records_processed++ > 100000) break; // do not read the entire data set for just an estimate
      if (done) break;
      num_reads++;
      if (read_record[1].length() < kmer_len) continue;
      num_kmers += read_record[1].length() - kmer_len + 1;
    }
    total_records_processed += records_processed;
    estimated_total_records += records_processed * fileSize / reads_file.zstr_tellg();
    progbar.done();
    barrier();
  }
  //auto all_records_processed = reduce_one(total_records_processed, op_fast_add, 0).wait();
  //SOUT("all records processed ", all_records_processed, "\n");
  double fraction = (double) total_records_processed / (double) estimated_total_records;
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), and found ", num_kmers, " kmers\n");
  auto all_num_lines = reduce_one(num_lines / fraction, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads / fraction, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_all(num_kmers / fraction, op_fast_add).wait();
  int percent = 100.0 * fraction;
  SOUT("Processed ", percent, " % of the estimated total of ", all_num_lines,
       " lines (", all_num_reads, " reads), and found a maximum of ", all_num_kmers, " kmers\n");
  int my_cardinality = all_num_kmers / rank_n();
  SOUT("Cardinality estimated as ", my_cardinality, "\n");
  return my_cardinality;
}

void count_kmers(unsigned kmer_len, int qual_offset, vector<string> reads_fname_list, dist_object<KmerDHT> &kmer_dht, PASS_TYPE pass_type)
{
  Timer timer(__func__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  string read_record[4];
  int max_read_len = 0;
  int64_t num_kmers = 0;
  string progbar_prefix = "";
  switch (pass_type) {
    case BLOOM_SET_PASS: progbar_prefix = "Pass 1: Parsing reads file to setup bloom filter"; break;
    case BLOOM_COUNT_PASS: progbar_prefix = "Pass 2: Parsing reads file to count kmers"; break;
    case NO_BLOOM_PASS: progbar_prefix = "Parsing reads file to count kmers"; break;
  };
  char special = qual_offset + 2;
  IntermittentTimer t_io("reads IO");
  for (auto const &reads_fname : reads_fname_list) {
    string merged_reads_fname = get_merged_reads_fname(reads_fname);
    int64_t bytes_read = 0;
    zstr::ifstream reads_file(merged_reads_fname);
    /*
    stringstream reads_file_buf;
    {
      t_io.start();
      zstr::ifstream reads_file(merged_reads_fname);
      reads_file_buf << reads_file.rdbuf();
      t_io.stop();
    }
    */
    ProgressBar progbar(merged_reads_fname, &reads_file, progbar_prefix);
    while (!reads_file.eof()) {
      bool done = false;
      for (int i = 0; i < 4; i++) {
        getline(reads_file, read_record[i]);
        if (read_record[i] == "") {
          done = true;
          break;
        }
        bytes_read += read_record[i].length();
        num_lines++;
      }
      progbar.update(bytes_read);
      if (done) break;
      string seq = move(read_record[1]);
      string quals = move(read_record[3]);
      num_reads++;
      if (seq.length() < kmer_len) continue;
      if (seq.length() > max_read_len) max_read_len = seq.length();

      // split into kmers
      auto kmers = Kmer::getKmers(seq);

      // disable kmer counting of kmers after a bad quality score (of 2) in the read
      // ... but allow extension counting (if an extention q score still passes the QUAL_CUTOFF)
      size_t foundBadQual = quals.find_first_of(special);
      if (foundBadQual == string::npos) foundBadQual = seq.length();   // remember that the last valid position is length()-1
      int foundBadQualKmer = foundBadQual - kmer_len + 1;
      assert( (int) kmers.size() >= foundBadQualKmer );

      // skip kmers that contain an N
      size_t foundN = seq.find_first_of('N');
      if (foundN == string::npos) foundN = seq.length();

      for (int i = 0; i < kmers.size(); i++) {
        // skip kmers that contain an N
        if (i + kmer_len > foundN) {
          i = foundN; // skip
          // find the next N
          foundN = seq.find_first_of('N', foundN+1);
          if (foundN == string::npos) foundN = seq.length();
          continue;
        }
        char left_base = '0';
        if (i > 0 && quals[i - 1] >= qual_offset + QUAL_CUTOFF) {
          left_base = seq[i - 1];
        }
        char right_base = '0';
        if (i + kmer_len < seq.length() && quals[i + kmer_len] >= qual_offset + QUAL_CUTOFF) {
          right_base = seq[i + kmer_len];
        }
        int count = (i < foundBadQualKmer) ? 1 : 0;
        kmer_dht->add_kmer(kmers[i], left_base, right_base, count, pass_type);
        DBG("kcount add_kmer ", kmers[i].to_string(), " count ", count, "\n");
        num_kmers++;
      }
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates(pass_type);
  }
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads), max read length ", max_read_len, "\n");
  auto all_num_lines = reduce_one(num_lines, op_fast_add, 0).wait();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto all_distinct_kmers = kmer_dht->get_num_kmers();
  auto all_max_read_len = reduce_one(max_read_len, op_fast_max, 0).wait();
  SOUT("Processed a total of ", all_num_lines, " lines (", all_num_reads, " reads), max read length ", all_max_read_len, "\n");
  if (pass_type != BLOOM_SET_PASS) SOUT("Found ", perc_str(all_distinct_kmers, all_num_kmers), " unique kmers\n");
}

void add_ctg_kmers(unsigned kmer_len, string ctgs_fname, dist_object<KmerDHT> &kmer_dht, bool use_bloom, int pass_num_mask)
{
  Timer timer(__func__);
  // first pass over ctgs file - count the number of kmers
  if ((pass_num_mask&1) == 1) {
    zstr::ifstream ctgs_file(ctgs_fname);
    ProgressBar progbar(ctgs_fname, &ctgs_file, "Parsing contigs for extra kmers: pass 1");
    string cname, seq;
    int64_t bytes_read = 0;
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      getline(ctgs_file, seq);
      if (seq == "") break;
      bytes_read += cname.length() + seq.length();
      int64_t num_mers = seq.length() - kmer_len - 1;
      if (num_mers < 0) num_mers = 0;
      if (use_bloom && seq.length() >= kmer_len) {
        auto kmers = Kmer::getKmers(seq);
        if (kmers.size() != seq.length() - kmer_len + 1)
          DIE("kmers size mismatch ", kmers.size(), " != ", (seq.length() - kmer_len + 1), " '", seq, "'");
        for (int i = 1; i < seq.length() - kmer_len; i++) {
          kmer_dht->add_kmer(kmers[i], seq[i - 1], seq[i + kmer_len], 1, CTG_BLOOM_SET_PASS);
        }
      }
      progbar.update(bytes_read);
      progress();
    }
    if (use_bloom) kmer_dht->flush_updates(CTG_BLOOM_SET_PASS);
    progbar.done();
  }
  barrier();
  // second pass over contigs file - insert the new kmers into the dht
  if ((pass_num_mask&2) == 2) {
    string cname, seq, buf;
    int64_t num_kmers = 0;
    int64_t num_ctgs = 0;
    int64_t num_prev_kmers = kmer_dht->get_num_kmers();
    int64_t bytes_read = 0;
    int64_t offset_in_buf = 0;
    zstr::ifstream ctgs_file(ctgs_fname);
    ProgressBar progbar(ctgs_fname, &ctgs_file, "Parsing contigs for extra kmers: pass 2");
    while (!ctgs_file.eof()) {
      getline(ctgs_file, cname);
      if (cname == "") break;
      // seq is folded, so read until we encounter '>' or eof
      seq = "";
      while (true) {
        getline(ctgs_file, buf);
        if (buf == "") break;
        seq += buf;
        auto next_ch = ctgs_file.peek();
        if (next_ch == (int)'>' || next_ch == EOF) break;
      }
      if (seq == "") break;
      bytes_read += cname.length() + seq.length();
      // depth is the last field in the cname
      size_t lastspace = cname.find_last_of(" ");
      if (lastspace == std::string::npos) DIE("Depth is missing from uutigs file ", ctgs_fname, " for contig ", cname);
      double depth = stod(cname.substr(lastspace));
      num_ctgs++;
      if (seq.length() >= kmer_len + 2) {
        auto kmers = Kmer::getKmers(seq);
        if (kmers.size() != seq.length() - kmer_len + 1)
          DIE("kmers size mismatch ", kmers.size(), " != ", (seq.length() - kmer_len + 1), " '", seq, "'");
        for (int i = 1; i < seq.length() - kmer_len; i++) {
          kmer_dht->add_kmer(kmers[i], seq[i - 1], seq[i + kmer_len], depth, CTG_KMERS_PASS);
          num_kmers++;
        }
      }
      progbar.update(bytes_read);
      progress();
    }
    progbar.done();
    kmer_dht->flush_updates(CTG_KMERS_PASS);
    DBG("This rank processed ", num_ctgs, " contigs and ", num_kmers , " kmers\n");
    auto all_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
    auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
    SOUT("Processed a total of ", all_num_ctgs, " contigs and ", all_num_kmers, " kmers\n");
    SOUT("Found ", perc_str(kmer_dht->get_num_kmers() - num_prev_kmers, all_num_kmers), " additional unique kmers\n");
  }
}

