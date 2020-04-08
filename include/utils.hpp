#ifndef __UTILS_H
#define __UTILS_H

#include <sys/stat.h>
#include <unistd.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <regex>
#include <thread>
#include <zlib.h>
#include <upcxx/upcxx.hpp>

#include "bytell_hash_map.hpp"
#include "colors.h"


using std::string;
using std::string_view;
using std::stringstream;
using std::ostringstream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::cerr;
using std::min;

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define __FILEFUNC__ (__FILENAME__ + string(":") + __func__)

#define ONE_MB (1024*1024)
#define ONE_GB (ONE_MB*1024)

#ifdef USE_BYTELL
#define HASH_TABLE ska::bytell_hash_map
#else
#define HASH_TABLE std::unordered_map
#endif

#define CLOCK_NOW std::chrono::high_resolution_clock::now


inline void find_and_replace(std::string& subject, const std::string& search, const std::string& replace) {
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
}

inline std::string_view substr_view(const std::string &s, size_t from, size_t len=string::npos) {
  if (from>=s.size()) return {};
  return std::string_view(s.data() + from, std::min(s.size() - from,len));
}

// this shouldn't really be defined here, but I didn't want yet another header file
enum class QualityLevel {
  SINGLE_PATH_ONLY,
  DEPTH_RESLN_ONLY,
  ALL
};

inline bool file_exists(const string &fname) {
  ifstream ifs(fname, std::ios_base::binary);
  return ifs.good();
}

extern ofstream _logstream;
extern bool _verbose;

inline void init_logger(bool verbose) {
  _verbose = verbose;
  if (!upcxx::rank_me()) {
    bool old_file = file_exists("mhmxx.log");
    _logstream.open("mhmxx.log", std::ofstream::out | std::ofstream::app);
    // the old file should only exist if this is a restart
    if (old_file) _logstream << "\n\n==========  RESTART  ==================\n\n";
  }
}

inline void logger(ostringstream &os) {}

template <typename T, typename... Params>
inline void logger(ostringstream &os, T first, Params... params) {
  os << first;
  logger(os, params ...);
}

template <typename T, typename... Params>
inline void logger(ostream &stream, bool fail, bool serial, bool flush, T first, Params... params) {
  if (serial && upcxx::rank_me()) return;
  ostringstream os;
  os << first;
  logger(os, params ...);
  if (fail) {
    std::cerr << "\n" << KNORM;
    throw std::runtime_error(os.str());
  }
  if (stream.rdbuf() != std::cout.rdbuf() && stream.rdbuf() != std::cerr.rdbuf()) {
    // strip out colors for log file
    string outstr = os.str();
#ifdef USE_COLORS
    for (auto c : COLORS) find_and_replace(outstr, c, "");
#endif
    stream << outstr;
  } else {
    stream << os.str();
  }
  if (flush) stream.flush();
}


#define SOUT(...) do {                                   \
    logger(cout, false, true, true, ##__VA_ARGS__);      \
  } while (0)
#define WARN(...) do {                                                  \
    logger(_logstream, false, false, true, "\n", KRED, "[", upcxx::rank_me(), "] <", __FILENAME__, ":", __LINE__, \
           "> WARNING: ", ##__VA_ARGS__, KNORM, "\n");                  \
    logger(cerr, false, false, true, "\n", KRED, "[", upcxx::rank_me(), "] <", __FILENAME__, ":", __LINE__, \
           "> WARNING: ", ##__VA_ARGS__, KNORM, "\n");                  \
  } while (0)
#define DIE(...) do {                                                   \
    logger(_logstream, false, false, true, "\n", KLRED, "[", upcxx::rank_me(), "] <", __FILENAME__ , ":", __LINE__, \
           "> ERROR: ", ##__VA_ARGS__, KNORM, "\n");                    \
    logger(cerr, true, false, true, "\n", KLRED, "[", upcxx::rank_me(), "] <", __FILENAME__ , ":", __LINE__, \
           "> ERROR: ", ##__VA_ARGS__, KNORM, "\n");                    \
  } while (0)
#define SWARN(...) do {                                                 \
    logger(cerr, false, true, true, "\n", KRED, "WARNING: ", ##__VA_ARGS__, KNORM, "\n\n"); \
    logger(_logstream, false, true, true, "\n", KRED, "WARNING: ", ##__VA_ARGS__, KNORM, "\n\n"); \
  } while (0)
#define SDIE(...) do {                                                  \
    logger(_logstream, false, true, true, "\n", KLRED, "[", upcxx::rank_me(), "] <", __FILENAME__ , ":", __LINE__, \
           "> ERROR: ", ##__VA_ARGS__, KNORM, "\n");                    \
    logger(cerr, true, true, true, "\n", KLRED, "[", upcxx::rank_me(), "] <", __FILENAME__ , ":", __LINE__, \
           "> ERROR: ", ##__VA_ARGS__, KNORM, "\n");                    \
  } while (0)


#define SLOG(...) do {                                              \
    logger(cout, false, true, true, ##__VA_ARGS__);                 \
    logger(_logstream, false, true, true, ##__VA_ARGS__);           \
  } while (0)

#define SLOG_VERBOSE(...) do {                                       \
    if (_verbose) logger(cout, false, true, true, ##__VA_ARGS__);    \
    logger(_logstream, false, true, true, ##__VA_ARGS__);           \
  } while (0)

#ifdef DEBUG
extern ofstream _dbgstream;
#define DBG(...) do {                                                   \
    if (_dbgstream) {                                                   \
      logger(_dbgstream, false, false, true, "<", __FILENAME__, ":", __LINE__, "> ", ##__VA_ARGS__); \
    }                                                                   \
  } while(0)
#else
#define DBG(...)
#endif

static double get_free_mem(void) {
  string buf;
  ifstream f("/proc/meminfo");
  double mem_free = 0;
  while (!f.eof()) {
    getline(f, buf);
    if (buf.find("MemFree") == 0 || buf.find("Buffers") == 0 || buf.find("Cached") == 0) {
      stringstream fields;
      string units;
      string name;
      double mem;
      fields << buf;
      fields >> name >> mem >> units;
      if (units[0] == 'k') mem *= 1024;
      mem_free += mem;
    }
  }
  return mem_free;
}

static string get_size_str(int64_t sz) {
  int64_t abs_sz = std::abs(sz);
  if (abs_sz < 1024) return to_string(abs_sz) + "B";
  double frac = 0;
  string units = "";
  if (abs_sz >= ONE_GB * 1024l) {
    frac = (double)abs_sz / (ONE_GB * 1024l);
    units = "TB";
  } else if (abs_sz >= ONE_GB) {
    frac = (double)abs_sz / ONE_GB;
    units = "GB";
  } else if (abs_sz >= ONE_MB) {
    frac = (double)abs_sz / ONE_MB;
    units = "MB";
  } else if (abs_sz >= 1024) {
    frac = (double)abs_sz / 1024;
    units = "KB";
  }
  ostringstream os;
  if (sz < 0) os << '-';
  os << std::fixed << std::setprecision(2) << frac << units;
  return os.str();
}

class IntermittentTimer {

  std::chrono::time_point<std::chrono::high_resolution_clock> t;
  double t_elapsed, t_interval;
  string name, interval_label;
public:
  IntermittentTimer(const string &name, string interval_label = "") : name{name}, interval_label{interval_label} {
    t_elapsed = 0;
    t_interval = 0;
  }

  void done_all() {
    auto max_t_elapsed = upcxx::reduce_one(t_elapsed, upcxx::op_fast_max, 0).wait();
    auto avg_t_elapsed = upcxx::reduce_one(t_elapsed, upcxx::op_fast_add, 0).wait() / upcxx::rank_n();
    SLOG_VERBOSE(KLCYAN, "--- ", name, " took ", std::setprecision(2), std::fixed, " avg ", avg_t_elapsed,
                 " s max ", max_t_elapsed, " s balance ", (avg_t_elapsed / max_t_elapsed), " ---", KNORM, "\n");
    DBG("--- ", name, " took ", std::setprecision(2), std::fixed, t_elapsed, " s ---\n");
  }

  void done() {
    SLOG_VERBOSE(KLCYAN, "--- ", name, " took ", std::setprecision(2), std::fixed, t_elapsed, " s ---", KNORM, "\n");
    DBG("--- ", name, " took ", std::setprecision(2), std::fixed, t_elapsed, " s ---\n");
  }

  string get_final() {
    ostringstream os;
    os << name << ": " << std::setprecision(2) << std::fixed << t_elapsed;
    return os.str();
  }

  double get_elapsed() {
    return t_elapsed;
  }

  void start() {
    if (!interval_label.empty() && !_verbose) SOUT(KBLUE, std::left, std::setw(40), interval_label + ":", KNORM);
    t = CLOCK_NOW();
  }

  void stop() {
    std::chrono::duration<double> interval = CLOCK_NOW() - t;
    t_interval = interval.count();
    t_elapsed += t_interval;
    if (!interval_label.empty() && !_verbose) SOUT(KBLUE, std::setprecision(2), std::fixed, t_interval, " s", KNORM, "\n");
  }

  double get_interval() {
    return t_interval;
  }
};


class Timer {
  std::chrono::time_point<std::chrono::high_resolution_clock> t;
  string name;
  double init_free_mem;
public:
  Timer(const string &name) {
    t = CLOCK_NOW();
    this->name = name;
    init_free_mem = (!upcxx::rank_me() ? get_free_mem() : 0);
    //if (always_show) SLOG(KLCYAN, "-- ", name, " (", init_free_mem, " GB free) --\n", KNORM);
    //else SLOG_VERBOSE(KLCYAN, "-- ", name, " (", init_free_mem, " GB free) --\n", KNORM);
  }

  ~Timer() {
    std::chrono::duration<double> t_elapsed = CLOCK_NOW() - t;
    DBG(KLCYAN, "-- ", name, " took ", std::setprecision(2), std::fixed, t_elapsed.count(), " s --", KNORM, "\n");
    upcxx::barrier();
    t_elapsed = CLOCK_NOW() - t;
    auto curr_free_mem = (!upcxx::rank_me() ? get_free_mem() : 0);
    auto used_mem = init_free_mem - curr_free_mem;
    SLOG_VERBOSE(KLCYAN, "-- ", name, " took ", std::setprecision(2), std::fixed, t_elapsed.count(), " s ",
                 "(memory used on node 0: ",  get_size_str(used_mem),
                 ", free ", get_size_str(curr_free_mem), ") --", KNORM, "\n");
  }
};


inline string tail(const string &s, int n) {
  return s.substr(s.size() - n);
}

inline string head(const string &s, int n) {
  return s.substr(0, n);
}

inline string perc_str(int64_t num, int64_t tot) {
  ostringstream os;
  os.precision(2);
  os << std::fixed;
  os << num << " (" << 100.0 * num / tot << "%)";
  return os.str();
}

inline std::vector<string> split(const string &s, char delim) {
  std::vector<string> elems;
  std::stringstream ss(s);
  string token;
  while (std::getline(ss, token, delim)) elems.push_back(token);
  return elems;
}

inline void replace_spaces(string &s) {
  for (int i = 0; i < s.size(); i++)
    if (s[i] == ' ') s[i] = '_';
}

inline string revcomp(const string &seq) {
  string seq_rc = "";
  seq_rc.reserve(seq.size());
  for (int i = seq.size() - 1; i >= 0; i--) {
    switch (seq[i]) {
      case 'A': seq_rc += 'T'; break;
      case 'C': seq_rc += 'G'; break;
      case 'G': seq_rc += 'C'; break;
      case 'T': seq_rc += 'A'; break;
      case 'N': seq_rc += 'N'; break;
      default:
        DIE("Illegal char in revcomp of '", seq, "'\n");
    }
  }
  return seq_rc;
}

inline char comp_nucleotide(char ch) {
  switch (ch) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case 'N': return 'N';
      case '0': return '0';
      default: DIE("Illegal char in comp nucleotide of '", ch, "'\n");
  }
  return 0;
}

// returns 1 when it created the directory, 0 otherwise, -1 if there is an error
inline int check_dir(const char *path) {
  if (0 != access(path, F_OK)) {
    if (ENOENT == errno) {
      // does not exist
      // note: we make the directory to be world writable, so others can delete it later if we
      // crash to avoid cluttering up memory
      mode_t oldumask = umask(0000);
      if (0 != mkdir(path, 0777) && 0 != access(path, F_OK)) {
        umask(oldumask);
        DIE("Could not create the directory: ", path, " (", strerror(errno), ")");
        return -1;
      }
      umask(oldumask);
    }
    if (ENOTDIR == errno) {
      // not a directory
      DIE("Expected ", path, " was not a directory");
      return -1;
    }
  } else {
    return 0;
  }
  return 1;
}

// replaces the given path with a rank based path, inserting a rank-based directory
// example:  get_rank_path("path/to/file_output_data.txt", rank) -> "path/to/per_rank/<rankdir>/<rank>/file_output_data.txt"
// of if rank == -1, "path/to/per_rank/file_output_data.txt"
inline bool get_rank_path(string &fname, int rank) {
  const int MAX_RANKS_PER_DIR = 1024;
  char buf[PATH_MAX];
  strcpy(buf, fname.c_str());
  int pathlen = strlen(buf);
  char newPath[PATH_MAX*2+50];
  char *lastslash = strrchr(buf, '/');
  int checkDirs = 0;
  int thisDir;
  char *lastdir = NULL;

  if (pathlen + 25 >= PATH_MAX) {
    WARN("File path is too long (max: ", PATH_MAX, "): ", buf, "\n");
    return false;
  }
  if (lastslash) {
    *lastslash = '\0';
  }
  if (rank < 0) {
    if (lastslash) {
      snprintf(newPath, PATH_MAX*2+50, "%s/per_rank/%s", buf, lastslash + 1);
      checkDirs = 1;
    } else {
      snprintf(newPath, PATH_MAX*2+50, "per_rank/%s", buf);
      checkDirs = 1;
    }
  } else {
    if (lastslash) {
      snprintf(newPath, PATH_MAX*2+50, "%s/per_rank/%08d/%08d/%s", buf, rank / MAX_RANKS_PER_DIR, rank, lastslash + 1);
      checkDirs = 3;
    } else {
      snprintf(newPath, PATH_MAX*2+50, "per_rank/%08d/%08d/%s", rank / MAX_RANKS_PER_DIR, rank, buf);
      checkDirs = 3;
    }
  }
  strcpy(buf, newPath);
  while (checkDirs > 0) {
    strcpy(newPath, buf);
    thisDir = checkDirs;
    while (thisDir--) {
      lastdir = strrchr(newPath, '/');
      if (!lastdir) {
        WARN("What is happening here?!?!\n");
        return false;
      }
      *lastdir = '\0';
    }
    check_dir(newPath);
    checkDirs--;
  }
  fname = buf;
  return true;
}

inline string get_current_time(bool fname_fmt=false) {
  auto t = std::time(nullptr);
  std::ostringstream os;
  if (!fname_fmt) os << std::put_time(localtime(&t), "%D %T");
  else os << std::put_time(localtime(&t), "%y%m%d%H%M%S");
  return os.str();
}

inline int hamming_dist(string_view s1, string_view s2, bool require_equal_len=true) {
  if (require_equal_len && s2.size() != s1.size())//abs((int)(s2.size() - s1.size())) > 1)
    DIE("Hamming distance substring lengths don't match, ", s1.size(), ", ", s2.size(), "\n");
  int d = 0;
  int min_size = min(s1.size(), s2.size());
  for (int i = 0; i < min_size; i++)
    d += (s1[i] != s2[i]);
  return d;
}

static string remove_file_ext(const string &fname) {
  size_t lastdot = fname.find_last_of(".");
  if (lastdot == std::string::npos) return fname;
  return fname.substr(0, lastdot);
}

static string get_basename(const string &fname) {
  size_t i = fname.rfind('/', fname.length());
  if (i != string::npos) return(fname.substr(i + 1, fname.length() - i));
  return fname;
}

static string remove_fname_extension(const string &fname) {
  size_t lastdot = fname.find_last_of(".");
  if (lastdot == std::string::npos) return fname;
  return fname.substr(0, lastdot);
}

static string get_merged_reads_fname(const string &reads_fname) {
  // always relative to the current working directory
  return remove_file_ext(get_basename(reads_fname)) + "-merged.fastq";
}

static int64_t get_file_size(string fname) {
  struct stat s;
  if (stat(fname.c_str(), &s) != 0) return -1;
  return s.st_size;
}

inline void switch_orient(int &start, int &stop, int &len) {
  int tmp = start;
  start = len - stop;
  stop = len - tmp;
}

#define IN_NODE_TEAM() (!(upcxx::rank_me() % upcxx::local_team().rank_n()))

class MemoryTrackerThread {
  std::thread *t = nullptr;
  double start_free_mem, min_free_mem;
  int ticks = 0;
  bool fin = false;
  std::unique_ptr<upcxx::team> node_team;

public:
  void start() {
    // create teams of one process per node
    node_team = std::make_unique<upcxx::team>(upcxx::world().split(IN_NODE_TEAM(), 0));
    if (!IN_NODE_TEAM()) return;
    start_free_mem = get_free_mem();
    auto all_start_mem_free = upcxx::reduce_one(start_free_mem, upcxx::op_fast_add, 0, *node_team).wait();
    SLOG("Initial free memory across all nodes: ", std::setprecision(3), std::fixed, get_size_str(all_start_mem_free), "\n");
    min_free_mem = start_free_mem;
    t = new std::thread([&] {
      while (!fin) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        double free_mem = get_free_mem();
        if (free_mem < min_free_mem) min_free_mem = free_mem;
        ticks++;
      }
    });
  }

  void stop() {
    if (IN_NODE_TEAM()) {
      if (t) {
        fin = true;
        t->join();
        delete t;
        auto peak_mem = start_free_mem - min_free_mem;
        auto all_peak_mem = upcxx::reduce_one(peak_mem, upcxx::op_fast_add, 0, *node_team).wait();
        SLOG("Peak memory used across all nodes: ", get_size_str(all_peak_mem), "\n");
      }
    }
    upcxx::barrier();
    node_team->destroy();
    upcxx::barrier();
  }
};

#endif
