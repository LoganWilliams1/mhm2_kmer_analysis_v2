#include <chrono>
#include <iomanip>
#include <iostream>
#include <upcxx/upcxx.hpp>

using namespace std;
using namespace upcxx;

using timepoint_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
using std::chrono::seconds;
using duration_seconds = std::chrono::duration<double>;

timepoint_t now() { return std::chrono::high_resolution_clock::now(); }

#define ASYNC

int main(int argc, char **argv) {
  upcxx::init();
  if (argc != 2) {
    if (!rank_me()) cout << "ERROR: need to pass buffer size in bytes\n";
    upcxx::finalize();
    return 1;
  }
  size_t buf_size = atoi(argv[1]);
  barrier();
  if (!rank_me()) cerr << "Running test of rgets with " << buf_size << " buffer size\n";
  barrier();
  dist_object<global_ptr<char>> buf(new_array<char>(buf_size));
  vector<global_ptr<char>> other_bufs(rank_n());
  barrier();
  for (int i = 0; i < rank_n(); i++) other_bufs[i] = buf.fetch(i).wait();
  barrier();
  char *receive_buf = new char[buf_size];
  size_t num_rounds = 10000000;
  vector<size_t> get_counts(rank_n(), 0);
  auto t = now();
  size_t num_gets = 0;
  auto f = make_future<>();
  for (size_t round = 0; round < num_rounds; round++) {
    int target = rand() % rank_n();
    get_counts[target]++;
    auto frget = rget(other_bufs[target], receive_buf, buf_size);
#ifdef ASYNC
    f = when_all(f, frget);
#else
    fget.wait();
#endif
    num_gets++;
  }
  f.wait();
  barrier();
  
  auto t_elapsed = ((duration_seconds)(now() - t)).count();
  double bandwidth = (double)(num_gets * buf_size) / 1024.0 / 1024.0 / t_elapsed;
  
 if (!rank_me()) {
    cout << "Test took " << t_elapsed << " seconds\n";
    cout << "Average rget time is " << fixed << setprecision(3)
         << (1000000.0 * (double)t_elapsed / num_gets) << " us\n";
    cout << "Bandwidth per rank " << fixed << setprecision(4) << bandwidth << " MB/s\n";
    double tot = 0;
    size_t max_count = 0;
    for (int i = 0; i < rank_n(); i++) {
      tot+= get_counts[i];
      max_count = max(max_count, get_counts[i]);
    }
    double avg = tot / rank_n();
    double balance = avg / max_count;
    cout << "Rank 0 received on average " << fixed << setprecision(2) << avg
         << " rgets per rank, max " << max_count << " (balance " << balance << ")" << endl;
  }
  barrier();
  upcxx::finalize();
  return 0;
}
