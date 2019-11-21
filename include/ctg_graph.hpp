// cgraph
// Steven Hofmeyr, LBNL, Aug 2018

#ifndef __CTG_GRAPH
#define __CTG_GRAPH

#include <iostream>
#include <map>
#include <fstream>
#include <stdarg.h>
#include <numeric>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "progressbar.hpp"
#include "utils.hpp"
#include "aggr_store.hpp"

using std::pair;
using std::string;
using std::vector;
using std::endl;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::shared_ptr;
using std::make_shared;
using std::tuple;
using std::get;
using std::make_tuple;
using std::min;
using std::max;

//#define DBG_BUILD DBG
#define DBG_BUILD(...)

//#define DBG_WALK DBG
#define DBG_WALK(...)

#define DBG_SPANS DBG
//#define DBG_SPANS(...)

using cid_t = int64_t;


template<typename T>
void serialize_vector(ostringstream &os, vector<T> &vec) 
{
  size_t sz = vec.size();
  os.write((char*)&sz, sizeof(sz));
  for (auto elem : vec) os.write((char*)&elem, sizeof(elem));
}


template<typename T>
void serialize_vector_of_vectors(ostringstream &os, vector<vector<T> > &vec_vec) 
{
  size_t sz = vec_vec.size();
  os.write((char*)&sz, sizeof(sz));
  for (auto vec : vec_vec) serialize_vector(os, vec);
}


template<typename T>
void deserialize_vector(istringstream &is, vector<T> &vec) 
{
  size_t sz;
  is.read((char*)&sz, sizeof(sz));
  vec.clear();
  for (auto i = 0; i < sz; i++) {
    T elem;
    is.read((char*)&elem, sizeof(elem));
    vec.push_back(elem);
  }
}


template<typename T>
void deserialize_vector_of_vectors(istringstream &is, vector<vector<T> > &vec_vec) 
{
  size_t sz;
  is.read((char*)&sz, sizeof(sz));
  vec_vec.clear();
  for (auto i = 0; i < sz; i++) {
    vector<T> vec;
    deserialize_vector(is, vec);
    vec_vec.push_back(vec);
  }
}


inline void serialize_string(ostringstream &os, const string &s) 
{
  size_t l = s.length();
  os.write((char*)&l, sizeof(l));
  if (l) os.write(s.c_str(), l);
}


inline void deserialize_string(istringstream &is, string &s) 
{
  size_t l;
  is.read((char*)&l, sizeof(l));
  if (l) {
    s.resize(l);
    is.read(&s[0], l);
  }
}


enum class Dirn {FORWARD, BACKWARD};
enum class Orient {NORMAL, REVCOMP};

inline string dirn_str(Dirn dirn) 
{
  return (dirn == Dirn::FORWARD ? "forward" : "backward");
}

inline string orient_str(Orient orient)
{
  return (orient == Orient::NORMAL ? "+" : "-");
}


inline Orient flip_orient(Orient orient)
{
  return (orient == Orient::NORMAL ? Orient::REVCOMP : Orient::NORMAL);
}


static const size_t MAX_RNAME_LEN = 80;

struct GapRead {
  char read_name[MAX_RNAME_LEN];
  // used for resolving positive splints
  int gap_start;
  // used for resolving positive spans
  //int rstart, rstop;
  char orient;
  cid_t cid;

  GapRead() {}
  
  //GapRead(const string &read_name, int gap_start, int rstart, int rstop, int orient, cid_t cid) {
  GapRead(const string &read_name, int gap_start, int orient, cid_t cid) {
    if (read_name.size() >= MAX_RNAME_LEN)
      DIE("Read name exceeds buffer size: ", read_name.size(), " > ", MAX_RNAME_LEN, "\n");
    strcpy(this->read_name, read_name.c_str());
    this->gap_start = gap_start;
    //this->rstart = rstart;
    //this->rstop = rstop;
    this->orient = orient;
    this->cid = cid;
  }

  bool operator==(const GapRead &other) const {
    return strcmp(read_name, other.read_name) == 0;
  }
};


struct CidPair {
  cid_t cid1, cid2;

  bool operator==(const CidPair &other) const {
    return (cid1 == other.cid1 && cid2 == other.cid2);
  }

  bool operator!=(const CidPair &other) const {
    return (cid1 != other.cid1 || cid2 != other.cid2);
  }

  friend ostream &operator<<(ostream &os, const CidPair &cids) {
    os << "(" << cids.cid1 << ", " << cids.cid2 << ")";
    return os;
  }
};

static const CidPair NULL_CIDS = CidPair{-1, -1};

namespace std {
  template <>
  struct hash<CidPair> {
    std::size_t operator()(const CidPair& cids) const {
      return (std::hash<cid_t>()(cids.cid1) ^ (std::hash<cid_t>()(cids.cid2) << 1));
    }
  };
}


enum class EdgeType {SPLINT, SPAN};

inline string edge_type_str(EdgeType edge_type) 
{
  switch (edge_type) {
    case EdgeType::SPLINT: return "splint";
    case EdgeType::SPAN: return "span";
    default: return "unknown";
  }
}


struct Edge {
  // the cids of the vertices connected by this edge. Note that the largest number cid is always first
  CidPair cids;
  // the ends correspond to the cids above. 
  int end1, end2;
  int gap;
  int support;
  int aln_len;
  // the best aln score that generated this edge
  int aln_score;
  EdgeType edge_type;
  // the sequence filling a positive gap - if the gap is non-positive, this is empty
  string seq;
  // these flags are set during graph construction to keep track of errors in edge construction
  bool mismatch_error, conflict_error, excess_error, short_aln;
  // contains information of reads that map to a positive gap - used for filling the gap
  vector<GapRead> gap_reads;

  string serialize() {
    ostringstream os(stringstream::binary);
    os.write((char*)&cids, sizeof(cids));
    os.write((char*)&end1, sizeof(end1));
    os.write((char*)&end2, sizeof(end2));
    os.write((char*)&gap, sizeof(gap));
    os.write((char*)&support, sizeof(support));
    os.write((char*)&aln_len, sizeof(aln_len));
    os.write((char*)&aln_score, sizeof(aln_score));
    os.write((char*)&edge_type, sizeof(edge_type));
    serialize_string(os, seq);
    os.write((char*)&mismatch_error, sizeof(mismatch_error));
    os.write((char*)&conflict_error, sizeof(conflict_error));
    os.write((char*)&excess_error, sizeof(excess_error));
    os.write((char*)&short_aln, sizeof(short_aln));
    serialize_vector(os, gap_reads);
    return os.str();
  }

  void deserialize(const string &s) {
    istringstream is(s, stringstream::binary);
    is.read((char*)&cids, sizeof(cids));
    is.read((char*)&end1, sizeof(end1));
    is.read((char*)&end2, sizeof(end2));
    is.read((char*)&gap, sizeof(gap));
    is.read((char*)&support, sizeof(support));
    is.read((char*)&aln_len, sizeof(aln_len));
    is.read((char*)&aln_score, sizeof(aln_score));
    is.read((char*)&edge_type, sizeof(edge_type));
    deserialize_string(is, seq);
    is.read((char*)&mismatch_error, sizeof(mismatch_error));
    is.read((char*)&conflict_error, sizeof(conflict_error));
    is.read((char*)&excess_error, sizeof(excess_error));
    is.read((char*)&short_aln, sizeof(short_aln));
    deserialize_vector(is, gap_reads);
  }
};


struct Vertex {
  cid_t cid;
  int clen;
  double depth;
  // set to true if visited in previous round
  bool visited;
  upcxx::global_ptr<char> seq_gptr;
  // the neighbors at the different ends
  vector<cid_t> end5;
  vector<cid_t> end3;
  // the merged series of nbs
  // FIXME: when using spans make sure these are valid
  vector<vector<cid_t> > end5_merged;
  vector<vector<cid_t> > end3_merged;
  // book-keeping fields for resolving walk conflicts between ranks - 
  // choose the walk with the longest scaffold, and if there is a tie, choose the highest rank
  int walk_score;
  int walk_rank;
  int walk_i;

  string serialize() {
    ostringstream os(stringstream::binary);
    os.write((char*)&cid, sizeof(cid));
    os.write((char*)&clen, sizeof(clen));
    os.write((char*)&depth, sizeof(depth));
    os.write((char*)&visited, sizeof(visited));
    os.write((char*)&seq_gptr, sizeof(seq_gptr));
    serialize_vector(os, end5);
    serialize_vector(os, end3);
    serialize_vector_of_vectors(os, end5_merged);
    serialize_vector_of_vectors(os, end3_merged);
    os.write((char*)&walk_score, sizeof(walk_score));
    os.write((char*)&walk_rank, sizeof(walk_rank));
    os.write((char*)&walk_i, sizeof(walk_i));
    return os.str();
  }

  void deserialize(const string &s) {
    istringstream is(s, stringstream::binary);
    is.read((char*)&cid, sizeof(cid));
    is.read((char*)&clen, sizeof(clen));
    is.read((char*)&depth, sizeof(depth));
    is.read((char*)&visited, sizeof(visited));
    is.read((char*)&seq_gptr, sizeof(seq_gptr));
    deserialize_vector(is, end5);
    deserialize_vector(is, end3);
    deserialize_vector_of_vectors(is, end5_merged);
    deserialize_vector_of_vectors(is, end3_merged);
    is.read((char*)&walk_score, sizeof(walk_score));
    is.read((char*)&walk_rank, sizeof(walk_rank));
    is.read((char*)&walk_i, sizeof(walk_i));
  }
};


class CtgGraph {
private:
  using vertex_map_t = HASH_TABLE<cid_t, Vertex>;
  using edge_map_t = HASH_TABLE<CidPair, Edge>;
  using reads_map_t = HASH_TABLE<string, string>;
  upcxx::dist_object<vertex_map_t> vertices;
  upcxx::dist_object<edge_map_t> edges;
  upcxx::dist_object<reads_map_t> read_seqs;
  HASH_TABLE<cid_t, shared_ptr<Vertex> > vertex_cache;
  HASH_TABLE<CidPair, shared_ptr<Edge> > edge_cache;
  const int MAX_CACHE_SIZE = 2000000;

  struct VertexDepthInfo {
    cid_t cid;
    double depth;
    int clen;
  };
    
  AggrStore<VertexDepthInfo> vertex_depth_store;

  struct EdgeGapReadInfo {
    CidPair cids;
    GapRead gap_read;
    int aln_len;
    int aln_score;
  };

  AggrStore<EdgeGapReadInfo> edge_gap_read_store;
  
  edge_map_t::iterator edge_iter;
  vertex_map_t::iterator vertex_iter;

  size_t get_vertex_target_rank(cid_t cid) {
    return std::hash<cid_t>{}(cid) % upcxx::rank_n();
  }

  size_t get_edge_target_rank(CidPair &cids) {
    return std::hash<CidPair>{}(cids) % upcxx::rank_n();
  }

  size_t get_read_target_rank(const string &r) {
    return std::hash<string>{}(r) % upcxx::rank_n();
  }

  struct UpdateDepthFunc {
    void operator()(VertexDepthInfo &vertex_depth_info, upcxx::dist_object<vertex_map_t> &vertices) {
      const auto it = vertices->find(vertex_depth_info.cid);
      if (it == vertices->end()) DIE("could not fetch vertex ", vertex_depth_info.cid, "\n");
      auto v = &it->second;
      v->depth += (vertex_depth_info.depth / vertex_depth_info.clen);
      if (v->clen != vertex_depth_info.clen)
        DIE("Mismatched clen for ctg ", vertex_depth_info.cid, ": ", vertex_depth_info.clen, " != ", v->clen, "\n");
    }
  };
  dist_object<UpdateDepthFunc> update_depth_func;


  struct AddEdgeGapReadFunc {
    void operator()(EdgeGapReadInfo &edge_gap_read_info, upcxx::dist_object<edge_map_t> &edges) {
      const auto it = edges->find(edge_gap_read_info.cids);
      if (it == edges->end()) DIE("SPAN edge not found ", edge_gap_read_info.cids.cid1, " ", edge_gap_read_info.cids.cid2, "\n");
      auto edge = &it->second;
      edge->gap_reads.push_back(edge_gap_read_info.gap_read);
      edge->aln_len = max(edge->aln_len, edge_gap_read_info.aln_len);
      edge->aln_score = max(edge->aln_score, edge_gap_read_info.aln_score);
    }
  };
  dist_object<AddEdgeGapReadFunc> add_edge_gap_read_func;
  
public:
  int max_read_len;

  CtgGraph() : vertices({}), edges({}), read_seqs({}), vertex_cache({}), edge_cache({}), vertex_depth_store({}),
               update_depth_func({}), edge_gap_read_store({}), add_edge_gap_read_func({}) {
    vertex_cache.reserve(MAX_CACHE_SIZE);
    edge_cache.reserve(MAX_CACHE_SIZE);
    vertex_depth_store.set_size("vertex depths store", ONE_MB);
    edge_gap_read_store.set_size("edge gaps store", ONE_MB);
  }

  void clear() {
    for (auto it = vertices->begin(); it != vertices->end(); ) {
      it = vertices->erase(it);
    }
    for (auto it = edges->begin(); it != edges->end(); ) {
      it = edges->erase(it);
    }
  }
  
  ~CtgGraph() {
    clear();
  }
  
  int64_t get_num_vertices(bool all = false) {
    if (!all) return upcxx::reduce_one(vertices->size(), upcxx::op_fast_add, 0).wait();
    else return upcxx::reduce_all(vertices->size(), upcxx::op_fast_add).wait();
  }
  
  int64_t get_local_num_vertices(void) {
    return vertices->size();
  }
  
  int64_t get_num_edges(bool all = false) {
    if (!all) return upcxx::reduce_one(edges->size(), upcxx::op_fast_add, 0).wait();
    else return upcxx::reduce_all(edges->size(), upcxx::op_fast_add).wait();
  }

  int64_t get_local_num_edges(void) {
    return edges->size();
  }

  shared_ptr<Vertex> get_vertex(cid_t cid) {
    size_t target_rank = get_vertex_target_rank(cid);
    if (target_rank == upcxx::rank_me()) {
      upcxx::progress();
      const auto it = vertices->find(cid);
      if (it == vertices->end()) return nullptr;
      return make_shared<Vertex>(it->second);
    }
    return upcxx::rpc(target_rank, 
                      [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid) {
                        const auto it = vertices->find(cid);
                        if (it == vertices->end()) return string("");
                        return it->second.serialize();
                      }, vertices, cid).then(
                        [](string s) -> shared_ptr<Vertex> {
                          if (s == "") return nullptr;
                          auto v = make_shared<Vertex>();
                          v->deserialize(s);
                          return v;
                        }).wait();
  }

  int get_vertex_clen(cid_t cid) {
    return upcxx::rpc(get_vertex_target_rank(cid),
                      [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid) {
                        const auto it = vertices->find(cid);
                        if (it == vertices->end()) return -1;
                        return it->second.clen;
                      }, vertices, cid).wait();
  }

  shared_ptr<Vertex> get_local_vertex(cid_t cid) {
    const auto it = vertices->find(cid);
    if (it == vertices->end()) return nullptr;
    return make_shared<Vertex>(it->second);
  }
  
  void set_vertex_visited(cid_t cid) {
    upcxx::rpc(get_vertex_target_rank(cid),
               [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid) {
                 const auto it = vertices->find(cid);
                 if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
                 auto v = &it->second;
                 v->visited = true;
               }, vertices, cid).wait();
  }

  void update_vertex_walk(cid_t cid, int walk_score, int walk_i) {
    upcxx::rpc(get_vertex_target_rank(cid),
               [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid, int walk_score, int walk_i, int myrank) {
                 const auto it = vertices->find(cid);
                 if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
                 auto v = &it->second;
                 if (myrank == v->walk_rank) {
                   // same rank, select in favor of highest score
                   if (walk_score > v->walk_score) {
                     v->walk_score = walk_score;
                     v->walk_i = walk_i;
                   }
                 } else {
                   // different rank, select highest score, break ties with highest rank
                   if ((walk_score == v->walk_score && myrank > v->walk_rank) || walk_score > v->walk_score) {
                     v->walk_score = walk_score;
                     v->walk_rank = myrank;
                     v->walk_i = walk_i;
                   }
                 }
               }, vertices, cid, walk_score, walk_i, upcxx::rank_me()).wait();
  }
  
  Vertex *get_first_local_vertex() {
    vertex_iter = vertices->begin();
    if (vertex_iter == vertices->end()) return nullptr;
    auto v = &vertex_iter->second;
    vertex_iter++;
    return v;
  }
  
  Vertex *get_next_local_vertex() {
    if (vertex_iter == vertices->end()) return nullptr;
    auto v = &vertex_iter->second;
    vertex_iter++;
    return v;
  }

  void add_vertex(Vertex &v, const string &seq) {
    v.clen = seq.length();
    v.seq_gptr = upcxx::allocate<char>(v.clen + 1);
    strcpy(v.seq_gptr.local(), seq.c_str());
    upcxx::rpc(get_vertex_target_rank(v.cid),
               [](upcxx::dist_object<vertex_map_t> &vertices, string s) {
                 Vertex v;
                 v.deserialize(s);
                 v.visited = false;
                 vertices->insert({v.cid, v});
               }, vertices, v.serialize()).wait();
  }

  void add_vertex_nb(cid_t cid, cid_t nb, char end) {
    upcxx::rpc(get_vertex_target_rank(cid),
               [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid, cid_t nb, int end) {
                 const auto it = vertices->find(cid);
                 if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
                 auto v = &it->second;
                 if (end == 5) v->end5.push_back(nb);
                 else v->end3.push_back(nb);
               }, vertices, cid, nb, end).wait();
  }

  string get_vertex_seq(upcxx::global_ptr<char> seq_gptr, int64_t seq_len) {
    char buf[seq_len + 1];
    upcxx::rget(seq_gptr, buf, seq_len + 1).wait();
    string s(buf);
    return s;
  }

  void mark_edge_short_aln(CidPair cids) {
    upcxx::rpc(get_edge_target_rank(cids),
               [](upcxx::dist_object<edge_map_t> &edges, CidPair cids) {
                 const auto it = edges->find(cids);
                 if (it == edges->end()) DIE("Can't find edge ", cids);
                 it->second.short_aln = true;
               }, edges, cids);
  }

  shared_ptr<Edge> get_edge(cid_t cid1, cid_t cid2) {
    CidPair cids = { .cid1 = cid1, .cid2 = cid2 };
    if (cid1 < cid2) std::swap(cids.cid1, cids.cid2);
    size_t target_rank = get_edge_target_rank(cids);
    if (target_rank == upcxx::rank_me()) {
      upcxx::progress();
      const auto it = edges->find(cids);
      if (it == edges->end()) return nullptr;
      return make_shared<Edge>(it->second);
    }
    return upcxx::rpc(target_rank,
                      [](upcxx::dist_object<edge_map_t> &edges, CidPair cids) -> string {
                        const auto it = edges->find(cids);
                        if (it == edges->end()) return string("");
                        return it->second.serialize();
                      }, edges, cids).then(
                        [](string s) -> shared_ptr<Edge> {
                          if (s == "") return nullptr;
                          auto edge = make_shared<Edge>();
                          edge->deserialize(s);
                          return edge;
                        }).wait();
  }

  Edge *get_first_local_edge() {
    edge_iter = edges->begin();
    if (edge_iter == edges->end()) return nullptr;
    auto edge = &edge_iter->second;
    edge_iter++;
    return edge;
  }
  
  Edge *get_next_local_edge() {
    if (edge_iter == edges->end()) return nullptr;
    auto edge = &edge_iter->second;
    edge_iter++;
    return edge;
  }

  void add_or_update_edge(Edge &edge) {
    upcxx::rpc(get_edge_target_rank(edge.cids),
               [](upcxx::dist_object<edge_map_t> &edges, string s) {
                 Edge new_edge;
                 new_edge.deserialize(s);
                 const auto it = edges->find(new_edge.cids);
                 if (it == edges->end()) {
                   // not found, always insert
                   edges->insert({new_edge.cids, new_edge});
                 } else {
                   auto edge = &it->second;
                   // always a failure
                   if (edge->mismatch_error) return;
                   if (edge->edge_type == EdgeType::SPLINT && new_edge.edge_type == EdgeType::SPAN) {
                     DBG_BUILD("span confirms splint: ", edge->cids, "\n");
                     edge->support++;
                   } else {
                     if (edge->edge_type == EdgeType::SPLINT && new_edge.edge_type == EdgeType::SPLINT) {
                       // check for mismatches in gap size if they're both splints
                       if (abs(new_edge.gap - edge->gap) > 2) {
                         DBG_BUILD("gap mismatch for ", new_edge.cids, " ", new_edge.gap, " != ", edge->gap, "\n");
                         edge->mismatch_error = true;
                         // this edge will be dropped
                         return;
                       }
                       edge->gap = min(new_edge.gap, edge->gap);
                     }
                     if (edge->edge_type == EdgeType::SPAN && new_edge.edge_type == EdgeType::SPAN)
                       edge->gap += new_edge.gap;
                     edge->support++;
                     edge->aln_len = max(edge->aln_len, new_edge.aln_len);
                     edge->aln_score = max(edge->aln_score, new_edge.aln_score);
                     // count conflicts
                     if (edge->end1 != new_edge.end1 || edge->end2 != new_edge.end2) edge->conflict_error = true;
                     if (new_edge.gap > 0) {
                       // add reads to positive gap for splints
                       //if (new_edge.gap_reads.size() == 0) DIE("trying to get a read from an empty array\n");
                       for (auto gap_read : new_edge.gap_reads) 
                         edge->gap_reads.push_back(gap_read);
                     }
                   }
                 }
               }, edges, edge.serialize()).wait();
  }

  void add_edge_gap_read(CidPair cids, GapRead gap_read, int aln_len, int aln_score) {
    auto target_rank = get_edge_target_rank(cids);
    EdgeGapReadInfo edge_gap_read_info = { cids, gap_read, aln_len, aln_score };
    edge_gap_read_store.update(target_rank, edge_gap_read_info, add_edge_gap_read_func, edges);
  }

  void flush_edge_gap_reads() {
    edge_gap_read_store.flush_updates(add_edge_gap_read_func, edges);
  }

  void purge_error_edges(int64_t *mismatched, int64_t *conflicts, int64_t *empty_spans) {
    for (auto it = edges->begin(); it != edges->end(); ) {
      auto edge = make_shared<Edge>(it->second);
      if (edge->mismatch_error) {
        (*mismatched)++;
        it = edges->erase(it);
      } else if (edge->edge_type == EdgeType::SPAN && edge->gap > 0 && !edge->gap_reads.size()) {
        // don't use positive span gaps without filler
        (*empty_spans)++;
        it = edges->erase(it);
      } else {
        // retain the conflicts to prevent falsely choosing a path when it should be a fork
        if (edge->conflict_error) (*conflicts)++;
        it++;
      }
    }
  }
  
  int64_t purge_excess_edges() {
    int64_t excess = 0;
    for (auto it = edges->begin(); it != edges->end(); ) {
      auto edge = make_shared<Edge>(it->second);
      if (edge->excess_error) {
        excess++;
        it = edges->erase(it);
      } else {
        it++;
      }
    }
    return excess;
  }

  void remove_nb(cid_t cid, int end, cid_t nb) {
    upcxx::rpc(get_vertex_target_rank(cid),
               [](upcxx::dist_object<vertex_map_t> &vertices, cid_t cid, int end, cid_t nb) {
                 const auto it = vertices->find(cid);
                 if (it == vertices->end()) DIE("could not fetch vertex ", cid, "\n");
                 auto v = &it->second;
                 if (end == 3) {
                   for (auto it = v->end3.begin(); it != v->end3.end(); it++) {
                     if (*it == nb) {
                       v->end3.erase(it);
                       return;
                     }
                   }
                 } else {
                   for (auto it = v->end5.begin(); it != v->end5.end(); it++) {
                     if (*it == nb) {
                       v->end5.erase(it);
                       return;
                     }
                   }
                 }
                 DIE("Could not find the nb to remove");
               }, vertices, cid, end, nb).wait();
  }
               
  int64_t purge_short_aln_edges() {
    int64_t num_short = 0;
    for (auto it = edges->begin(); it != edges->end(); ) {
      auto edge = make_shared<Edge>(it->second);
      if (edge->short_aln) {
        num_short++;
        remove_nb(edge->cids.cid1, edge->end1, edge->cids.cid2);
        remove_nb(edge->cids.cid2, edge->end2, edge->cids.cid1);
        it = edges->erase(it);
      } else {
        it++;
      }
    }
    return num_short;
  }
  
  void add_pos_gap_read(const string &read_name) {
    upcxx::rpc(get_read_target_rank(read_name),
               [](upcxx::dist_object<reads_map_t> &read_seqs, string read_name) {
                 read_seqs->insert({read_name, ""});
               }, read_seqs, read_name).wait();
  }

  bool update_read_seq(const string &read_name, const string &seq) {
    return upcxx::rpc(get_read_target_rank(read_name),
                      [](upcxx::dist_object<reads_map_t> &read_seqs, string read_name, string seq) {
                        auto it = read_seqs->find(read_name);
                        if (it == read_seqs->end()) return false;
                        (*read_seqs)[read_name] = seq;
                        return true;
                      }, read_seqs, read_name, seq).wait();
  }
  
  string get_read_seq(const string &read_name) {
    return upcxx::rpc(get_read_target_rank(read_name),
                      [](upcxx::dist_object<reads_map_t> &read_seqs, string read_name) -> string {
                        const auto it = read_seqs->find(read_name);
                        if (it == read_seqs->end()) return string("");
                        return it->second;
                      }, read_seqs, read_name).wait();
  }
  
  size_t get_num_read_seqs(bool all = false) {
    if (!all) return upcxx::reduce_one(read_seqs->size(), upcxx::op_fast_add, 0).wait();
    else return upcxx::reduce_all(read_seqs->size(), upcxx::op_fast_add).wait();
  }
  
  int get_other_end(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge = nullptr) {
    if (!edge) edge = get_edge_cached(v1->cid, v2->cid);
    // the cids in the edge are organized as (largest, smallest), so the edges are determined that way too
    return v1->cid >= v2->cid ? edge->end2 : edge->end1;
  }

  int get_other_end_local(Vertex* v1, shared_ptr<Vertex> v2, shared_ptr<Edge> edge  = nullptr) {
    if (!edge) edge = get_edge(v1->cid, v2->cid);
    // the cids in the edge are organized as (largest, smallest), so the edges are determined that way too
    return v1->cid >= v2->cid ? edge->end2 : edge->end1;
  }
  
  void clear_caches() {
    vertex_cache.clear();
    edge_cache.clear();
  }

  shared_ptr<Vertex> get_vertex_cached(cid_t cid) {
    auto it = vertex_cache.find(cid);
    if (it != vertex_cache.end()) {
      upcxx::progress();
      return it->second;
    }
    auto v = get_vertex(cid);
    // load factor around .5
    if (vertex_cache.size() < MAX_CACHE_SIZE / 2) vertex_cache[cid] = v;
    //else DBG("vertex cache is full\n");
    return v;
  }

  shared_ptr<Edge> get_edge_cached(cid_t cid1, cid_t cid2) {
    CidPair cids = { .cid1 = cid1, .cid2 = cid2 };
    if (cid1 < cid2) std::swap(cids.cid1, cids.cid2);
    auto it = edge_cache.find(cids);
    if (it != edge_cache.end()) {
      upcxx::progress();
      return it->second;
    }
    auto edge = get_edge(cids.cid1, cids.cid2);
    if (edge_cache.size() < MAX_CACHE_SIZE / 2) edge_cache[cids] = edge;
    //else DBG("edge cache is full\n");
    return edge;
  }

  void print_stats(string graph_fname="") {
    Timer timer(__func__);
    auto get_avg_min_max = [](vector<int64_t> &vals) {
      int64_t total = std::accumulate(vals.begin(), vals.end(), 0);
      int64_t max_val = *std::max_element(vals.begin(), vals.end());
      int64_t min_val = *std::min_element(vals.begin(), vals.end());
      int64_t all_min_val =  upcxx::reduce_one(min_val, upcxx::op_fast_min, 0).wait();
      int64_t all_max_val =  upcxx::reduce_one(max_val, upcxx::op_fast_max, 0).wait();
      double all_total = upcxx::reduce_one(total, upcxx::op_fast_add, 0).wait();
      size_t all_nvals = upcxx::reduce_one(vals.size(), upcxx::op_fast_add, 0).wait();
      ostringstream os;
      os.precision(2);
      os << std::fixed;
      os << (all_total / all_nvals) << " [" << all_min_val << ", " << all_max_val << "]";
      return os.str();
    };

    ofstream graph_stream;
    if (!graph_fname.empty()) {
      get_rank_path(graph_fname, rank_me());
      graph_stream.open(graph_fname);
    }
    vector<int64_t> depths;
    depths.reserve(get_local_num_vertices());
    vector<int64_t> clens;
    clens.reserve(get_local_num_vertices());
    for (auto v = get_first_local_vertex(); v != nullptr; v = get_next_local_vertex()) {
      if (v->clen >= ASSM_CLEN_THRES) {
        depths.push_back(round(v->depth));
        clens.push_back(v->clen);
      }
      if (!graph_fname.empty()) graph_stream << v->cid << "," << v->clen << "," << v->depth << endl;
    }
    vector<int64_t> supports;
    supports.reserve(get_local_num_edges());
    vector<int64_t> aln_lens;
    aln_lens.reserve(get_local_num_edges());
    vector<int64_t> aln_scores;
    aln_scores.reserve(get_local_num_edges());
    vector<int64_t> gaps;
    gaps.reserve(get_local_num_edges());
    {
      ProgressBar progbar(get_local_num_edges(), "Compute graph stats");
      for (auto edge = get_first_local_edge(); edge != nullptr; edge = get_next_local_edge()) {
        aln_lens.push_back(edge->aln_len);
        aln_scores.push_back(edge->aln_score);
        auto clen1 = get_vertex_clen(edge->cids.cid1);
        auto clen2 = get_vertex_clen(edge->cids.cid2);
        if (clen1 >= ASSM_CLEN_THRES || clen2 >= ASSM_CLEN_THRES) {
          supports.push_back(edge->support);
          gaps.push_back(edge->gap);
        }
        if (!graph_fname.empty()) {
          graph_stream << edge->cids.cid1 << "," << edge->cids.cid2 << "," << edge->end1 << "," << edge->end2 << "," << edge->gap
                       << "," << edge->support << "," << edge->aln_len << "," << edge->aln_score << ","
                       << edge_type_str(edge->edge_type) << endl;
        }
        progbar.update();
      }
      progbar.done();
    }

    auto num_vertices = get_num_vertices();
    auto num_edges = get_num_edges();
    SLOG("Graph statistics:\n");
    SLOG("    vertices:  ", num_vertices, "\n");
    SLOG("    edges:     ", num_edges, "\n");
    SLOG("    degree:    ", (double)num_edges / num_vertices, "\n");
    SLOG("    aln_len:   ", get_avg_min_max(aln_lens), "\n");
    SLOG("    aln_score: ", get_avg_min_max(aln_scores), "\n");
    SLOG("  for contigs >= ", ASSM_CLEN_THRES, " length:\n");
    SLOG("    depth:     ", get_avg_min_max(depths), "\n");
    SLOG("    clen:      ", get_avg_min_max(clens), "\n");
    SLOG("    support:   ", get_avg_min_max(supports), "\n");
    SLOG("    gap:       ", get_avg_min_max(gaps), "\n");
  }

};

#endif
