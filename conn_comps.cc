#include <htslib/sam.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

class BamRecord {
  using deleter = void(*)(bam1_t*);
  private:
    unique_ptr<bam1_t, deleter> unq;
    static unique_ptr<bam1_t, deleter> alloc() {
      if (auto p = bam_init1()) {
        return unique_ptr<bam1_t, deleter>(p, bam_destroy1);
      } else {
        throw bad_alloc();
      }
    }
  public:
    bool eof{false};
    BamRecord() : unq(alloc()) {};
    bam1_t* ptr() {
      return unq.get();
    }
    string qname() {
      return string(bam_get_qname(unq.get()));
    }
    int32_t ref_id() {
      return unq.get()->core.tid;
    }
};

class BamHeader {
  using deleter = void(*)(bam_hdr_t*);
  private:
    shared_ptr<bam_hdr_t> sptr;
    static shared_ptr<bam_hdr_t> alloc() {
      if (auto p = bam_hdr_init()) {
        return shared_ptr<bam_hdr_t>(p, bam_hdr_destroy);
      } else {
        throw bad_alloc();
      }
    }
  public:
    BamHeader() {};
    BamHeader(bam_hdr_t *header) : sptr(header, bam_hdr_destroy) {};
    bam_hdr_t* ptr() const {
      return sptr.get();
    }
    bool const operator!(){
      return !sptr;
    }
};

enum class BamMode { Read, Write };

class BamFile {
  using deleter = int(*)(htsFile*);
  private:
    unique_ptr<samFile, deleter> unq{nullptr, hts_close};
    BamHeader header;
    void read_header() {
      if (!header) {
        bam_hdr_t *hdr = sam_hdr_read(unq.get());
        if (!hdr) {
          throw runtime_error("Could not read header");
        }
        header = BamHeader(hdr);
      }
    }
  public:
    BamFile(const char *path, BamMode mode) {
      const char *mode_str;
      switch (mode) {
        case BamMode::Read: mode_str = "r"; break;
        case BamMode::Write: mode_str = "wb"; break;
      }
      samFile *file = sam_open(path, mode_str);
      if (!file) {
        stringstream msg;
        msg << "Could not open input file: " << strerror(errno);
        throw runtime_error(msg.str());
      }
      unq = unique_ptr<samFile, deleter>(file, hts_close);
    };
    samFile* ptr() {
      return unq.get();
    }
    BamHeader get_header() {
      read_header();
      return header;
    }
    BamRecord next_record() {
      BamRecord rec;

      // ensure that the header is read
      read_header();

      int r = sam_read1(unq.get(), header.ptr(), rec.ptr());
      if (r == -1)
        rec.eof = true;
      if (r < -1)
        throw runtime_error("The file appears corrupted");
      return rec;
    }
};

class ConnectedComponents {
  public:
    // component_of[i] is the component number of the element i
    vector<int> component_of;

    ConnectedComponents(const Graph &g) : component_of(num_vertices(g)) {
      connected_components(g, &component_of[0]);
    }
    void write(const char *path, const BamHeader &bh) {
      ofstream out(path);
      out << "seqid\tseqname\tcomponent\n";
      bam_hdr_t *hdr = bh.ptr();
      for (unsigned i = 0; i < component_of.size(); i++) {
        out << i << "\t" << hdr->target_name[i] << "\t" << component_of[i] << "\n";
      }
    }
};

Graph construct_graph(BamFile &file) {
  unordered_map<string,vector<int>> map;
  BamRecord rec;

  while (rec = file.next_record(), !rec.eof) {
    map[rec.qname()].push_back(rec.ref_id());
  }

  Graph g;
  for (auto pair : map) {
    auto vec = pair.second;
    for (size_t i = 1; i < vec.size(); i++) {
      add_edge(vec[0], vec[i], g);
    }
  }

  return g;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "USAGE: cc input.bam output.tsv" << endl;
    exit(1);
  }

  try {
    BamFile bam(argv[1], BamMode::Read);
    Graph g = construct_graph(bam);
    ConnectedComponents comps(g);
    comps.write(argv[2], bam.get_header());
  } catch (const std::exception &exc)
  {
    std::cerr << exc.what() << endl;
    exit(1);
  }
}
