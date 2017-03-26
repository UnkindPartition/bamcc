#include <htslib/sam.h>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

class BamRecord {
  using deleter = void(*)(bam1_t*);
  private:
    unique_ptr<bam1_t, deleter> unq;
  public:
    bool eof;
    BamRecord() : unq(bam_init1(), bam_destroy1) {};
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
  public:
    BamHeader() {};
    BamHeader(bam_hdr_t *header) : sptr(header, bam_hdr_destroy) {};
    bam_hdr_t* ptr() {
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
        throw runtime_error("The file appears corrupt");
      return rec;
    }
};

Graph construct_graph(BamFile file) {
  string last_qname;
  int last_isoform;
  bool first = true;
  Graph g;

  BamRecord rec;
  while (rec = file.next_record(), !rec.eof) {
    string this_qname = rec.qname();
    int this_isoform = rec.ref_id();
    if (!first && last_qname == this_qname)
      add_edge(last_isoform, this_isoform, g);

    last_qname = this_qname;
    last_isoform = this_isoform;
    first = false;
  }

  return g;
}

void write_bam_files(vector<int> component) {
  return;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "USAGE: cc file.bam" << endl;
    exit(1);
  }

  vector<int> component;
  int num_components;

  // TODO check sorting order
  {
    Graph g = construct_graph(BamFile(argv[1], BamMode::Read));
    component.resize(num_vertices(g));
    num_components = connected_components(g, &component[0]);
  } // graph g should be released here

  write_bam_files(component);
}
