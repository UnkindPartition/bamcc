#include <htslib/sam.h>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS> Graph;

Graph construct_graph(samFile *bam, bam_hdr_t *header) {
  string last_qname;
  int last_isoform;
  bool first = true;
  Graph g;

  bam1_t *b = bam_init1();

  while (sam_read1(bam, header, b) >= 0) {
    string this_qname = bam_get_qname(b);
    int this_isoform = b->core.tid;
    if (!first && last_qname == this_qname)
      add_edge(last_isoform, this_isoform, g);

    last_qname = this_qname;
    last_isoform = this_isoform;
    first = false;
  }

  bam_destroy1(b);
  return g;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "USAGE: cc file.bam" << endl;
    exit(1);
  }

  samFile *bam = sam_open(argv[1], "r");
  if (!bam) {
    cerr << "Could not open input file: " << strerror(errno) << endl;
    exit(1);
  }

  bam_hdr_t *header = sam_hdr_read(bam);
  if (!bam) {
    cerr << "Could not read header" << endl;
    exit(1);
  }

  // TODO check sorting order
  Graph g = construct_graph(bam, header);

  sam_close(bam);

  vector<int> component(num_vertices(g));
  int num = connected_components(g, &component[0]);
  // Calculate max comp size
  map<int,int> count;
  for (int i : component) {
    count[i]++;
  }

  for (auto c : count) {
    cout << c.second << endl;
  }
}
