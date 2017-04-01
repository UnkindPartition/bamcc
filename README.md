# bamcc

## Description

A set of sequence alignments (represented as a BAM/SAM/CRAM file) defines
a bipartite graph in which the nodes are reference and template sequences
and edges are alignments between a template and a reference.

bamcc splits this graph into connected components and reports the mapping from
reference sequences to their components.

This is useful e.g. when analyzing multireads in an RNA-Seq experiment,
where connected components can be analyzed independently.

## Usage

```
bamcc input.bam output.tsv
```

The `output.tsv` file will look like this:

```
seqid  seqname      component
0      FBtr0005088  0
1      FBtr0006151  1
2      FBtr0070000  2
3      FBtr0070002  0
4      FBtr0070003  2
5      FBtr0070006  3
```

where:

1. `seqid` is the 0-based number of the reference sequence.
    Sequence numbers are defined by the input file and are stable.
2.  `seqname` is the reference sequence name.
3.  `component` is the 0-based component number to which the reference has been
    assigned.

To extract, say, the 17th component into a separate bam file, run

```
samtools view -bh -o example.17.bam example.sorted.bam \
  $(awk 'BEGIN{ORS=" "} NR>1 && $3==17 {print $2}' rsem_orig.tsv)
```

## Building

### Dependencies

* A C++-14 compiler
* [htslib](http://www.htslib.org/)
* [Boost](http://www.boost.org/)

### Compilation

```
CXXFLAGS=-O2 make
```

This will create an executable `bamcc` in the current directory.

## Testing

First, run `make`.

Then, run `./test`.
This will update all files `test_files/example*.tsv`.
Failures may appear either as messages from `bamcc` or differences in output
files reported by `git diff`.

## Security

Do not run this program on untrusted or potentially malformed input files.
