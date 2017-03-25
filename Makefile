LDFLAGS=-lhts -ltcmalloc
CXXFLAGS=-O2 -gdwarf
LINK.o = $(LINK.cc)

conn_comps: conn_comps.o
