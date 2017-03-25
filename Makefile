LDFLAGS=-lhts
CXXFLAGS=-O2 -std=c++14
LINK.o = $(LINK.cc)

conn_comps: conn_comps.o
