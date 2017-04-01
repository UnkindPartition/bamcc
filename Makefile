LDFLAGS += -lhts
CXXFLAGS += -std=c++14 -Wall
LINK.o = $(LINK.cc)

bamcc: bamcc.o
