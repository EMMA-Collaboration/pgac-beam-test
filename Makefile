# Makefile

CXXFLAGS = -g -O2 -Wall -Wuninitialized -Wno-write-strings -std=c++0x 

# optional ZLIB library

CXXFLAGS += -DHAVE_ZLIB

# ROOT libraries

ifdef ROOTSYS
#ROOTLIBS  = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)  -lXMLParser -lThread
ROOTGLIBS = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lThread
#RPATH    += -Wl,-rpath,$(ROOTSYS)/lib
CXXFLAGS += -DHAVE_ROOT $(shell $(ROOTSYS)/bin/root-config --cflags)
endif

ifdef MIDASSYS
MIDASLIBS = $(MIDASSYS)/linux/lib/libmidas.a -lutil -lrt
CXXFLAGS += -DHAVE_MIDAS -DOS_LINUX -Dextname -I$(MIDASSYS)/include

#UNAME=$(shell uname)
#ifeq ($(UNAME),Darwin)
#CXXFLAGS += -DOS_LINUX -DOS_DARWIN
#MIDASLIBS = $(MIDASSYS)/darwin/lib/libmidas.a
#RPATH=
#endif

endif

ifdef ROOTSYS
CXXFLAGS += -DHAVE_LIBNETDIRECTORY -IlibNetDirectory
endif

ROOTANAINC = -I../include
ROOTANALIBS = ../lib/librootana.a

OBJS:=
OBJS += TV792Histogram.o TV1190Histogram.o
OBJS += TL2249Histogram.o TAgilentHistogram.o

all: $(OBJS) ana.exe anaDisplay.exe ana_bare.exe ana_tree.exe

ana.exe: ana.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl -lrt -lutil

anaDisplay.exe: anaDisplay.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl -lrt -lutil

ana_bare.exe: ana_bare.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl -lrt -lutil

ana_tree.exe: ana_tree.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl -lrt -lutil

%Dict.o: %Dict.cxx
	$(CXX) $(CXXFLAGS) $(ROOTANAINC) -c $<

%Dict.cxx: %.h %_LinkDef.h
	LD_LIBRARY_PATH=$(ROOTSYS)/lib $(ROOTSYS)/bin/rootcint -f $@ -c -p $(CXXFLAGS) $(ROOTANAINC) $^


%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(ROOTANAINC) -c $<

dox:
	doxygen

clean::
	-rm -f *.o *.a *.exe 

# end

