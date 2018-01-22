CXX = g++
CXXFLAGS = -g -O2 -Wall -fPIC

# --- ROOT --------------------------------------------------------------
CXXFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lThread

# --- EXOSoftware -------------------------------------------------------
#CXXFLAGS +=-I/home/exodaq/EXOSoftware/include
#EXOLIBS = -L/home/exodaq/EXOSoftware/lib -lEXOUtilities

INCLUDE :=-I/nfs/slac/g/exo/software/builds/current/include
CXXFLAGS += $(INCLUDE)
EXOLIBS = -L/nfs/slac/g/exo/software/builds/current/lib -lEXOUtilities

all: AlphaSpectroscopy

AlphaSpectroscopy: AlphaSpectroscopy.o
	@ echo "Linking $@..."
	@ ${CXX} ${CXXFLAGS} -o $@ $^ ${ROOTLIBS} ${EXOLIBS}

AlphaSpectroscopy.o: AlphaSpectroscopy.cc
	@ ${CXX} ${CXXFLAGS} -c $^ -o $@

.PHONY : clean

clean:
	@echo "Cleaning up..."
	@ rm -f AlphaSpectroscopy AlphaSpectroscopy.o
