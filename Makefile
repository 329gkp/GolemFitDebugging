CFLAGS=-O0 -fPI

CXXFLAGS+=-std=c++11 -fPIC
CXXFLAGS+=-I. -O3 -Wall -Wunused-but-set-variable -Wunused-variable

ifeq (${PREFIX},)
PREFIX=/usr/local
endif

ifeq (${SROOT},)
SROOT=/usr/local
endif
CXXFLAGS+=-I${SROOT}/include
LDFLAGS+=-L${SROOT}/lib
LDFLAGS+=-L${SROOT}/lib64

ifeq (${GOLEMBUILDPATH},)
GOLEMBUILDPATH=/usr/local
endif
CXXFLAGS+=-I${GOLEMBUILDPATH}/include/
CXXFLAGS+=-I${PREFIX}/include/
CXXFLAGS+=-DPHOTOSPLINE_INCLUDES_SPGLAM
LDFLAGS+=-L${GOLEMBUILDPATH}/lib/
LDFLAGS+=-L${PREFIX}/lib/
CXXFLAGS+=-DPHOTOSPLINE_INCLUDES_SPGLAM
LDFLAGS+=-L${GOLEMBUILDPATH}/lib/
#LDFLAGS+=-L/net/cvmfs_users/carguelles/lib

LDFLAGS+=-lgsl -lgslcblas -lm
LDFLAGS+=-lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LDFLAGS+=-lhdf5 -lhdf5_hl
LDFLAGS+=-pthread

# the following should be found in the include and lib directories
# of the SNOTBUILDPATH
LDFLAGS+=-lNewNuFlux
LDFLAGS+=-lPhysTools
LDFLAGS+=-lSQuIDS -lnuSQuIDS
LDFLAGS+=-lLeptonWeighter
LDFLAGS+=-lphotospline -lcfitsio
LDFLAGS+=-lnuFATE
LDFLAGS+=-lhealpix_cxx -lcxxsupport
LDFLAGS+=-lpthread

GOLEM_HELPER_OBJECTS = Event.o analysisWeighting.o utils.o compactIO.o GolemFit.o oversizeWeight.o self_veto.o GolemParameters.o GolemTools.o GolemEnumDefinitions.o DarkMatterInteractions.o Earth.o
GOLEM_OBJECT = GolemFitTest.o

OS_NAME=$(shell uname -s)

ifeq (${OS_NAME},Linux)
DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))
endif

ifeq (${OS_NAME},Darwin)
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)
endif

# Project files
NAME=GolemFit
VERSION=1.0.0
STAT_PRODUCT=../lib/lib$(NAME).a
DYN_PRODUCT=../lib/lib$(NAME)$(DYN_SUFFIX)
OBJECTS= ${GOLEM_HELPER_OBJECTS}

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT) GolemFitTestHESE GolemFitTestSTERILE

.PHONY: all clean

clean:
	rm -rf GolemFitTestHESE GolemFitTestSTERILE $(GOLEM_HELPER_OBJECTS) $(DYN_PRODUCT) $(STAT_PRODUCT) *.o

Earth.o : Earth.cpp Earth.h
	$(CXX) $(CXXFLAGS) Earth.cpp -c -o Earth.o

DarkMatterInteractions.o : DarkMatterInteractions.cpp DarkMatterInteractions.h
	$(CXX) $(CXXFLAGS) DarkMatterInteractions.cpp -c -o DarkMatterInteractions.o

analysisWeighting.o : analysisWeighting.cpp analysisWeighting.h
	$(CXX) $(CXXFLAGS) analysisWeighting.cpp -c -o analysisWeighting.o

GolemParameters.o : GolemParameters.cpp GolemParameters.h
	$(CXX) $(CXXFLAGS) GolemParameters.cpp -c -o GolemParameters.o

compactIO.o : compactIO.h compactIO.cpp Event.h analysisWeighting.h
	$(CXX) $(CXXFLAGS) compactIO.cpp -c -o compactIO.o

Event.o : Event.cpp Event.h json.hpp
	$(CXX) $(CXXFLAGS) Event.cpp -c -o Event.o

GolemTools.o : GolemTools.cpp GolemTools.h
	$(CXX) $(CXXFLAGS) GolemTools.cpp -c -o GolemTools.o

GolemEnumDefinitions.o : GolemEnumDefinitions.cpp GolemEnumDefinitions.h
	$(CXX) $(CXXFLAGS) GolemEnumDefinitions.cpp -c -o GolemEnumDefinitions.o

self_veto.o : self_veto.cpp self_veto.h
	$(CXX) $(CXXFLAGS) self_veto.cpp -c -o self_veto.o

utils.o : utils.cpp utils.h
	$(CXX) $(CXXFLAGS) utils.cpp -c -o utils.o

oversizeWeight.o : oversizeWeight.cpp oversizeWeight.h
	$(CXX) $(CXXFLAGS) oversizeWeight.cpp -c -o oversizeWeight.o

GolemFit.o : GolemFit.cpp Event.h analysisWeighting.h compactIO.h GolemFit.h oversizeWeight.h GolemParameters.h GolemEnumDefinitions.h DarkMatterInteractions.h FastMode.h
	$(CXX) $(CXXFLAGS) GolemFit.cpp -c -o GolemFit.o

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking dynamic library $(DYN_PRODUCT)
	@mkdir -p ../lib/
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking static library $(STAT_PRODUCT)
	@mkdir -p ../lib/
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

GolemFitTestHESE: $(GOLEM_HELPER_OBJECTS) GolemFitTestHESE.o
	$(CXX) $(CXXFLAGS) GolemFitTestHESE.cpp $(LDFLAGS) $(GOLEM_HELPER_OBJECTS) -o GolemFitTestHESE

GolemFitTestSTERILE: $(GOLEM_HELPER_OBJECTS) GolemFitTestSTERILE.o
	$(CXX) $(CXXFLAGS) GolemFitTestSTERILE.cpp $(LDFLAGS) $(GOLEM_HELPER_OBJECTS) -o GolemFitTestSTERILE

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/GolemFit
	@mkdir -p $(PREFIX)/include/GolemFit
	@cp ./*.h $(PREFIX)/include/GolemFit
	@cp ./*.hpp $(PREFIX)/include/GolemFit
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp golemfit.pc $(PREFIX)/lib/pkgconfig

