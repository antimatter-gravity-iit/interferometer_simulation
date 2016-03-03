#	Makefile
#
#	Assumes Root is installed in $HOME/root.
#
#	TARGETS:
#		all	build the usual program (default)
#		sptest	build a test program for SimplePlot
#
#	NOTE: This has only been tested on Mac OS X 10.9.5 (Mavericks).
CC = gcc-5
CXX = g++-5

# list source files here 
SRC = *.c

#// for MAC:  ROOT = $(HOME)/Physics/root
ROOT = $(HOME)/Downloads/root

CFLAGS = $(shell $(ROOT)/bin/root-config --cflags)
LIBS =   $(shell $(ROOT)/bin/root-config --libs)
RPATH = -Wl,-rpath,$(ROOT)/lib

all:	program sptest

program: $(SRC) SimplePlot.cpp SimplePlot.h Gratings.c Gratings.h BeamParams.c BeamParams.h PhaseShifts.c PhaseShifts.h Misc.c Misc.h
	$(CXX) $(CFLAGS) -o program  $(SRC) SimplePlot.cpp Misc.c BeamParams.c Gratings.c PhaseShifts.c $(RPATH) $(LIBS)
#//-Wall
sptest:	SimplePlot.h SimplePlot.cpp
	$(CXX) $(CFLAGS) -DTEST -o sptest SimplePlot.cpp $(RPATH) $(LIBS)

clean:
	rm -f *.o program sptest
