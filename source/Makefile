#	Makefile
#
#	Assumes Root is installed in $HOME/root.
#
#	TARGETS:
#		all	build the usual program (default)
#		sptest	build a test program for SimplePlot
#
#	NOTE: This has only been tested on Mac OS X 10.9.5 (Mavericks).

# list source files here 
SRC = *.c

ROOT = $(HOME)/Physics/root

CFLAGS = $(shell $(ROOT)/bin/root-config --cflags)
LIBS =   $(shell $(ROOT)/bin/root-config --libs)
RPATH = -Wl,-rpath,$(ROOT)/lib

all:	program sptest

program: $(SRC) SimplePlot.cc SimplePlot.hh Gratings.cc Gratings.hh BeamParams.cc BeamParams.hh PhaseShifts.cc PhaseShifts.hh Misc.cc Misc.hh
	g++ $(CFLAGS) -o program  $(SRC) SimplePlot.cc Misc.cc BeamParams.cc Gratings.cc PhaseShifts.cc $(RPATH) $(LIBS)
#//-Wall
sptest:	SimplePlot.hh SimplePlot.cc
	g++ $(CFLAGS) -DTEST -o sptest SimplePlot.cc $(RPATH) $(LIBS)


clean:
	rm -f *.o program sptest
