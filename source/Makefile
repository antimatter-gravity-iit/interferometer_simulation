#	Makefile
#
#	Assumes Root is installed in $HOME/root.
#
#	TARGETS:
#		all	build the usual program (default)
#		sptest	build a test program for SimplePlot Misc.c BeamParams.c Gratings.c PhaseShifts.c
#
#	NOTE: This has only been tested on Mac OS X 10.9.5 (Mavericks).
CC = gcc-5
CXX = g++-5

# list source files here 
SRC = *.c

# for MAC computer:  
# ROOT = $(HOME)/Physics/root
ROOT = $(HOME)/Downloads/root

CFLAGS = $(shell $(ROOT)/bin/root-config --cflags)
LIBS =   $(shell $(ROOT)/bin/root-config --libs)
RPATH = -Wl,-rpath,$(ROOT)/lib

all:	program sptest

program: $(SRC) SimplePlot.cpp SimplePlot.h Gratings.c Gratings.h BeamParams.c BeamParams.h PhaseShifts.c PhaseShifts.h Misc.c Misc.h
	$(CXX) $(CFLAGS) -fext-numeric-literals -o program  $(SRC) SimplePlot.cpp $(RPATH) $(LIBS) 	
#-Wall -Wextra  //extra debug code 

sptest:	SimplePlot.h SimplePlot.cpp
	$(CXX) $(CFLAGS) -DTEST -fext-numeric-literals -o sptest SimplePlot.cpp $(RPATH) $(LIBS)

clean:
	rm -f *.o program sptest
