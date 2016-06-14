#	Makefile
#
#	IMPORTANT: The ROOT symbolic constant has to be set to the path where ROOT is installed.
#
#	TARGETS:
#		all		builds the usual program plus the test program (default)
#		program 	builds the usual program only
#		sptest		builds a test program for SimplePlot Misc.c BeamParams.c Gratings.c PhaseShifts.c
#		debug_program	builds usual program with added debugging information (necessary for gdb) and extra warnings
#		debug_sptest	builds test program with added debugging information and extra warnings
#		debug		builds all with debugging information
#		clean		removes output files (executables)
#
CC = gcc-5
CXX = g++-5

# list source files here 
SRC = *.c

# Lucas Rios's computer
# ROOT = $(HOME)/Downloads/build

# Yuri Tonin's computer
# ROOT = $(HOME)/InstalledRoot

# Lucas Abrantes's computer
ROOT = /opt/root/build

CFLAGS = $(shell $(ROOT)/bin/root-config --cflags)
LIBS =   $(shell $(ROOT)/bin/root-config --libs)
RPATH = -Wl,-rpath,$(ROOT)/lib

all:	program sptest

program: $(SRC) SimplePlot.cpp SimplePlot.h Gratings.c Gratings.h BeamParams.c BeamParams.h PhaseShifts.c PhaseShifts.h Misc.c Misc.h
	$(CXX) $(CFLAGS) -fext-numeric-literals -o program  $(SRC) SimplePlot.cpp $(RPATH) $(LIBS) 	

sptest:	SimplePlot.h SimplePlot.cpp
	$(CXX) $(CFLAGS) -DTEST -fext-numeric-literals -o sptest SimplePlot.cpp $(RPATH) $(LIBS)

# -Wall and -Wextra are flags that enable "warnings about constructions that some users consider questionable".
# See https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html

debug_program: $(SRC) SimplePlot.cpp SimplePlot.h Gratings.c Gratings.h BeamParams.c BeamParams.h PhaseShifts.c PhaseShifts.h Misc.c Misc.h
	$(CXX) -g $(CFLAGS) -fext-numeric-literals -Wall -Wextra -o program_with_debugging_info  $(SRC) SimplePlot.cpp $(RPATH) $(LIBS)

debug_sptest: SimplePlot.h SimplePlot.cpp
	$(CXX) -g $(CFLAGS) -DTEST -fext-numeric-literals -Wall -Wextra -o sptest_with_debugging_info SimplePlot.cpp $(RPATH) $(LIBS)

debug: debug_program debug_sptest

clean:
	rm -f *.o program sptest program_with_debugging_info sptest_with_debugging_info
