// PhaseShifts.hh - contains the function(s) that contribute to different phase shifts in this system

#ifndef PHASESHIFTS_HH
#define PHASESHIFTS_HH

double ( * ReTandImTgenerator(double ReTorImTar[], double energy, int elecOrAtom, int ReTorImT, double vel, double width, double abszloc, int accountGrav)); // the function responsible for calculating the real and imaginary phase shifts of various physical effects


#endif // end of file
