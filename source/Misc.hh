// Misc.hh - contains functions of other miscellaneous tasks, such as taking care of boundary conditions, checking the max value in an array, etc.

#ifndef MISC_HH
#define MISC_HH

double sinc(double x); // returns sin(x)/x

int x2pnts(int value, int  * arr); // checks for a certain x-position in the array

double maximumvalue(double arr[], int rows); // obtains the max value in an array

double ( * ixgenerator(double a[], double zloc, int logchoice, int rows)); // accounts for decay (if applicable) and normalizing the intensity scale

#endif
