// Gratings.h - defines the diffraction depending on which grating you're at

#ifndef GRATINGS_H
#define GRATINGS_H

double ( * gp0(double z,double r0,double el0, double w0, double Grat3x[], double Grat3I[],double energy, int rows, int col, int xpnts));//prototype

double ( * gp1(double z12,double r0,double el0, double w0, double Grat3x[], double Grat3I[],double energy, int rows, int col, int elecOrAtom, double vel, int xpnts, double width, double abszloc, int accountGrav));//prototype

double ( * gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x, double Grat3x[], double Grat3I[],double energy, int rows, int col, int elecOrAtom, double vel, int xpnts, double width, double abszloc, int accountGrav));//prototype



#endif // end Gratings.hh
