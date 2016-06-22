// BeamParams.c - contains the functions defining the GSM (Gaussian-Schell Model) beam behavior

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string> 

#include "BeamParams.h"
#include "PhaseShifts.h"
#include "Misc.h"
#include "Gratings.h"

double zp(double z, double v) {
// compute magnification factor due to wavefront curvature
    double zp;
    zp = (v * z)/(z + v); 
    // v = GSM radius of wavefront curvature, z = current z location; currently, zp is only called with v = r0 (defined as a global variable) above. 
    return(zp);
}

/* The following are the same variables that are defined in the simulation parameters structure in the file Misc.h:
 * r0 = initial_radius_of_wavefront_curvature
 * el0 = initial_coherence_width
 * w0 = initial_beamwidth
 */
double calculate_width(double z, double r0, double el0, double w0, double which_width) {
    // Computes GSM beam width and beam coherence width (GSM = Gaussian-Schell Model of gratings)
    double lambda = sqrt(1.5e-18 / sp.energy);
    double width;
    width = which_width * fabs(z / zp(z,r0)) * sqrt(1 + pow(((lambda * zp(z,r0))/(el0 * w0)),2)); 
    return(width);
}

double v(double z,double r0, double el0, double w0) {
    // compute GSM radius of wavefront curvature
    double lambda = sqrt((1.5 * pow(10,-18))/(sp.energy)); 
    // value containing wavelength; this is approx. wavelength of xrays.
    double v;
    v=(z)/(1-zp(z,r0)/(z * (1 + pow(((lambda * zp(z,r0)/(el0 * w0))),2)))); 
    return(v);
}
