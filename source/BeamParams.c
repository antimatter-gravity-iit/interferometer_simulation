// BeamParams.c - contains the functions defining the GSM (Gaussian-Schell Model) beam behavior

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string> 

#include "BeamParams.h"
#include "PhaseShifts.h"
#include "Misc.h"
#include "Gratings.h"

double zp(double z, double v) 
{
// compute magnification factor due to wavefront curvature
    double zp;
    zp = (v * z)/(z + v); 
    // v = GSM radius of wavefront curvature, z = current z location; currently, zp is only called with v = r0 (defined as a global variable) above. 
    return(zp);
}

double w(double z,double r0, double el0, double w0, double energy) 
{
    // Compute GSM beam width (GSM = Gaussian-Schell Model of gratings)
    double lambda = sqrt((1.5 * pow(10,-18))/(energy)); 
    // variable containing value of the wavelength; wavelength of what?
    double w;
    w = (el0) * (fabs((z)/(zp(z,r0))) * ((sqrt((1 + (pow(((lambda * zp(z,r0))/(el0 * w0)),2))))))); 
    // what does this correspond to?
    return(w);
}

double el(double z,double r0, double el0, double w0, double energy) 
{
    // GSM beam coherence width
    double lambda = sqrt((1.5 * pow(10,-18))/(energy)); 
    // again, really small value; should this be global? value containing the wavelength; wavelength of what?
    double w;
    w = (el0) * (fabs((z)/(zp(z,r0))) * ((sqrt((1 + (pow(((lambda * zp(z,r0))/(el0 * w0)),2))))))); 
    return(w);
}

double v(double z,double r0, double el0, double w0, double energy) 
{
    // compute GSM radius of wavefront curvature
    double lambda = sqrt((1.5 * pow(10,-18))/(energy)); 
    // value containing wavelength; this is approx. wavelength of xrays.
    double v;
    v=(z)/(1-zp(z,r0)/(z * (1 + pow(((lambda * zp(z,r0)/(el0 * w0))),2)))); 
    return(v);
}
