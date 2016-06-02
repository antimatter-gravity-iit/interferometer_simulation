// BeamParams.h hosts the functions that correspond to the GSM Model presented in Dr. Ben McMorran's thesis work, located: 

#ifndef BEAMPARAMS_H
#define BEAMPARAMS_H


double zp(double z, double v); //// compute magnification factor due to wavefront curvature

//double w(double z,double r0, double el0, double w0, double energy);  // Compute GSM beam width (GSM = Gaussian-Schell Model of gratings)
double w(double z,double r0, double el0, double w0);  // Compute GSM beam width (GSM = Gaussian-Schell Model of gratings)

//double el(double z,double r0, double el0, double w0, double energy);  // computes GSM beam coherence width
double el(double z,double r0, double el0, double w0);  // computes GSM beam coherence width

//double v(double z,double r0, double el0, double w0, double energy); // compute GSM radius of wavefront curvature
double v(double z,double r0, double el0, double w0); // compute GSM radius of wavefront curvature


#endif // end of the file
