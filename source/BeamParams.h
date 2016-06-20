// BeamParams.h hosts the functions that correspond to the GSM Model presented in Dr. Ben McMorran's thesis work, located: 

#ifndef BEAMPARAMS_H
#define BEAMPARAMS_H


double zp(double z, double v); 				// Compute magnification factor due to wavefront curvature

/* Computes GSM beamwidth or coherence width, depending on the last argument (GSM = Gaussian-Schell model of gratings).
 * The following are the same variables that are defined in the simulation parameters structure in the file Misc.h:
 * r0 = initial_radius_of_wavefront_curvature 
 * el0 = initial_coherence_width
 * w0 = initial_beamwidth
 */
double calculate_width(double z, double r0, double el0, double w0, double which_width);

double v(double z,double r0, double el0, double w0); 	// Compute GSM radius of wavefront curvature


#endif // end of the file
