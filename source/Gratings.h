// Gratings.h - defines the diffraction depending on which grating you're at

#ifndef GRATINGS_H
#define GRATINGS_H

void ( * get_initial_intensity(double z, double Grat3x[], double intensity_array[]));

/* The following parameters have the same names as the variables defined in the file Misc.h:
 * r0 = sp.initial_radius_of_wavefront_curvature 
 * el0 = sp.initial_coherence_width
 * w0 = sp.initial_beamwidth
 */
void ( * intensity_after_1st_grating(double current_z_position,double r0,double el0, double w0, double Grat3x[], double intensity_array[]));

void ( * intensity_after_2nd_grating(double current_z_position, double el1x, double w1x, double r1x, double Grat3x[], double intensity_array[]));


#endif // end Gratings.h
