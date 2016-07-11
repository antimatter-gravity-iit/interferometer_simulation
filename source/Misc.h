/* 
 * Misc.h
 * Copyright (C) 2016 Antimatter Gravity Interferometer Group, Illinois Institute of Technology (IIT). 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * 				*
 *
 * Code inspired by thesis by Dr. Benjamin McMorran
 * Electron Diffraction and Interferometry Using Nanostructures
 * http://gradworks.umi.com/33/52/3352633.html
 *
 * For a list of collaborators see README.md and CREDITS.
 *				
 *				*
 * 
 * DESCRIPTION:
 * Header file for Misc.c. Contains structure declarations.
 *
 */

#ifndef MISC_H
#define MISC_H

#include "PhaseShifts.h"
#include "Gratings.h"
#include "BeamParams.h"

// global sp structure, contains all of the simulation parameters that can be modified. 
typedef struct {
int account_gravity;				// account for gravitational forces [1] yes, [0] no.
int account_van_der_waals;	
double particle_velocity;			// Velocity of beam particles.
int simulation_option;				// choose whether or not to have a full simulation, or a final simulation interference pattern  
int logchoice;					// whether or not to scale interference pattern with a log base in order to see smaller intensities
double grating1_open_fraction;                  // G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it.
double grating2_open_fraction;	                // G2 open fraction; how open the second grating is.
double grating_period;                  	// 'Grating period' means the distance between two consecutive slit openings.
double initial_radius_of_wavefront_curvature;   // Corresponds to r0 in older versions.
double initial_coherence_width;                 // Corresponds to el0 in older versions.
double initial_beamwidth;                      	// Corresponds to w0 in older versions. Muonium beam width.
double z_position_1st_grating;                  // In m. Distance of 1st grating to z axis origin (z=0).
double z_position_2nd_grating;                  // In m. Distance of 2nd grating to z axis origin (z=0).
double G2_x;		                   	// Initial lateral offset of G2.
double twist_angle;		                   	// could be 0.05 or more. This is the twist between 1st and second gratings, in degrees. 2nd and 3rd grating are fixed to same rotational twist.
double grating_thickness;		        // In m. Used for Van der Waals interactions: see PhaseShifts.c.
double wedge_angle;	                     	// Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
double tilt_angle;                            	// A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
double resolution;	                        // Plot resolution in pixels.
double z_start;		                      	// defining bounds of the grating structure? This definitely corresponds to a location, probably the bottom of the grating.
double z_end;	                          	// probably the end of the grating.
double x_start;		                   	// x position start, -200 microns
double x_end;					// x position ends, 200 microns
// double y_start;		                // y position start, 110 microns
// double y_end;				// y position ends, 110 microns
double slit_height;				// maximum z height
double wavelength;				// beam's wavelength
double intensity_cutoff;			// =.000001  Intensity values below cutoff value are truncated or neglected (treated as 0). Can also be 5e-5 like in McMorran thesis. Or 0.001.
int number_of_rows_fourier_coefficient_array; 	// Rows of real_part_fourier_coefficient_array and imaginary_part_fourier_coefficient_array arrays; used to calculate phase shift.    
}simparam;

extern simparam sp;

double maximumvalue(double arr[], int array_size); 
// obtains the max value in an array

double sinc(double x); 
// returns sin(x)/x

double ( *ixgenerator(double a[], double current_z_position, int logchoice)); 
// accounts for decay (if applicable) and normalizing the intensity scale

int find_element_position_in_array(int value, int *arr); 
// checks for a certain x-position in the array

#endif
