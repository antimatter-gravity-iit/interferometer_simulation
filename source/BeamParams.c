/* 
 * BeamParams.c
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
 * This file contains the functions defining the GSM (Gaussian-Schell Model) beam behavior.
 *
 */

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
    double width;
    width = which_width * fabs(z / zp(z,r0)) * sqrt(1 + pow(((sp.wavelength * zp(z,r0))/(el0 * w0)),2)); 
    return(width);
}

double calculate_wavefront_radius(double z,double r0, double el0, double w0) {
    // compute GSM radius of wavefront curvature
    double v;
    v=(z)/(1-zp(z,r0)/(z * (1 + pow(((sp.wavelength * zp(z,r0)/(el0 * w0))),2)))); 
    return(v);
}
