/* 
 * BeamParams.h
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
 * Header file for BeamParams.c.
 *
 */ 

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

double calculate_wavefront_radius(double z,double r0, double el0, double w0); 	// Compute GSM radius of wavefront curvature

#endif
