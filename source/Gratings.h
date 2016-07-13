/* 
 * Gratings.h
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
 * Header file for Gratings.c.
 *
 */

#ifndef GRATINGS_H
#define GRATINGS_H

void ( * get_initial_intensity(double z, double x_positions_array[], double intensity_array[]));

/* The following parameters have the same names as the variables defined in the file Misc.h:
 * r0 = sp.initial_radius_of_wavefront_curvature 
 * el0 = sp.initial_coherence_width
 * w0 = sp.initial_beamwidth
 */
void ( * intensity_after_1st_grating(double current_z_position,double r0,double el0, double w0, double x_positions_array[], double intensity_array[]));

void ( * intensity_after_2nd_grating(double current_z_position, double el1x, double w1x, double r1x, double x_positions_array[], double intensity_array[]));

#endif
