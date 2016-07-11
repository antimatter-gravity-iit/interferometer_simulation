/* 
 * Gratings.c
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
 * This file holds the functions that compute the diffraction after each of the gratings in the system.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cmath>
#include <time.h>

#include "Gratings.h"
#include "BeamParams.h"
#include "Misc.h"
#include "PhaseShifts.h"

#include <complex.h>
#include <limits.h> 		//fix for islimit

void ( * get_initial_intensity(double z, double x_positions_array[], double intensity_array[]))
// get intensity profile 
{
	clock_t start, end; //starting a timer to get the time spent in function get_initial_intensity
	start = clock();
	double diff=0;


	double w1;

	// Width of incoming beam, calculated from initial parameters.
	w1 = calculate_width(z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_beamwidth); 

	for(int i=0; i<sp.resolution; i++) {
		intensity_array[i]=exp(-(M_PI * pow((x_positions_array[i]/w1),2))); 		// a[i][1] is the intensity of the beam at the xposition at step i.
	}

	end = clock();
	diff =((double) (end - start))/ CLOCKS_PER_SEC;
	printf("Time elapsed: %f seconds\n",diff); 		//printing the time spent in get_initial_intensity.
}


void ( * intensity_after_1st_grating(double current_z_position,double el1, double w1, double r1, double x_positions_array[], double intensity_array[]))
// get intensity profile after one grating
{
	clock_t start, end; //starting a timer to get the time spent in function intensity_after_1st_grating
	start = clock();
	double diff=0;

	double current_z_distance_to_1st_grating    = current_z_position - sp.z_position_1st_grating;	//z location between 1st and 2nd gratings
	double coefficient; 
    	double diffraction_orders  = 5;	
	/*
	 * Explanation of variables:
	 * w2 = GSM width of beam after the first grating.
	 * r2 = radius of GSM wavefront curvature after the first grating
	 * el2 = GSM beam coherence width after the first grating.
    	 */
	double w2  = calculate_width(current_z_distance_to_1st_grating, r1, el1, w1, w1);		// width of beam between z1 and z2 (after grating 1)
	double r2  = calculate_wavefront_radius(current_z_distance_to_1st_grating, r1, el1, w1);	// radius of wavefront curvature between grating 1 and 2
	double el2 = calculate_width(current_z_distance_to_1st_grating, r1, el1, w1, el1);		// beam coherence width
    	
	int pos[41]={0};

	for (int i=0; i<sp.number_of_rows_fourier_coefficient_array; i++) // since sp.number_of_rows_fourier_coefficient_array is currently 41, pos[i] = -20 to 20.
	{
	pos[i]=i-((sp.number_of_rows_fourier_coefficient_array-1)/2);
	}

	double real_part_fourier_coefficient_array[41]={0};
	real_and_imaginary_arrays_generator(real_part_fourier_coefficient_array, 1, current_z_position); // calculates phase shift where 1 is to consider real components
		
	double imaginary_part_fourier_coefficient_array[41]={0};
	real_and_imaginary_arrays_generator(imaginary_part_fourier_coefficient_array, 2, current_z_position); // calculates phase shift where 2 is to consider real components

	for (int i=0; i<sp.resolution; i++) { 
		for (int n1=-diffraction_orders; n1<=diffraction_orders; n1++) {
			for (int n2=-diffraction_orders; n2<=diffraction_orders; n2++) {
				double delta_n =n1-n2; 
				// TODO: explain what average_n, average_m, delta_m, delta_n are. LR, Y 
				double average_n = (n1 + n2)/2;

				if (sp.account_gravity == 0 && sp.account_van_der_waals == 0)
				{ // if useimagecharge = 0, ignore image charge effects at G1. 
					coefficient = sinc(sp.grating1_open_fraction * M_PI * n1)  *  (sinc(sp.grating1_open_fraction * M_PI * n2) * pow((sp.grating1_open_fraction), 2));
				}

				else
				{ // if usechargeimage = 1, don't ignore image charge effects at G1.
					coefficient = real_part_fourier_coefficient_array[find_element_position_in_array(n1, (int * )pos)] * real_part_fourier_coefficient_array[find_element_position_in_array(n2,(int * )pos)] + imaginary_part_fourier_coefficient_array[find_element_position_in_array(n1,(int * )pos)] * imaginary_part_fourier_coefficient_array[find_element_position_in_array(n2,(int * )pos)];
				}
				
				coefficient = coefficient * exp(-M_PI * pow((delta_n * sp.wavelength * current_z_distance_to_1st_grating)/(sp.grating_period * el2),2));
				// added isfinite macro in order to avoid inf values

				if (std::isfinite(coefficient)==0 || coefficient < sp.intensity_cutoff) { // if coefficient is infinite, then:
					coefficient=0;
				}
			      	else { // if coefficient ends up larger than cutoff value, add the values to the current a[i][1]'s intensities.
					intensity_array[i] +=   coefficient * exp(-M_PI * pow(((x_positions_array[i]-average_n * sp.wavelength * current_z_distance_to_1st_grating/sp.grating_period)/w2),2)) * cos(2 * M_PI * (delta_n/sp.grating_period) * (x_positions_array[i]-average_n * sp.wavelength * current_z_distance_to_1st_grating/sp.grating_period) * (1-current_z_distance_to_1st_grating/r2));
				    // Since a[i][1] etc. is actually the ix array, and arrays essentially get passed by reference, this is modifying the ix array.
				    continue;
				}
			}
		}
	}

	end = clock();
	diff =((double) (end - start))/ CLOCKS_PER_SEC;
	printf("Time elapsed: %f seconds\n",diff); //printing the time spent in intensity_after_1st_grating.

}

void ( * intensity_after_2nd_grating(double current_z_position, double el1x, double w1x, double r1x, double x_positions_array[], double intensity_array[]))
{	
	clock_t start, end; //starting a timer to get the time spent in function intensity_after_2nd_grating
	start = clock();
	double diff=0;

	// get intensity profile after two grating
	double G2_x    = sp.G2_x; //x position after second grating
	double r1y     = r1x;
	double w1y     = w1x;
	double el1y    = el1x;
	double current_z_distance_to_2nd_grating     = current_z_position - sp.z_position_2nd_grating;
    	double d1     = sp.grating_period;			// period = period of gratings
    	double d2     = sp.grating_period;
    	double current_z_distance_to_1st_grating    = current_z_position - sp.z_position_1st_grating; // z distance between grating 1 and 3
    	double diffraction_orders    = 5;
    	double _Complex coefficient;

    	/* THIS FUNCTION IS USING GSM MODEL FROM MCMORRAN, CRONIN 2008
    	 * IT MUST BE CHANGED IF YOU WANT TO USE OLDER MODEL FROM BREZGER 2003
	 */

    	double delta_n = 0;
    	double delta_m = 0;
    	double average_m  = 0;
    	double average_n  = 0;
    	int    central_index_1 = 0;
    	int    central_index_2 = 0;
    	int    central_index_3 = 0;
    	int    central_index_4 = 0;
        double argument_d;
        double argument_f;
        double argument_p;
        double argument_v;
	double argument_f_p;
        double function_d;
        double function_v;

	/* Explanation of variables:
	 * w3 = GSM width of beam after the second grating.
	 * r3 = radius of GSM wavefront curvature after the second grating
	 * el3 = GSM beam coherence width after the second grating.
    	 */
    	double el3x = calculate_width(current_z_distance_to_1st_grating, r1x, el1x, w1x, el1x);	// current_z_distance_to_1st_grating == G2z - G1z + z_start + 0 * zres; GSM coherence width in x-axis
    	double w3x  = calculate_width(current_z_distance_to_1st_grating, r1x, el1x, w1x, w1x); 	// Beam width in x-axis
    	double v3x  = calculate_wavefront_radius(current_z_distance_to_1st_grating,r1x,el1x,w1x);	// Gaussian-Schell Model (GSM) radius of wavefront curvature in x-axis
    	double el3y = calculate_width(current_z_distance_to_1st_grating, r1y, el1y, w1y, el1y); 	// Coherence width in y-axis
    	double w3y  = calculate_width(current_z_distance_to_1st_grating, r1y, el1y, w1y, w1y); 	// Beam width in y-axis
    	double v3y  = calculate_wavefront_radius(current_z_distance_to_1st_grating,r1y,el1y,w1y);	// radius of wavefront curvature in y-axis
    
    	int pos[41]={0}; // array of 41 elements

    	for (int i=0; i<sp.number_of_rows_fourier_coefficient_array; i++){     // sp.number_of_rows_fourier_coefficient_array = rows of real_part_fourier_coefficient_array and imaginary_part_fourier_coefficient_array arrays = 41 for now
    		pos[i]=i-((sp.number_of_rows_fourier_coefficient_array-1)/2);  // so this goes from -20 to 20
    	}
    
	double real_part_fourier_coefficient_array[41]={0};
	real_and_imaginary_arrays_generator(real_part_fourier_coefficient_array, 1, current_z_position); // calculates phase shift where 1 is to consider real components
		
	double imaginary_part_fourier_coefficient_array[41]={0};
	real_and_imaginary_arrays_generator(imaginary_part_fourier_coefficient_array, 2, current_z_position); // calculates phase shift where 2 is to consider real components 
	

	for (int i=0; i<sp.resolution; i++) {
		for (int m1 = -diffraction_orders; m1 <= diffraction_orders; m1++) {
			for (int m2 = -diffraction_orders; m2 <= diffraction_orders; m2++) {
				for (int n1 = -diffraction_orders; n1 <= diffraction_orders; n1++) {
					for (int n2 = -diffraction_orders; n2 <= diffraction_orders; n2++) {
						
						delta_n = n1-n2;
						average_n  = ((double)(n1 + n2))/2;
						delta_m = m1-m2;
						average_m  = ((double)(m1 + m2))/2;

						central_index_1 = ( find_element_position_in_array( m1, (int *)pos ) );
						central_index_2 = ( find_element_position_in_array( m2, (int *)pos ) );
						central_index_3 = ( find_element_position_in_array( n1, (int *)pos ) );
						central_index_4 = ( find_element_position_in_array( n2, (int *)pos ) );

						// 0 means ignore image charge effects, 1 means include image charge effects

						if (sp.account_gravity == 0 && sp.account_van_der_waals == 0) {
							coefficient = 	       sinc(sp.grating1_open_fraction * M_PI * m1) +  0 * _Complex_I;
							coefficient = coefficient * (sinc(sp.grating1_open_fraction * M_PI * m2)) +  0 * _Complex_I;

						}
						else { // assumes G1 is identical to G2
							coefficient = 	      (real_part_fourier_coefficient_array[central_index_1] + imaginary_part_fourier_coefficient_array[central_index_1] * _Complex_I); // 
							coefficient = coefficient * (real_part_fourier_coefficient_array[central_index_2] - imaginary_part_fourier_coefficient_array[central_index_2] * _Complex_I);
						}
					
							coefficient = coefficient * (real_part_fourier_coefficient_array[central_index_3] + imaginary_part_fourier_coefficient_array[central_index_3] * _Complex_I);
							coefficient = coefficient * (real_part_fourier_coefficient_array[central_index_4] - imaginary_part_fourier_coefficient_array[central_index_4] * _Complex_I);

						//argument_d corresponds to the argument of equation 18b from McMorran & Cronin 2008. Note that y=0
                        			argument_d = -M_PI*(pow( x_positions_array[i]-sp.wavelength*current_z_distance_to_2nd_grating*(average_n*cos(sp.twist_angle)/d2 + average_m*current_z_distance_to_1st_grating/(d1*current_z_distance_to_2nd_grating) ),2 )/pow(w3x,2) + pow((average_n*sin(sp.twist_angle)*sp.wavelength)/(d2*w3y),2));
                        			//argument_f corresponds to the argument of equation 18c from McMorran & Cronin 2008. Note that y=0
                        			argument_f = -2*M_PI * x_positions_array[i] * ((delta_n*cos(sp.twist_angle)/d2)*(1-current_z_distance_to_2nd_grating/v3x) + (delta_m/d1)*(1-current_z_distance_to_1st_grating/v3x));
                        			//argument_p corresponds to the argument of equation 18d from McMorran & Cronin 2008
                        			argument_p = (2*M_PI*sp.wavelength*current_z_distance_to_1st_grating*delta_m/d1)*(average_n*cos(sp.twist_angle)/d2 + average_m/d1)*(1-current_z_distance_to_1st_grating/v3x) +     (2*M_PI*sp.wavelength*current_z_distance_to_2nd_grating*delta_n/d2)*( (average_m*cos(sp.twist_angle)/d1) * (1-current_z_distance_to_1st_grating/v3x) - (average_n*current_z_distance_to_2nd_grating/d2)*(pow(cos(sp.twist_angle),2)/v3x)+pow(sin(sp.twist_angle),2)/v3y ) ;
                        			//argument_v corresponds to the argument of equation 18e from McMorran & Cronin 2008
                        			argument_v = -M_PI* pow(sp.wavelength*current_z_distance_to_2nd_grating* (delta_n*cos(sp.twist_angle)/d2 + delta_m*current_z_distance_to_1st_grating/(d1*current_z_distance_to_2nd_grating)),2)/pow(el3x,2) -M_PI*pow(delta_n*sin(sp.twist_angle)*sp.wavelength*current_z_distance_to_2nd_grating/(d2*el3y),2);

						
						if (((__real__ coefficient) >= sp.intensity_cutoff) || ((__imag__ coefficient) >= sp.intensity_cutoff)) {
						    
                        				function_d   = exp(argument_d);                      
                        				function_v   = exp(argument_v);
                        				argument_f_p = argument_f + argument_p;
                        
                        				intensity_array[i] +=  ((__real__ coefficient) * cos(argument_f_p) - (__imag__ coefficient) * sin(argument_f_p))*function_d*function_v ;
						}
					}
				}
	    		}
		}
		
		
		
    	}


	end = clock();
	diff =((double) (end - start))/ CLOCKS_PER_SEC;
	printf("Time elapsed: %f seconds\n",diff); //printing the time spent in intensity_after_2nd_grating.
}
