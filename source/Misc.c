/* 
 * Misc.c
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
 * This file contains functions responsible for miscellaneous tasks, such as taking care of boundary conditions,
 * checking which is the maximum value in an array, etc.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "Misc.h"

double maximumvalue(double arr[], int array_size){    
    // finds the max value of the array
    double m =0;
    for (int i = 0; i < array_size; i++ ){
        if (m<arr[i]){
            m = arr[i];// finds the maximum value of the array
        }
    }
    return m;
}


double sinc(double x){
    // this function avoids a division by x=0, sinc is defined to be sin(x)/x which will definied as 1 in such case
    double sinc_value=1;
    if (x!=0){
        sinc_value = (sin(x))/x;
    }
    return(sinc_value);
}


/* 
 * The function 'ixgenerator' is the x direction intensity calculator. It normalizes and compares intensities
 * to cutoff value, then determines which value to input to the array of x intensitites. Its arguments are:
 * 	intensity at the time of calculation (an initially empty array);
 * 	zlocation;
 * 	choice of scale for the plot.
 * The function only modifies the array fed to it; it doesn't return any value. 
 */
double *ixgenerator(double a[], double current_z_position, int logchoice) 
{
    double mu_lifetime = 0.0000022; 
    // half life average decay time of a muon, in seconds
    double max = maximumvalue(a, sp.resolution); 
    // maximum intensity value for a slice in z
    double min = sp.intensity_cutoff;
    // the minimum value the intensity can be at before being set to 0; cutoff is a global variable

    for (int i=0; i<sp.resolution; i++) {
        /*
	 * HALF-LIFE DECAY - means normalizing the intensity at the point i to what the intensity should be after a certain time t.
	 * This will be approximated by substituting distance, zloc, for t; x = vt, v = sp.particle_velocity, so t = x/sp.particle_velocity.
	 * In this case it's z instead of x. Since the electrons should become more and more as the antimuons decay, we actually want
	 * the electron presence to start out small, then go up.
	 */

        // Modeling electrons:
        // a[i][1] = a[i][1] * (1 - pow(M_E, (-1 * current_z_position / sp.particle_velocity)/mu_lifetime)); // why is mu_lifetime included for an electron modeling array?

        // Modeling muonium: 
        a[i] = a[i] * (pow(M_E, (-1 * current_z_position / sp.particle_velocity)/mu_lifetime));
        
        if (a[i]<min){ // could also have min = 0
            a[i]=0;
        }

        // a[i][1] = a[i][1]/max;//divides each element of the array by max, normalize it
        if (logchoice == 1 && a[i] > 0)
        {
            a[i] = log(a[i]/min); // it tries to make sure all elements are above 1, then takes the logarithm to scale the stuff better
        }
    }
}

int find_element_position_in_array(int value, int  * arr)
{
    for (int i = 0; i < sp.number_of_rows_fourier_coefficient_array; i++)
    {
        if (value ==  * (arr + i))
        {
            return(i);            
        }
    }
    printf("Error! n or m value does not match with any value of x2pnt array\n");
}
