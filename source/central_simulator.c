/* 
 * Central Simulator
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
 * 				*
 *
 * This central program takes user input for:
 *     inclusion or exclusion of gravitational acceleration,
 *     resolution of output graph,
 *     velocity of particles,
 *     pitch of gratings,
 *     choice of intensity profile or final simulation paths,
 *     and choice of logarithmic or linear intensity scale for the plot (only if the full simulation is requested),
 * and calculates the necessary parameters in order to perform the simulation with the user's desired arguments.
 *
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cfloat>
#include <ctype.h>
#include <complex.h>
#include <string> 
#include <complex.h>
// Fix for islimit error.
#include <limits.h> 		

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
// This program uses Root's TApplication.h and creates a plot based on values passed into the function.
#include "SimplePlot.h" 
// Contains functions of other miscellaneous tasks, such as taking care of boundary conditions, checking the max value in an array, etc.
#include "Misc.h"
// This program contains the functions defining the Gaussian Schell-model beam (GSM beam) behavior.
#include "BeamParams.h"      
// Contains functions that compute the intensity profiles before the first grating, after the 1st, and the second grating.
#include "Gratings.h"
// Contains functions that compute various phase shifts due to different physical effects
#include "PhaseShifts.h"

// sp is the simulation parameters structure that contains all of the simulation dependent variables. The struct is located in Misc.h.
simparam sp;

/*
 * The function 'main' is the main program that contains all of the global data and function calls to other functions
 * in order to simulate wave particles through a diffraction grating. It takes an integer and an array as arguments
 * and outputs a ROOT graph display.
 *
 * The main simulation parameters are initialized by argument values. This means that the user has to input
 * these arguments on the command line when executing the program. An example would be "./program 1 300 6400 100 1 1".
 *
 * These are the simulation parameters that the user has to provide, in the format "<position>. <explanation>":
 * 1. Account for gravity? 1 = True, 0 = False.
 * 2. Account for Van der Waals effects? 1 = True, 0 = False.
 * 3. Resolution in pixels [300-400 recommended]. 
 * 4. Velocity of particles in m/s.
 * 5. Pitch of gratings in nm.
 * 6. Output the total simulation [1]? or the final interference pattern [2]?
 * 7. If previous argument is 1, use a logarithmic scale for plotting intensities [1] or a normal scale [0]?
 *
 * Remember that argv[] is the standard argument vector used in C for inputting command-line arguments to the program.
 * The "simparam" structure is defined in Misc.h.
 */
int main(int argc, char *argv[])
{
	// Account for gravity?	1 = True, 0 = False.
	sp.account_gravity		= atoi(argv[1]);
	// Account for Van der Waals? 1 = True, 0 = False.
	sp.account_van_der_waals	= atoi(argv[2]);
	// Resolution.
	sp.resolution			= atoi(argv[3]);
	// Velocity of beam particles.
	sp.particle_velocity		= atof(argv[4]);
	// The period is inputted in nanometers, but the program converts it to meters for internal usage.
	sp.grating_period		= atof(argv[5]) / 1.0e9;
	// Output the total simulation [1]? or the final interference pattern [2]?
	sp.simulation_option		= atoi(argv[6]);
	/*
	 * If the user asks for the total simulation (sp.simulation_option is 1), the program needs to know
	 * if the intensity scale of the plot is going to be linear or logarithmic (the latter helps
	 * to see where more of the particles go). sp.logchoice == 0 means linear scale, == 1 means
	 * logarithmic scale.
	 */
	if(sp.simulation_option == 1)
		sp.logchoice = atoi(argv[7]);
	else
		sp.logchoice = 0;
	
	// Gaussian Schell-model beam parameters.
	sp.initial_beamwidth				= 3.0e-5;	// In m.
	sp.initial_coherence_width			= 1.0e-6;	// In m.
	sp.initial_radius_of_wavefront_curvature	= -4.04;	// In m.
	
	/*
	 * de Broglie wavelength for the matter waves.
	 * de Broglie equation: wavelength = Planck constant / mass * velocity,
	 * where we approximate the mass as the sum of the muon and electron masses.
	 */
	sp.wavelength 			= 6.626068e-34 / (1.8926409e-28 * sp.particle_velocity);

	// Grating geometry parameters.
	sp.z_position_1st_grating	= 1.0e-6;		// In m.
	sp.z_position_2nd_grating	= 4.5e-2;		// In m.
	// Height of each slit is calculated as half the grating period (distance between gratings).
	sp.slit_height 			= sp.grating_period / 2;
	sp.grating1_open_fraction	= 0.4;			// TODO LAcomment: explain.
	sp.grating2_open_fraction	= 0.4;			// TODO LAcomment: explain.
	sp.G2_x				= 5e-8;			// In m. TODO LAcomment: explain based on Dr McMorran's info.
	sp.theta			= 1e-6;			// In degrees. This is the relative rotation angle between the two gratings.
	sp.grating_thickness		= 1.0e-6;		// In m. Used for the VdW effect for atoms: see PhaseShifts.c.
	sp.wedge_angle			= 0* M_PI/180;		  	// In degrees. Grating wedge angle.
	sp.tilt_angle			= 0* M_PI/180; 	      		// Tilt.

	/* 
	 * Spatial parameters.
	 *
	 * These parameters control 'the size of the grating' in the x direction:
	 * delta x = (x_end - x_start) should equal the height of one grating.
	 *
	 * The total distance the beam travels is delta z = (z_end - z_start),
	 * such that the distance between gratings is controlled by both the
	 * z_position_1st_grating and z_position_2nd_grating variables above,
	 * and z_start and z_end below.
	 */	
	sp.z_start    			= -1.0e-2;	      	// In m. z start position
	sp.z_end      			= 1.0e-1;	      	// In m. z end position
 	sp.x_start    			= -5.0e-3;     		// In m. x start position
	sp.x_end      			= 5.0e-3;      		// In m. x end position
	/* 
	 * TODO Ycomment: y_start and y_end are not being used in the code. 
	 * We didn't delete it for now in case the value they assume help us understand something later
	 */
	//sp.y_start  = -5.0e-3;      				// y start position
	//sp.y_end    = 5.0e-3;       				// y end position

	sp.intensity_cutoff				= 1e-10;	// Point at which the intensity cuts off and is treated as 0.
	sp.number_of_rows_fourier_coefficient_array	= 41;  		// Rows of real_part_fourier_coefficient_array and imaginary_part_fourier_coefficient_array arrays; used to calculate phase shift.

	/* 
	 * The program prints a standard message before proceeding to the simulation. It includes a copyright notice
	 * and displays the values the user has given as input for the simulation.
	 */
	printf("Interferometer Simulation 1.0\n");
	printf("Copyright (C) 2016 Antimatter Gravity Interferometer Group, Illinois Institute of Technology (IIT).\n");
	printf("License: GNU GPL version 2\n----------\n");
	//printf("Interferometer Simulation comes with ABSOLUTELY NO WARRANTY; for details type 'show warranty'.\n");
	//printf("This is free software, and you are welcome to redistribute it under certain conditions; type 'show copying' for details.\n\n");
	printf("The simulation will run as follows.\n");
	printf("Is gravity being considered? ");
	if (sp.account_gravity == 1)
		printf("yes\n");
	else
		printf("no\n");
	printf("Is Van der Waals interaction being considered? ");
	if (sp.account_van_der_waals == 1)
		printf("yes\n");
	else
		printf("no\n");
	printf("Resolution of simulation plot: %3.0f pixels\n", sp.resolution);
	printf("Velocity of particles: %4.1f m/s\n", sp.particle_velocity);
	printf("'Period' of gratings (distance between two successive slits): %3.1f nm\n", sp.grating_period * 1.0e9);
	printf("Computing the full simulation or just the final interference pattern (using relative intensities)? ");
	if (sp.simulation_option == 1) {
		printf("full simulation\n");
		printf("Is the intensity being plotted in a logarthmic scale? ");
		if (sp.logchoice == 1)
			printf("yes\n");
		else
			printf("no\n");
	}
	else
		printf("final interference pattern\n");
	printf("Press 'Enter' to continue.");
	while (getchar() != '\n')
		;

	/*
	 * The intensity_array will be used in the calculation of the intensity of each row of pixels along the z direction.
	 * The x_positions_array contains the positions (in meters) of the pixels along the x direction and does not change as the code runs,
	 * given that the beam propagates in the z direction. That is, with each step in z the intensity_array holds the intensity
	 * values associated with each position in x.
	 */
	double *intensity_array;						// Intensity array.
    	double *x_positions_array;						// Array of x position of intensity.
	intensity_array = (double*) calloc(sp.resolution, sizeof(double)); 	// Allocate dynamic memory for intensity array.
	x_positions_array = (double*) calloc(sp.resolution, sizeof(double)); 	// Allocate dynamic memory for horizontal position array.
	int initial_z_position;							// Where z position begins.
	double max;								// Stores computed max value of intensity at a specific x location.
   	int total_number_of_pixels = sp.resolution * sp.resolution;  		// Pixels on full simulation graph.
    	double *pixel_array_memory = (double*) calloc(total_number_of_pixels, sizeof(double)); 	// Allocating dynamic memory for pixel array.
	double z_resolution = (sp.z_end-sp.z_start)/sp.resolution; 		// Step resolution used in computation.

     
	/*
	 * The following functions are used to calculate Gaussian Schell-model (GSM) values at the first grating:
	 * 	1. 'calculate_width' computes and outputs either the Gaussian Schell-model (GSM) beam width (in meters),
	 *   	or the GSM beam coherence width. Both quantities evolve according to the same formula, so the last parameter
	 * 	with which this function is called decides which quantity is being calculated. See the function's definition
	 *	in BeamParams.c and the program documentation for more information.
	 *
	 * 	2. 'v' computes and outputs the GSM radius of wavefront curvature.
	 *
	 * The function 'v' has the following parameters :
	 * 	location of first grating [set to 1 micron],
	 *	initial radius of wavefront curvature,
	 *	initial coherence width,
	 *	initial beam width.
	 *
	 * The function 'calculate_width' is similar, but it has to be called with one extra parameter at the end of the parameter list
	 * in order to calculate the correct width. If the last parameter is the initial beamwidth, then the function will calculate that
	 * beam parameter's evolution; if it's the initial coherence width, that is the beam parameter being calculated.
	 *
	 * The functions 'calculate_width' and 'v' output a double.
	 */
	double w1  = calculate_width(sp.z_position_1st_grating, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_beamwidth);
	double el1 = calculate_width(sp.z_position_1st_grating, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_coherence_width);
	double r1  = calculate_wavefront_radius(sp.z_position_1st_grating, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth);

	// Developer: follow this indent structure.
	if (sp.simulation_option == 1) {
		initial_z_position = 0;
	}
	else if (sp.simulation_option == 2) {
		initial_z_position = sp.resolution - 1;
	}

	for (int i=0; i<sp.resolution; i++)
		x_positions_array[i] = sp.x_start + (i) * ((sp.x_end-sp.x_start)/(sp.resolution-1));
	
	for (int i=(initial_z_position); i<sp.resolution; i++) {
		// Each time the loop repeats, you reset the array's positions and intensities to zero. 
		memset(intensity_array, 0, sp.resolution * sizeof(double));
		
		// Where you are with respect to z.
		double current_z_position = sp.z_start  +  i * z_resolution;		 

		/*
		 * These control structures determine where you are on zloc: depending on that, you interact with different gratings.
		 * The functions 'gpM' compute the intensity profile after the M-th grating. Thus:
		 * 	'intensity_after_2nd_grating' computes the intensity profile after the beam hits the second grating; 
		 * 	'intensity_after_1st_grating' computes the intensity profile after the beam hits the first grating;
		 * 	'get_initial_intensity' computes the intensity profile before the beam hits any grating.
		 * Note that the intensity profile is an array of positions and their respective intensities. 
		 * Also, keep in mind these functions modify the intensity and x positions arrays, but don't return any values.
		 * 
		 * Their arguments are:
		 * 	For 'intensity_after_2nd_grating':
		 *		 z position after gratings,
		 *		 beam coherence width,
		 *		 beam width,
		 *		 radius of curvature,
		 *		 X positions,
		 *		 x intensity profile.
		 * 	For 'intensity_after_1st_grating':
		 * 		 z position after gratings,
		 * 		 beam coherence width,
		 *		 beam width,
		 *		 radius of curvature,
		 *		 X positions,
		 *		 x intensity profile.
		 * 	For 'get_initial_intensity':
		 * 		 z position,
		 *		 X positions,
		 *		 x intensity profile.
		 */
		if (current_z_position > sp.z_position_2nd_grating) { 
			printf("Entering intensity_after_2nd_grating for row z = %d\n",i); //checking if the looping is working
			// If the location is above z_position_2nd_grating [which is currently 1]:
			intensity_after_2nd_grating(current_z_position, el1, w1, r1, x_positions_array, intensity_array); 
			/* 
			 * The function 'maximumvalue' outputs the maximum value in a given array.
			 * Its arguments are an array and the length of that array [integer].
			 * Here it gives the largest intensity value.
		     	 */
			max = maximumvalue(intensity_array, sp.resolution);
			/* 
			 * The function 'ixgenerator' is the x direction intensity calculator. It normalizes and compares intensities
			 * to cutoff value, then determines which value to input to the array of x intensitites. Its arguments are:
			 * 	intensity at the time of calculation (an initially empty array);
			 * 	zlocation;
			 * 	choice of scale for the plot.
			 * The function only modifies the array fed to it; it doesn't return any value. 
			 */  
		    	ixgenerator(intensity_array, current_z_position, sp.logchoice); 
		}
		else if (current_z_position > sp.z_position_1st_grating) {
			// If interacting with the first grating, calculates intensity profile.
			printf("Entering intensity_after_1st_grating for row z = %d\n",i); //checking if the looping is working
			intensity_after_1st_grating(current_z_position, el1, w1, r1, x_positions_array, intensity_array); 
			// Max value of intensity calculated here.
			max = maximumvalue(intensity_array, sp.resolution); 
			// As before.		
			ixgenerator(intensity_array, current_z_position, sp.logchoice); 
		}
		else {
			// Simple GSM propagation until it hits the first grating.
			printf("Entering get_initial_intensity for row z = %d\n",i); //checking if the looping is working
			get_initial_intensity(current_z_position,x_positions_array, intensity_array);
		
			// If at the origin?
			max = maximumvalue(intensity_array, sp.resolution); 
			// As before.
			ixgenerator(intensity_array, current_z_position, sp.logchoice); 
		}   

		for (int j=0; j<sp.resolution; j++ ) {
			// TODO LAcomment: update? "Resolution * i + j; still keeping track of location."
			int f = sp.resolution * i + j; 
		    	// The f-th element of pixel_array_memory is set to be the intensity of the beam at the j-th point.	
		    	pixel_array_memory[f] = intensity_array[j]; 
		}
	}
	
	if (sp.simulation_option == 1) {
		// Free the memory used by this array; since the simulation is over, pixel_array_memory has all the data.
		free(x_positions_array); 
		// Same as above.
		free(intensity_array); 
		// Using ROOT to plot pixel_array_memory.
		SimplePlot::twoD("Intensity graph as particles diffract through gratings",pixel_array_memory,sp.z_start,sp.z_end,sp.x_start,sp.x_end,sp.resolution,sp.resolution); 
	}
	else if (sp.simulation_option == 2) {
		// Using ROOT to plot intensity vs. position at end of interferometer.
		SimplePlot::graph("Relative Intensity Along Final Grating", x_positions_array, intensity_array, sp.resolution);		
		// Free the memory used by this array, since the simulation is over.
		free(x_positions_array); 
		// Same as above.
		free(intensity_array); 
		  
	}

	// Free up the space used by the pixel_array_memory array.
	free(pixel_array_memory);
}
