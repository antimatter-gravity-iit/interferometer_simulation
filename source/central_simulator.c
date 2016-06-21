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
 *  Collaborators:
 *     Arthur Romero, 16 July to 24 August 2015
 *         -- Original author
 *     Adam Denchfield, 20 October 2015 to December 2015
 *         -- Additional parameters and functionality
 *     Melanie Cornelius (Dooley), 1 Nov 2015 to present
 *         -- Optimizations, standards, and readability edits
 *         -- Comments tagged with Mcomment
 *     Isaac Gewarges, January 2016 to April 28 2016
 *         -- Optimizations, standards, readability edits, and variable organization.
 *     Lucas Maia Rios, 23 May 2016 to present
 *         -- Comments tagged with LRcomment
 *     Lucas Neves Abrantes, 23 May 2016 to present
 *         -- Comments tagged with LAcomment
 *     Yuri Rossi Tonin, 6 May 2016 to present
 *         -- Comments tagged with Ycomment
 *
 * 				*
 *
 * This central program takes user input for:
 *     electronic or atomic beam type,
 *     inclusion or exclusion of gravitational acceleration,
 *     resolution of output graph,
 *     velocity of particles,
 *     pitch of gratings,
 *     and choice of intensity profile or final simulation paths,
 * and calculates simulation parameters in order to perform the simulation with the user's desired arguments.
 *
 * Top-down view of gratings:
 *
 *
 *                                    ---------------------   -> at G2_z = 1,         z2
 *
 *                                    ---------------------   -> at G1_z = 0.000001,  z1
 *
 *                                       (beam going up)
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

// Rows of ReT and ImT arrays; used to calculate phase shift.
int rowsT = 41;                              

/*
 * The function 'main' is the main program that contains all of the global data and function calls to other functions
 * in order to simulate wave particles through a diffraction grating. It takes an integer and an array as arguments
 * and outputs a ROOT graph display.
 *
 * TODO LAcomment: Understand why this is a char array when all the elements are integers.
 * argv[x] is the argument array [vector] used for changing the inputs to the program.
 *
 * List of elements of argv[x], in the format "<position>. <explanation>": 
 * 1. Account for gravity? 1 = True, 0 = False.
 * 2. Electron beam or atom beam? Electron beam = 1, Atom beam = 2. ELECTRON BEAM NOT MODELED!!!
 * 3. Resolution [300-400 recommended]. 
 * 4. Velocity of particles in m/s.
 * 5. Pitch of gratings [in nm] 
 * 6. Output the total simulation [1]? or the final interference pattern [2]?
 * 7. (if argv[6] == 1), logscale [1] or normal scale [0]?
 */
 
int main(int argc, char *argv[]){ 

	// TODO LAcomment: rename variables and if necessary add brief comments. Verify if existing comments are needed.
	// Initializing all simulator parameters by either default values or argument values. Variable definitions and comments are in misc.h.
	sp.accountGrav = atoi(argv[1]);
	sp.elecOrAtom = atoi(argv[2]);
	// Velocity of particle.
	sp.vel = atoi(argv[4]);
	// Why was it defined like this? This equation seems to come from DeBroglie's model.
	sp.energy = 1.5e-18 / pow(1e-11,2) * (1);
	// sp.energy = 1.5e4
	sp.simchoice = atoi(argv[6]);
	// TODO LAcomment: add comment explaining why this 'if' exists.
	if(sp.simchoice == 1)
	sp.logchoice = atoi(argv[7]);
	sp.useimagecharge = 0;
	sp.eta1 = 0.4;
	sp.eta2 = 0.4;
	sp.g_period = atoi(argv[5]);
	sp.initial_radius_of_wavefront_curvature = -4.04;
	sp.initial_coherence_width = 1.0e-6;
        sp.initial_beamwidth = 3.0e-5;
	sp.G1_z = 1.0e-6;
	sp.G2_z = 1.0;
	sp.G2_x = 5e-8;
	sp.theta = 1e-6;
	sp.thick = 1.4e-8;
	sp.Gthick = 1.0e3;
	// Wedge angle.
	sp.wedgeangle =	0;
	// Tilt.
	sp.tilt = 0; 
	// Resolution.
	sp.res = atoi(argv[3]);
	sp.zstart = -0.1;
	sp.zend = 2.1;
	sp.xstart = -2.0e-4;
	// xend.
	sp.xend = 2.0e-4;
	// ystart.
	sp.ystart = -1.1e-4;
	// yend.
	sp.yend = 1.1e-4;
	sp.height = (sp.g_period / 2) / 1.0e9;
	// Point at which the intensity cuts off and is treated as 0. Can also be 5e-5 like in McMorran thesis, or 0.001.
	sp.cutoff = 1e-6;
	// Scaled logarithmically so they can see where more of the particles go [0] = no, [1] = yes.
	//sp.logchoice = 0;
	
	/*
	 * TODO LAcomment: make sense of this. Remember Mcomment: "I tried to fix this, but what is it actually saying?
	 * It's the intensities of the x-positions and the intensities?  That doesn't make sense to me."
	 * Original description: "Initializing two arrays to contain the intensities and xpositions of each intensity."	
	 */
	
	double *Grat3I;							// Intensity array.
    	double *Grat3x;							// Array of x position of intensity.
	Grat3I = (double*) calloc(sp.res, sizeof(double)); 		// Allocate dynamic memory for intensity array.
	Grat3x = (double*) calloc(sp.res, sizeof(double)); 		// Allocate dynamic memory for horizontal position array.
	int zlocstart;							// Where z position begins.
	int rows = sp.res;						// Numbers of horizontal component arrays of full simulation graph.
	double max;							// Stores computed max value of intensity at a specific x location.
		
    	// TODO LAcomment: make sense of these and rename variables accordingly.

   	int izxnumels = rows * rows;  					// Pixels on full simulation graph.
    	double izxsize = izxnumels * sizeof(double); 			// Size of pixels array.
    	double *izx = (double*) calloc(izxnumels, sizeof(double)); 	// Allocating dynamic memory for pixel array.
    	double zres = (sp.zend-sp.zstart)/sp.res; 			// Step resolution used in computation.
     
	/*
	 * The following functions are used to calculate Gaussian Schell-model (GSM) values at the first grating:
 	 * 	'w' computes and outputs the Gaussian Schell-model (GSM) beam width in meters [m];
         * 	'v' computes and outputs the GSM radius of wavefront curvature; and
	 * 	'el' computes and outputs the GSM beam coherence width.	
	 *
	 * All of them take the same arguments:
	 * 	height of first grating [set to 1 micron],
	 *	initial radius of wavefront curvature,
	 *	initial coherence width,
	 *	initial beam width.
	 *
	 *  Note that they all output a double.
	 */
		
	double w1 = calculate_width(sp.G1_z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_beamwidth);
	double el1 = calculate_width(sp.G1_z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_coherence_width);
	double r1 = v(sp.G1_z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth);

	// Follow this indent structure in the future.
	if (sp.simchoice == 1) {
		sp.logchoice = atoi(argv[7]);
		zlocstart = 0;
	}
	else if (sp.simchoice == 2) {
		zlocstart = sp.res - 1;
	}

	for (int i=(zlocstart); i<rows; i++) {
		// TODO: LAcomment: said "i=299 is just to get last row of z." What?
		memset(Grat3x, 0, rows * sizeof(double)); 
		// Each time the loop repeats, you reset the array's positions and intensities to zero. 
		memset(Grat3I, 0, rows * sizeof(double));

	     	// Where you are with respect to z.
		double zloc = sp.zstart  +  i * zres;		 

		/*
		 * These control structures determine where you are on zloc: depending on that, you interact with different gratings.
		 * The functions 'gpM' compute the intensity profile after the M-th grating. Thus:
		 * 	'gp2' computes the intensity profile after the beam hits the second grating; 
		 * 	'gp1' computes the intensity profile after the beam hits the first grating;
		 * 	'gp0' computes the intensity profile before the beam hits any grating.
		 * Note that the intensity profile is an array of positions and their respective intensities. 
		 * Also, keep in mind these functions modify the Grat3I and Grat3x arrays, but don't return any values.
		 * 
		 * Their arguments are:
		 * 	For 'gp2':
		 *		 z position after gratings,
		 *		 beam coherence width,
		 *		 beam width,
		 *		 radius of curvature,
		 *		 X positions,
		 *		 x intensity profile.
		 * 	For 'gp1':
		 * 		 z position after gratings,
		 *		 radius of curvature,
		 * 		 beam coherence width,
		 *		 beam width,
		 *		 X positions,
		 *		 x intensity profile.
		 * 	For 'gp0':
		 * 		 z position,
		 *		 X positions,
		 *		 x intensity profile.
		 */
		if (zloc > sp.G2_z) { 
			//printf("Entering gp2 for row equal to %d\n",i); //checking if the looping is working
			// If the location is above G2_z [which is currently 1]:
			gp2(zloc, el1, w1, r1, Grat3x, Grat3I); 
			/* 
			 * The function 'maximumvalue' outputs the maximum value in a given array.
			 * Its arguments are an array and the length of that array [integer].
			 * Here it gives the largest intensity value.
		     	 */
			max = maximumvalue(Grat3I, rows);
			/* 
			 * The function 'ixgenerator' is the x direction intensity calculator. It normalizes and compares intensities
			 * to cutoff value, then determines which value to input to the array of x intensitites. Its arguments are:
			 * 	intensity at the time of calculation (an initially empty array);
			 * 	zlocation;
			 * 	choice of scale for the plot; and
			 * 	number of rows.
			 * The function only modifies the array fed to it; it doesn't return any value. 
			 */  
		    	ixgenerator(Grat3I, zloc, sp.logchoice, rows); 
		}
		else if (zloc > sp.G1_z) {
			// If interacting with the first grating, calculates intensity profile.
			//printf("Entering gp1 for row equal to %d\n",i); //checking if the looping is working
			gp1(zloc, r1, el1, w1, Grat3x, Grat3I); 
			// Max value of intensity calculated here.
			max = maximumvalue(Grat3I, rows); 
			// As before.		
			ixgenerator(Grat3I, zloc, sp.logchoice, rows); 
		}
		else {
			// Simple GSM propagation until it hits the first grating.
			//printf("Entering gp0 for row equal to %d\n",i); //checking if the looping is working
			gp0(zloc,Grat3x, Grat3I);
		
			// If at the origin?
			max = maximumvalue(Grat3I, rows); 
			// As before.
			ixgenerator(Grat3I, zloc, sp.logchoice, rows); 
		}   

		for (int j=0; j<rows; j++ ) {
			// TODO LAcomment: update? "Resolution * i + j; still keeping track of location."
			int f = rows * i + j; 
		    	// The f-th element of izx is set to be the intensity of the beam at the j-th point.	
		    	izx[f] = Grat3I[j]; 
		}
		}
	
	if (sp.simchoice == 1) {
		// Free the memory used by this array; since the simulation is over, izx has all the data.
		free(Grat3x); 
		// Same as above.
		free(Grat3I); 
		// Using ROOT to plot izx.
		SimplePlot::twoD("Intensity graph as particles diffract through gratings",izx,-200,200,0.0,220,rows,rows); 
	}
	else if (sp.simchoice == 2) {
		// Using ROOT to plot intensity vs. position at end of interferometer.
		SimplePlot::graph("Relative Intensity Along Final Grating", Grat3x, Grat3I, rows);  
		// Free the memory used by this array, since the simulation is over.
		free(Grat3x); 
		// Same as above.
		free(Grat3I); 
	}

	// Free up the space used by the izx array.
	free(izx);
}
