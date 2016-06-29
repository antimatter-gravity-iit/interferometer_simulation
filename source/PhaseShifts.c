// PhaseShifts.c - contains the functions responsible for implementing the phase shifts of different physical effects

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>


#include "PhaseShifts.h"
#include "BeamParams.h"
#include "Misc.h"
#include "Gratings.h"

#include <complex.h>

double ( * ReTandImTgenerator(double ReTorImTar[], int ReTorImT, double current_z_position))
{
	double gravity_acceleration = -9.8;    // acceleration due to gravity. 
  	double tilt = sp.tilt;  //0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
  	double eta1 = sp.eta1; //.4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it. Varname could be changed to better represent it.
  	double eta2 = sp.eta2; //.4; //G2 open fraction; how open the second grating is.
    
	//values not included in simparam structure.  can be moved there but not entirely necessary 
	double chargeratio =0.0; //strength of image charge (units of e, electron charge); values of 0.03, 0.05, or more can be had //not used
    	double C3 = 2.0453e-2; 			// the VdW coefficient for hydrogen (assumed to be the same for muonium). In meV * nm^3.
    	double hbar = 6.58212e-13; 		// Planck's reduced constant in meV * s.
  
  	double eta = sp.slit_height/sp.grating_period; // ratio of 'height' of slit/windows in gratings to the period of the gratings
    	double alpha = sp.wedgeangle * M_PI/180; // depends on wedgeangle above, which is a relatively free parameter. Appears to be bend of 'window' (slits in grating), if they bend forward or not.
    	double beta = tilt * M_PI; // depends on tilt angle, = 0 if beam is normal to gratings
    	double exnmleft; // how many nm from the left side of each slit are we?
    	double exnmright; // how many nm from the right side of each slit are we?
    	long double xmin; // beginning of path of wave through the slit
    	long double xmax; // end of path of wave/beam through the slit
    	double fc;
    	long double phase_van_der_waals; // phase shift if dealing with neutral atoms/molecules
    	double phase_gravity; // phase shift due to gravity.
    	double ex;
    	double time_free_fall;
    	int j;

	// if the beam is not normal/perpendicular to the gratings it encounters
	if (beta>=0) {
		/*
		 * TODO: LAcomment and LRcomment: make sense of this and explain.
 		 * Previous comment: "minimum distance beam travels through slit; or maybe it's when the 1st order
		 * diffracted beams are going in diagonally, what is the min x?
		 */
		xmin= sp.slit_height * (1/sp.resolution - cos(beta)/2); 
		// if the beam is very orthogonal to gratings (almost 90 degrees), or wedge angle is significant
		if (beta<=alpha) {
			xmax=(sp.slit_height * cos(beta))/2-sp.slit_height/sp.resolution;
		}
		// if beam is not very perpendicular to gratings, then it travels through the slit diagonally, covering more distance, more image charge interaction, etc.	
		else {
		xmax= sp.slit_height  *  cos(beta)/2 - sp.slit_height/sp.resolution  +  sp.thick  *  (tan(alpha)-tan(beta));
		}
    	}
	// if beta < 0; this time xmin changes, xmax is the same
	else {
		xmax = (sp.slit_height * cos(beta)/2)-sp.slit_height/sp.resolution;
		/*
		 * fabsl is for long doubles and returns a long double absolute value; once again,
		 * if the tilt isn't that bad, one bound (this time xmin) is just sp.slit_height * cos(beta)/2  +  sp.slit_height/res.
		 */
		if (fabsl(beta)<=alpha) {
		  xmin = -((sp.slit_height * cos(beta))/2) + sp.slit_height/sp.resolution; 
		}
		else { // if the beam is far from perpendicular to grating slits
		    xmin = -((sp.slit_height * cos(beta))/2) + sp.slit_height/sp.resolution - sp.thick * (tan(alpha)-tan(beta));
		}
	}  
    	for (int n=-((sp.rowsT-1)/2);n<=((sp.rowsT-1)/2);n++) {
       		/*
		 * TODO: LAcomment: make sense of this and explain. 
		 * Original comment stated the code below was somewhat copied from McMorran's 'NOT IMPLEMENTED YET' section.
		 * "resolution = step resolution in x-axis. ex += height of window / steps = (40 nm / 1000)"
		 */
		for (ex=xmin; ex<xmax; ex +=sp.slit_height/sp.resolution) {
			// ex is how far you are from the grating 'wall'
			// exnm is how far from the wall in nanometers.
			exnmleft = ex * 1.0e9;
			exnmright = (xmax - ex) * 1.0e9; 
			
			/* TODO LAcomment : rewrite. 
			 * "fc is another electron thing; or the diffraction pattern? 2pi*n*x/period?
			 * so the first fc = 2  *  pi  *  -20  *  xmin / period, last fc = 2  *  pi  *  20  *  xmax / period.
			 * Looks like fc is a quantity proportional to distance from bottom/top of each 'window'/slit"
			 */
			fc = 2 * M_PI * n * ex/sp.grating_period; 
			
			// TODO LAcomment: rewrite. "How much time has passed for a point in this system?" 
			time_free_fall = current_z_position / sp.particle_velocity;

			/* Atoms will fall due to gravity.
			 * According to Dr. Daniel Kaplan's paper at arxiv.org/ftp/arxiv/papers/1308/1308.0878.pdf,
			 * the phase shift caused is 2 * pi * g * t^2 / d, where t is the time in free fall and d is the period of the gratings.
			 */
			if (sp.account_gravity == 1)
				// phase shift due to gravity on particles
				phase_gravity = (2 * M_PI * gravity_acceleration * pow(time_free_fall, 2)) / sp.grating_period;
			else
				phase_gravity = 0;
			  
			if (sp.account_van_der_waals == 1) {
				// phase_van_der_waals is phase shift on Muonium/other neutral molecules due to Van der Waals effects through the gratings.
				if (exnmleft == 0 || exnmright == 0)
					phase_van_der_waals = 0;
				else
					phase_van_der_waals = -C3 * sp.grating_thickness / (hbar * sp.particle_velocity * pow(exnmleft, 3)) -  -C3 * sp.grating_thickness / (hbar * sp.particle_velocity * pow(exnmright, 3));
			}
			else {
				phase_van_der_waals = 0;
			}	

			// j goes from 0 to 40 (right now sp.rowsT = 41)
			j = n + ((sp.rowsT-1)/2);
				
			// if it's the ReT array
			if (ReTorImT == 1)										
				// so fc and phase_van_der_waals are both phase shifts; angles.
				ReTorImTar[j]  += cos(phase_van_der_waals + fc + phase_gravity);
			// if it's the ImT array
			else if (ReTorImT == 2)
				// so fc and phase_van_der_waals are both phase shifts; angles.
				ReTorImTar[j]  += sin(phase_van_der_waals + fc + phase_gravity); 
		}
	}
	
	for (int i=0; i<sp.rowsT; i++) {	
		// TODO LAcomment: explain. "So is this some sort of normalization?"	
		ReTorImTar[i] = ReTorImTar[i]/sp.resolution;
	}
}
