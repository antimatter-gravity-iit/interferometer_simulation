// Gratings.c -> Holds the functions that determine the diffraction after each of the gratings in the system

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
	double jj;

	// Width of incoming beam, calculated from initial parameters.
	w1 = calculate_width(z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_beamwidth); 

	for(int i=0; i<sp.resolution; i++) {
		jj = pow((x_positions_array[i]/w1),2); 		// jj = (xpos/beamwidth)^2 
		intensity_array[i]=exp(-(M_PI * jj)); 		// a[i][1] is the intensity of the beam at the xposition at step i.
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

	double z12    = current_z_position - sp.z_position_1st_grating;	//z location between 1st and 2nd gratings
	/* 
	 * A free parameter. Beta variable below depends on this. If beam is perpendicular to gratings, then tilt (and thus Beta) is 0.
	 * This is the twist about the x-axis.
	 */

    	double coef; 
    	double lim    = 5;
    	double tilt   = sp.tilt; 
	double eta1   = sp.eta1; 				// G1 open fraction; how open the first grating is. With 0.4 open, a little over than half the muonium should pass through.
    	double eta2   = sp.eta2; 				// G2 open fraction; how open the second grating is.		 
	double eta    = sp.slit_height/sp.grating_period; 	// ratio of window 'height' to period of grating
    	double alpha  = sp.wedgeangle * M_PI/180; 		// alpha and beta have been defined in almost every other function. Global variables? 
    	double beta   = tilt * M_PI; 				// defined in other functions too, same purpose.	
	/*
	 * Explanation of variables:
	 * w2 = GSM width of beam after the first grating.
	 * r2 = radius of GSM wavefront curvature after the first grating
	 * el2 = GSM beam coherence width after the first grating.
    	 */
	double w2  = calculate_width(z12, r1, el1, w1, w1);		// width of beam between z1 and z2 (after grating 1)
	double r2  = calculate_wavefront_radius(z12, r1, el1, w1);	// radius of wavefront curvature between grating 1 and 2
	double el2 = calculate_width(z12, r1, el1, w1, el1);		// beam coherence width
    	
	int pos[41]={0};

	for (int i=0; i<sp.rowsT; i++) // since sp.rowsT is currently 41, pos[i] = -20 to 20.
	{
	pos[i]=i-((sp.rowsT-1)/2);
	}

	double ReT[41]={0};
	int RealorIm = 1;
	ReTandImTgenerator(ReT, RealorIm, current_z_position); // calculates phase shift
		
	double ImT[41]={0};
	RealorIm = 2;
	ReTandImTgenerator(ImT, RealorIm, current_z_position); // calculates phase shift for imaginary part

	for (int i=0; i<sp.resolution; i++) { 
		for (int n=-lim; n<=lim; n++) {
			for (int m=-lim; m<=lim; m++) {
				double dn =n-m; 
				// TODO: explain what n, m, dm, dn are. LR, Y 
				double dm = (m + n)/2;

				coef = ReT[x2pnts(n, (int * )pos)] * ReT[x2pnts(m,(int * )pos)] + ImT[x2pnts(n,(int * )pos)] * ImT[x2pnts(m,(int * )pos)];
				
				coef = coef * exp(-M_PI * pow((dn * sp.wavelength * z12)/(sp.grating_period * el2),2));
				// added isfinite macro in order to avoid inf values

				if (std::isfinite(coef)==0 || coef < sp.intensity_cutoff) { // if coef is infinite, then:
					coef=0;
				}
			      	else { // if coef ends up larger than cutoff value, add the values to the current a[i][1]'s intensities.
					intensity_array[i] +=   coef * exp(-M_PI * pow(((x_positions_array[i]-dm * sp.wavelength * z12/sp.grating_period)/w2),2)) * cos(2 * M_PI * (dn/sp.grating_period) * (x_positions_array[i]-dm * sp.wavelength * z12/sp.grating_period) * (1-z12/r2));
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
	double z12     = sp.z_position_2nd_grating - sp.z_position_1st_grating;
	double z23     = current_z_position - sp.z_position_2nd_grating;
	double mytheta = sp.theta;
	/* 
	 * A free parameter. Beta variable below depends on this. If beam is perpendicular to gratings, then tilt (and thus Beta) is 0.
	 * This is the twist about the x-axis.
	 */
    	double tilt   = sp.tilt; 
	double eta1   = sp.eta1;				// G1 open fraction; how open the first grating is. With 0.4 open, a little over than half the muonium should pass through. 
    	double eta2   = sp.eta2;				// G2 open fraction; how open the second grating is.
    	double eta    = sp.slit_height/sp.grating_period; 	// ratio of slit window 'height' to the period of the gratings
        double alpha  = sp.wedgeangle * M_PI/180;
    	double beta   = tilt * M_PI; 				// 0 if beam is normal to gratings
    	double theta  = M_PI * mytheta/180;
    	double d1     = sp.grating_period;			// period = period of gratings
    	double d2     = sp.grating_period;
    	double z13    = z12  +  z23; 				// z distance between grating 1 and 3
    	double phi    = 0;
    	double lim    = 5;
    	double resolution = sp.resolution; 			// This is the resolution we want this graph at.
    	double _Complex coef;

    	/* THIS FUNCTION IS USING GSM MODEL FROM MCMORRAN, CRONIN 2008
    	 * IT MUST BE CHANGED IF YOU WANT TO USE OLDER MODEL FROM BREZGER 2003
	 */

    	double dn = 0;
    	double dm = 0;
    	double m  = 0;
    	double n  = 0;
    	int    a5 = 0;
    	int    b5 = 0;
    	int    c5 = 0;
    	int    d5 = 0;
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
    	double el3x = calculate_width(z13, r1x, el1x, w1x, el1x);	// z13 == G2z - G1z + z_start + 0 * zres; GSM coherence width in x-axis
    	double w3x  = calculate_width(z13, r1x, el1x, w1x, w1x); 	// Beam width in x-axis
    	double v3x  = calculate_wavefront_radius(z13,r1x,el1x,w1x);	// Gaussian-Schell Model (GSM) radius of wavefront curvature in x-axis
    	double el3y = calculate_width(z13, r1y, el1y, w1y, el1y); 	// Coherence width in y-axis
    	double w3y  = calculate_width(z13, r1y, el1y, w1y, w1y); 	// Beam width in y-axis
    	double v3y  = calculate_wavefront_radius(z13,r1y,el1y,w1y);	// radius of wavefront curvature in y-axis
    
    	int pos[41]={0}; // array of 41 elements

    	for (int i=0; i<sp.rowsT; i++){     // sp.rowsT = rows of ReT and ImT arrays = 41 for now
    		pos[i]=i-((sp.rowsT-1)/2);  // so this goes from -20 to 20
    	}
    
    	double ReT[41]={0};     	// array of 41 0's for now.
    	int RealorIm = 1; 
    	ReTandImTgenerator(ReT, RealorIm, current_z_position); 
    	
	double ImT[41]={0};
    	RealorIm = 2;
    	ReTandImTgenerator(ImT, RealorIm, current_z_position); 
	

	for (int i=0; i<sp.resolution; i++) {
		for (int m1=-lim; m1<=lim; m1++) {
			for (int m2=-lim; m2<=lim; m2++) {
				for (int n1=-lim; n1<=lim; n1++) {
					for (int n2=-lim; n2<=lim; n2++) {
						
						dn = n1-n2;
						n  = ((double)(n1 + n2))/2;
						dm = m1-m2;
						m  = ((double)(m1 + m2))/2;
						a5 = ( x2pnts( m1, (int *)pos ) );
						b5 = ( x2pnts( m2, (int *)pos ) );
						c5 = ( x2pnts( n1, (int *)pos ) );
						d5 = ( x2pnts( n2, (int *)pos ) );

						//argument_d corresponds to the argument of equation 18b from McMorran & Cronin 2008. Note that y=0
                        			argument_d = -M_PI*(pow( x_positions_array[i]-sp.wavelength*z23*(n*cos(theta)/d2 + m*z13/(d1*z23) ),2 )/pow(w3x,2) + pow((n*sin(theta)*sp.wavelength)/(d2*w3y),2));
                        			//argument_f corresponds to the argument of equation 18c from McMorran & Cronin 2008. Note that y=0
                        			argument_f = -2*M_PI * x_positions_array[i] * ((dn*cos(theta)/d2)*(1-z23/v3x) + (dm/d1)*(1-z13/v3x));
                        			//argument_p corresponds to the argument of equation 18d from McMorran & Cronin 2008
                        			argument_p = (2*M_PI*sp.wavelength*z13*dm/d1)*(n*cos(theta)/d2 + m/d1)*(1-z13/v3x) +     (2*M_PI*sp.wavelength*z23*dn/d2)*( (m*cos(theta)/d1) * (1-z13/v3x) - (n*z23/d2)*(pow(cos(theta),2)/v3x)+pow(sin(theta),2)/v3y ) ;
                        			//argument_v corresponds to the argument of equation 18e from McMorran & Cronin 2008
                        			argument_v = -M_PI* pow(sp.wavelength*z23* (dn*cos(theta)/d2 + dm*z13/(d1*z23)),2)/pow(el3x,2) -M_PI*pow(dn*sin(theta)*sp.wavelength*z23/(d2*el3y),2);

						coef =        (ReT[a5] + ImT[a5] * _Complex_I); 
						coef = coef * (ReT[b5] - ImT[b5] * _Complex_I);						
						coef = coef * (ReT[c5] + ImT[c5] * _Complex_I);
						coef = coef * (ReT[d5] - ImT[d5] * _Complex_I);
						
						if (((__real__ coef) >= sp.intensity_cutoff) || ((__imag__ coef) >= sp.intensity_cutoff)) {
						    
                        				function_d   = exp(argument_d);                      
                        				function_v   = exp(argument_v);
                        				argument_f_p = argument_f + argument_p;
                        
                        				intensity_array[i] +=  ((__real__ coef) * cos(argument_f_p) - (__imag__ coef) * sin(argument_f_p))*function_d*function_v ;
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
