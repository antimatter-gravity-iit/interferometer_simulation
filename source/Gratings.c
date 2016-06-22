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

void ( * gp0(double z, double Grat3x[], double Grat3I[]))
// get intensity profile 
{
	time_t start, end; //starting a timer to get the time spent in function gp0
	start = clock();
	double diff=0;

	double xstart = sp.xstart;
	double xend = sp.xend;
	double w1;
	double jj;

	// Width of incoming beam, calculated from initial parameters.
	w1 = calculate_width(z, sp.initial_radius_of_wavefront_curvature, sp.initial_coherence_width, sp.initial_beamwidth, sp.initial_beamwidth); 

	for(int i=0; i<sp.res; i++) {
		Grat3x[i]= xstart + (i) * ((xend-xstart)/(sp.res-1)); 	// current x-position at step i is put into a[i][0]
		jj = pow((Grat3x[i]/w1),2); 				// jj = (xpos/beamwidth)^2 
		Grat3I[i]=exp(-(M_PI * jj)); 				// a[i][1] is the intensity of the beam at the xposition at step i.
	}

	end = clock();
	diff =((double) (end - start));
	//printf("%f\n",diff); //printing the time spent in gp0.
}


void ( * gp1(double zloc,double r1,double el1, double w1, double Grat3x[], double Grat3I[]))
// get intensity profile after one grating
{
	time_t start, end; //starting a timer to get the time spent in function gp1
	start = clock();
	double diff=0;

	double z12 = zloc - sp.G1_z;		//z location between 1st and 2nd gratings
	double energy = sp.energy;
	int rows = sp.res;
	int elecOrAtom = sp.elecOrAtom;
	double vel = sp.vel;
	int xpnts = rows;
	double width = rows;
	double abszloc = sp.height; 		//z position
	int accountGrav = sp.accountGrav;
	int rowsT =41;				// rows of ReT and ImT array
	double xstart = sp.xstart; 
    	double xend = sp.xend;
    	//double period = sp.g_period;
	double period = 0.000000100;		// period of grating - 100 nanometers.
	// Grating wedge angle. Variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
    	double wedgeangle = sp.wedgeangle;
    	int useimagecharge = sp.useimagecharge;	// whether or not to consider image charge effects. 0 for False.
	/* 
	 * A free parameter. Beta variable below depends on this. If beam is perpendicular to gratings, then tilt (and thus Beta) is 0.
	 * This is the twist about the x-axis.
	 */
    	double tilt =sp.tilt; 
    	double cutoff = sp.cutoff; 		// The point at which the intensity cuts off and is treated as 0.
	// G1 open fraction; how open the first grating is. With 0.4 open, a little over than half the muonium should pass through.
	double eta1 = sp.eta1; 
	// G2 open fraction; how open the second grating is.
    	double eta2 = sp.eta2; 		
    	double coef; 
    	double lim=5;
    	double lambda = sqrt((1.5 * pow(10,-18))/(energy)); 
	// wavelength of what particles/waves we're working with; 
	double eta = width/period; 		// ratio of window 'height' to period of grating
	//double vel = pow(2 * energy * e_charge/e_mass,1/2); electron velocity
    	double alpha = wedgeangle * M_PI/180; 	// alpha and beta have been defined in almost every other function. Global variables? 
    	double beta = tilt * M_PI; 		// defined in other functions too, same purpose.
    	/*
	 * Explanation of variables:
	 * w2 = GSM width of beam after the first grating.
	 * r2 = radius of GSM wavefront curvature after the first grating
	 * el2 = GSM beam coherence width after the first grating.
    	 */
	double w2=calculate_width(z12, r1, el1, w1, w1); // width of beam between z1 and z2 (after grating 1)
	double r2 = v(z12, r1, el1, w1); // radius of wavefront curvature between grating 1 and 2
	double el2 = calculate_width(z12, r1, el1, w1, el1); // beam coherence width
    	int pos[41]={0};

	for (int i=0; i<rowsT; i++) // since rowsT is currently 41, pos[i] = -20 to 20.
	{
	pos[i]=i-((rowsT-1)/2);
	}

	double ReT[41]={0};
	int RealorIm = 1;
	ReTandImTgenerator(ReT,energy, elecOrAtom, RealorIm, vel, width, abszloc, accountGrav); // calculates phase shift
	RealorIm = 2;
	double ImT[41]={0};
	ReTandImTgenerator(ImT,energy, elecOrAtom, RealorIm, vel, width, abszloc, accountGrav); // calculates phase shift for imaginary part

	for (int i=0; i<rows; i++) {
	// a[i][0] just comprises all the xpositions you are thinking of as the beam passes through a grating.
	Grat3x[i] = xstart + (i) * ((xend-xstart)/(xpnts-1));
	}
    
	for (int i=0; i<rows; i++) { 
		for (int n=-lim; n<=lim; n++) {
			for (int m=-lim; m<=lim; m++) {
				double dn =n-m; 
				// TODO: explain what n, m, dm, dn are. LR, Y 
				double dm = (m + n)/2;

				if (useimagecharge==0) { // if useimagecharge = 0, ignore image charge effects at G1. 
					coef = sinc(eta1 * M_PI * n)  *  (sinc(eta1 * M_PI * m) * pow((eta1), 2));
				}

				else { // if usechargeimage = 1, don't ignore image charge effects at G1.
					coef = ReT[x2pnts(n, (int * )pos)] * ReT[x2pnts(m,(int * )pos)] + ImT[x2pnts(n,(int * )pos)] * ImT[x2pnts(m,(int * )pos)];
				}

				// lambda is the wavelength of our particles/waves
				coef = coef * exp(-M_PI * pow((dn * lambda * z12)/(period * el2),2));
				// added isfinite macro in order to avoid inf values

				if (std::isfinite(coef)==0 || coef < cutoff) { // if coef is infinite, then:
					coef=0;
				}
			      	else { // if coef ends up larger than cutoff value, add the values to the current a[i][1]'s intensities.
					Grat3I[i] = Grat3I[i]  +   coef * exp(-M_PI * pow(((Grat3x[i]-dm * lambda * z12/period)/w2),2)) * cos(2 * M_PI * (dn/period) * (Grat3x[i]-dm * lambda * z12/period) * (1-z12/r2));
				    // Since a[i][1] etc. is actually the ix array, and arrays essentially get passed by reference, this is modifying the ix array.
				    continue;
				}
			}
		}
	}

	end = clock();
	diff =((double) (end - start));
	//printf("%f\n",diff); //printing the time spent in gp1.

}

void ( * gp2(double zloc, double el1x, double w1x, double r1x, double Grat3x[], double Grat3I[]))
{	
	time_t start, end; //starting a timer to get the time spent in function gp2
	start = clock();
	double diff=0;

	// get intensity profile after two grating
	double G2_x = sp.G2_x; //x position after second grating
	double r1y = r1x;
	double w1y = w1x;
	double el1y = el1x;
    
	int elecOrAtom = sp.elecOrAtom;
	double vel = sp.vel;
	int rows = sp.res;
	int xpnts = rows;
	double width = sp.height;
	double abszloc = zloc;
	int accountGrav = sp.accountGrav;
	double z12 = ((sp.G2_z) - (sp.G1_z));
	double z23 = (zloc -(sp.G2_z));
	double mytheta = sp.theta;
	double energy = sp.energy;
    	int rowsT =41;// rows of ReT and ImT array
    	double xstart = sp.xstart; // what is this? It's negative 200 microns.
    	double xend = sp.xend;
	double period = 0.0000001;// period of grating - 100 nanometers.
    	double wedgeangle = sp.wedgeangle;
    	int useimagecharge = sp.useimagecharge;	// whether or not to consider image charge effects. 0 for False.
	/* 
	 * A free parameter. Beta variable below depends on this. If beam is perpendicular to gratings, then tilt (and thus Beta) is 0.
	 * This is the twist about the x-axis.
	 */
    	double tilt =sp.tilt; 
    	double cutoff = sp.cutoff; 		// The point at which the intensity cuts off and is treated as 0.
	// G1 open fraction; how open the first grating is. With 0.4 open, a little over than half the muonium should pass through.
	double eta1 = sp.eta1; 
	// G2 open fraction; how open the second grating is.
    	double eta2 = sp.eta2;
    	double lambda = sqrt((1.5 * pow(10,-18))/(energy)); // wavelength we're working with of particles/waves
    	double res = sp.res; // This is the resolution we want this graph at.
    	double eta = width/period; // ratio of slit window 'height' to the period of the gratings
        double alpha = wedgeangle * M_PI/180;
    	double beta = tilt * M_PI; // 0 if beam is normal to gratings
    	double theta = M_PI * mytheta/180;
    	double d1=period; // period = period of gratings
    	double d2=period;
    	double z13 = z12  +  z23; // z distance between grating 1 and 3
    	double phi = 0;
    	double lim =5;
    	double _Complex coef;

    	/* THIS FUNCTION IS USING GSM MODEL FROM MCMORRAN, CRONIN 2008
    	 * IT MUST BE CHANGED IF YOU WANT TO USE OLDER MODEL FROM BREZGER 2003
	 */

    	double dn = 0;
    	double dm =0;
    	double m=0;
    	double n=0;
    	int a5 =0;
    	int b  =0;
    	int c5 =0;
    	int d5=0;
        double argument_d;
        double argument_f;
        double argument_p;
        double argument_v;
        double function_d;
        double function_v;

	/* Explanation of variables:
	 * w3 = GSM width of beam after the second grating.
	 * r3 = radius of GSM wavefront curvature after the second grating
	 * el3 = GSM beam coherence width after the second grating.
    	 */
    	double el3x = calculate_width(z13, r1x, el1x, w1x, el1x);	// z13 == G2z - G1z + zstart + 0 * zres; GSM coherence width in x-axis
    	double w3x = calculate_width(z13, r1x, el1x, w1x, w1x); 	// Beam width in x-axis
    	double v3x = v(z13,r1x,el1x,w1x); 				// Gaussian-Schell Model (GSM) radius of wavefront curvature in x-axis
    	double el3y = calculate_width(z13, r1y, el1y, w1y, el1y); 	// Coherence width in y-axis
    	double w3y = calculate_width(z13, r1y, el1y, w1y, w1y); 	// Beam width in y-axis
    	double v3y = v(z13,r1y,el1y,w1y); 				// radius of wavefront curvature in y-axis
    
    	int pos[41]={0}; // array of 41 elements

    	for (int i=0; i<rowsT; i++){ // rowsT = rows of ReT and ImT arrays = 41 for now
    		pos[i]=i-((rowsT-1)/2); // so this goes from -20 to 20
    	}
    
    	double ReT[41]={0}; 
    	// array of 41 0's for now.
    	// array of 41 0's for now.
    	int RealorIm = 1; 
    	ReTandImTgenerator(ReT,energy, elecOrAtom, RealorIm, vel, width, abszloc, accountGrav); 
    	double ImT[41]={0};
    	RealorIm = 2;
    	ReTandImTgenerator(ImT,energy, elecOrAtom, RealorIm, vel, width, abszloc, accountGrav); 
    
    	for (int i=0; i<rows; i++) { 
		// these are just the x-positions; a[i][1] will have the intensities of the beam at the xpositions.
        	Grat3x[i]= xstart + (i) * ((xend-xstart)/(xpnts-1));
    	}
    
    	double *phix;
    	double *phiI;
    	phix = (double*) calloc(rows, sizeof(double)); // it has same x-positions as Grat3x as shown below
    	phiI = (double*) calloc(rows, sizeof(double)); // it will affect Grat3I array
    
	for (int i=0; i<rows; i++) {
	phix[i]= xstart + (i) * ((xend-xstart)/(xpnts-1));
	}
    	
	

	for (int i=0; i<rows; i++) {
		for (int m1=-lim; m1<=lim; m1++) {
			for (int m2=-lim; m2<=lim; m2++) {
				for (int n1=-lim; n1<=lim; n1++) {
					for (int n2=-lim; n2<=lim; n2++) {
						dn =n1-n2;
						n = ((double)(n1 + n2))/2;
						dm = m1-m2;
						m = ((double)(m1 + m2))/2;
						a5 = (x2pnts(m1, (int  * )pos));
						b = (x2pnts(m2, (int  * )pos));
						c5 = (x2pnts(n1, (int  * )pos));
						d5 = (x2pnts(n2, (int  * )pos));

						//argument_d corresponds to the argument of equation 18b from McMorran & Cronin 2008. Note that y=0
                        			argument_d = -M_PI*(pow( phix[i]-lambda*z23*(n*cos(theta)/d2 + m*z13/(d1*z23) ),2 )/pow(w3x,2) + pow((n*sin(theta)*lambda)/(d2*w3y),2));
                        			//argument_f corresponds to the argument of equation 18c from McMorran & Cronin 2008. Note that y=0
                        			argument_f = -2*M_PI * phix[i] * ((dn*cos(theta)/d2)*(1-z23/v3x) + (dm/d1)*(1-z13/v3x));
                        			//argument_p corresponds to the argument of equation 18d from McMorran & Cronin 2008
                        			argument_p = (2*M_PI*lambda*z13*dm/d1)*(n*cos(theta)/d2 + m/d1)*(1-z13/v3x) +     (2*M_PI*lambda*z23*dn/d2)*( (m*cos(theta)/d1) * (1-z13/v3x) - (n*z23/d2)*(pow(cos(theta),2)/v3x)+pow(sin(theta),2)/v3y ) ;
                        			//argument_v corresponds to the argument of equation 18e from McMorran & Cronin 2008
                        			argument_v = -M_PI* pow(lambda*z23* (dn*cos(theta)/d2 + dm*z13/(d1*z23)),2)/pow(el3x,2) -M_PI*pow(dn*sin(theta)*lambda*z23/(d2*el3y),2);

						// 0 means ignore image charge effects, 1 means include image charge effects
						if (useimagecharge==0) {
						    coef = sinc(eta1 * M_PI * m1) +  0 * _Complex_I;
						    coef = coef * (sinc(eta1 * M_PI * m2)) +  0 * _Complex_I;
						}
						else { // assumes G1 is identical to G2
						    coef = ReT[a5]  +  ImT[a5] * _Complex_I; // 
						    coef = coef  *  ( (ReT[b]-ImT[b] * _Complex_I) );
						}
						
						coef = coef * (ReT[c5]  +  ImT[c5] * _Complex_I);
						coef = coef * (ReT[d5]  -  ImT[d5] * _Complex_I);
						
						if (((__real__ coef)>=cutoff) || ((__imag__ coef)>=cutoff)) {
						    
                        			function_d = exp(argument_d);                      
                        			function_v = exp(argument_v);
                        			phiI[i] = argument_f + argument_p;
                        
                        			Grat3I[i] = Grat3I[i]  +  ((__real__ coef) * cos(phiI[i]) - (__imag__ coef) * sin(phiI[i]))*function_d*function_v ;
						}
					}
				}
	    		}
		}
		
		
		
    	}

	free(phix);
	free(phiI);

	end = clock();
	diff =((double) (end - start));
	//printf("%f\n",diff); //printing the time spent in gp2.
}
