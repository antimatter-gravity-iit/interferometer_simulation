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
	double gravAccel = -9.8;    // acceleration due to gravity. 
  	double tilt = sp.tilt;  //0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
  	double eta1 = sp.eta1; //.4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it. Varname could be changed to better represent it.
  	double eta2 = sp.eta2; //.4; //G2 open fraction; how open the second grating is.
    
	//values not included in simparam structure.  can be moved there but not entirely necessary 
	double chargeratio =0.0; //strength of image charge (units of e, electron charge); values of 0.03, 0.05, or more can be had //not used
    	float C3 = 2.0453e-2; // the VdW coefficient for hydrogen (assumed to be the same for muonium) //
    	double e_charge = 1.6021765e-19; // electric charge in Coulombs of electron (abs. value)
    	double Coulomb = 8.98755179e-9; // force; m^2/(Coulomb^-2)
    	double difPlancks = 6.58212e-13; // hbar in mev * s
    	double Plancks = 6.626068e-34; // Planck's constant
  
  	double eta = sp.slit_height/sp.grating_period; // ratio of 'height' of slit/windows in gratings to the period of the gratings
  	double nmvel = sp.particle_velocity * 1e9;  // converting a m/s velocity to nm/s.
    	double alpha = sp.wedgeangle * M_PI/180; // depends on wedgeangle above, which is a relatively free parameter. Appears to be bend of 'window' (slits in grating), if they bend forward or not.
    	double beta = tilt * M_PI; // depends on tilt angle, = 0 if beam is normal to gratings
    	double exnmleft; // how many nm from the left side of each slit are we?
    	double exnmright; // how many nm from the right side of each slit are we?
    	long double xmin; // beginning of path of wave through the slit
    	long double xmax; // end of path of wave/beam through the slit
    	float fc;
    	float phE; // phase shift if dealing with electrons
    	long double phM; // phase shift if dealing with neutral atoms/molecules
    	double phGrav; // phase shift due to gravity.
    	float ex;
    	double timeFreefall;
    	int j;

    if (beta>=0){ // if the beam is not normal/perpendicular to the gratings it encounters
        
        
        
      xmin= sp.slit_height * (1/sp.resolution - cos(beta)/2); // minimum distance beam travels through slit; or maybe it's when the 1st order diffracted beams are going in diagonally, what is the min x?

        if (beta<=alpha) { // if the beam is very orthogonal to gratings (almost 90 degrees), or wedge angle is significant
            xmax=(sp.slit_height * cos(beta))/2-sp.slit_height/sp.resolution;
        }
        else { // if beam is not very perpendicular to gratings, then it travels through the slit diagonally, covering more distance, more image charge interaction, etc.
            xmax= sp.slit_height  *  cos(beta)/2 - sp.slit_height/sp.resolution  +  sp.thick  *  (tan(alpha)-tan(beta));
        }
    }
    else { // if beta < 0; this time xmin changes, xmax is the same
        xmax = (sp.slit_height * cos(beta)/2)-sp.slit_height/sp.resolution;

    	if (fabsl(beta)<=alpha) { // fabsl is for long doubles and returns a long double absolute value; once again, if the tilt isn't that bad, one bound (this time xmin) is just sp.slit_height * cos(beta)/2  +  sp.slit_height/res.
          xmin = -((sp.slit_height * cos(beta))/2) + sp.slit_height/sp.resolution; 
        }

    	else { // if the beam is far from perpendicular to grating slits
            xmin = -((sp.slit_height * cos(beta))/2) + sp.slit_height/sp.resolution - sp.thick * (tan(alpha)-tan(beta));
        }
        
    }
    
    for(int n=-((sp.rowsT-1)/2);n<=((sp.rowsT-1)/2);n++) {
        for(ex=xmin; ex<xmax; ex +=sp.slit_height/sp.resolution) { //copied from above; resolution = step resolution in x-axis. ex += height of window / steps = (40 nm / 1000)

          ////// THIS IS NOT IMPLEMENTED YET, ORIGINALLY IGOR PRO CODE THAT WAS IN HERE

          //// TODO: IMPLEMENT EFFECTS OF GRATING TILT ON INTERACTION
//            ph = 0
//            for(i = 0; i <=1; i  += 1)
//                beta *= -1
//                ex *= -1
//                ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                    /(2 * muonium_vel * Plancks) * thick/((sp.slit_height/2 + ex)
//                     * cos(alpha) * cos(beta))
//                if(sin(alpha + beta) != 0)
//                    ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                        /(2 * muonium_vel * Plancks * sin(alpha + beta))
//                         * log((sp.slit_height/2 + ex + thick * (tan(alpha) + tan(beta)))
//                        /(sp.slit_height/2 + ex))
//                else
//                    ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                        /(2 * muonium_vel * Plancks) * (sp.slit_height/2 + ex)
//                         * thick * cos(alpha)/cos(beta)
//                endif
//            endfor

          ////// END OF CODE THAT ISN'T IMPLEMENTED YET

          // ex is how far you are from the grating 'wall'
    
          // exnm is how far from the wall in nanometers.
          exnmleft = ex * 1.0e9;
          exnmright = (xmax - ex) * 1.0e9; 
          // fc is another electron thing; or the diffraction pattern? 2pi*n*x/period?
          fc = 2 * M_PI * n * ex/sp.grating_period; // so the first fc = 2  *  pi  *  -20  *  xmin / period, last fc = 2  *  pi  *  20  *  xmax / period. Looks like fc is a quantity proportional to distance from bottom/top of each 'window'/slit

          // How much time has passed for a point in this system? x = vt, so t = x/v. x in this case is approximated by the z-location. We know v = vel.
          timeFreefall = current_z_position/ sp.particle_velocity;

          // Both electrons and atoms will fall due to gravity. According to Dr. Daniel Kaplan's paper at arxiv.org/ftp/arxiv/papers/1308/1308.0878.pdf, the phase shift caused is 2 * pi * g * t^2 / d, where t is the time in free fall and d is the period of the gratings.
	if (sp.account_gravity == 1)
		phGrav = (2 * M_PI * gravAccel * pow(timeFreefall, 2)) / sp.grating_period;  // phase shift due to gravity on particles

	else
		phGrav = 0;
         
          
          
	if (exnmleft == 0 | exnmright == 0)
		phM = 0;
 // phM is phase shift on Muonium/other neutral molecules due to Van der Waals effects through the gratings.
              	
	else
		phM = -C3 * sp.grating_thickness / (difPlancks * nmvel * pow(exnmleft, 3)) -  -C3 * sp.grating_thickness / (difPlancks * nmvel * pow(exnmright, 3));
                
	j=n + ((sp.rowsT-1)/2); // j goes from 0 to 40 (right now sp.rowsT = 41)
                
                
	if (ReTorImT == 1) // if it's the ReT array										
		ReTorImTar[j]  += cos(phM + fc + phGrav); // so fc and phM are both phase shifts; angles.

	else if (ReTorImT == 2) // if it's the ImT array
		ReTorImTar[j]  += sin(phM + fc + phGrav); // so fc and phM are both phase shifts; angles. 
  
          

            
        }
        
        
        
    }
    
    for (int i=0; i<sp.rowsT; i++) { // must optimize all these for's that can probably be inside other for loops.
        
      ReTorImTar[i] = ReTorImTar[i]/sp.resolution; // So is this some sort of normalization?
      
    }
    
    
    
}

