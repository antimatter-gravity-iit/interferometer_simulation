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

double ( * ReTandImTgenerator(double ReTorImTar[], double energy, int ReTorImT, double vel, double width, double abszloc, int accountGrav))
{
	double xstart = sp.xstart; // negative 200 microns
	double xend = sp.xend; //0.00020;      // positive 200 microns
  	double period = sp.grating_period; //0.0000001;// period of grating - 100 nanometers.
	double gravAccel = -9.8;    // acceleration due to gravity. 
  	double wedgeangle = sp.wedgeangle; //0; //Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
  	int useimagecharge = sp.useimagecharge; //0; // whether or not to consider image charge effects. 0 for False.
  	double tilt = sp.tilt;  //0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
  	double cutoff = sp.cutoff;// 0.000001; // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.
  	double eta1 = sp.eta1; //.4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it. Varname could be changed to better represent it.
  	double eta2 = sp.eta2; //.4; //G2 open fraction; how open the second grating is.
  	double thick = sp.thick; // 0.000000014; // 14 nanometers. Not (real) thickness of gratings, most likely. Gratings are actually 1 micrometer thick. This is used for the electron part of the code.
	double Gthick = sp.Gthick;// 1000; // thickness of gratings; 1 micrometer = 1000 nm, this is in nm on purpose (see function ReTgenerator) Varname could be better. Right now Gthick is used for the VdW effect for atoms.
	int rowsT = 41;// rows of ReT and ImT array
    	double resolution = sp.resolution;//1000; // This is the resolution we want this graph at. Varname could be better. MUST BE SAME AS IN MAIN!
    
	//values not included in simparam structure.  can be moved there but not entirely necessary 
	double chargeratio =0.0; //strength of image charge (units of e, electron charge); values of 0.03, 0.05, or more can be had //not used
    	float C3 = 2.0453e-2; // the VdW coefficient for hydrogen (assumed to be the same for muonium) //
    	double e_charge = 1.6021765e-19; // electric charge in Coulombs of electron (abs. value)
    	double Coulomb = 8.98755179e-9; // force; m^2/(Coulomb^-2)
    	double difPlancks = 6.58212e-13; // hbar in mev * s
    	double Plancks = 6.626068e-34; // Planck's constant
  
  	double eta = width/period; // ratio of 'height' of slit/windows in gratings to the period of the gratings
  	double nmvel = vel * 1e9;  // converting a m/s velocity to nm/s.
    	double alpha = wedgeangle * M_PI/180; // depends on wedgeangle above, which is a relatively free parameter. Appears to be bend of 'window' (slits in grating), if they bend forward or not.
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
        
        
        
      xmin= width * (1/resolution - cos(beta)/2); // minimum distance beam travels through slit; or maybe it's when the 1st order diffracted beams are going in diagonally, what is the min x?

        if (beta<=alpha) { // if the beam is very orthogonal to gratings (almost 90 degrees), or wedge angle is significant
            xmax=(width * cos(beta))/2-width/resolution;
        }
        else { // if beam is not very perpendicular to gratings, then it travels through the slit diagonally, covering more distance, more image charge interaction, etc.
            xmax= width  *  cos(beta)/2 - width/resolution  +  thick  *  (tan(alpha)-tan(beta));
        }
    }
    else { // if beta < 0; this time xmin changes, xmax is the same
        xmax = (width * cos(beta)/2)-width/resolution;

    	if (fabsl(beta)<=alpha) { // fabsl is for long doubles and returns a long double absolute value; once again, if the tilt isn't that bad, one bound (this time xmin) is just width * cos(beta)/2  +  width/res.
          xmin = -((width * cos(beta))/2) + width/resolution; 
        }

    	else { // if the beam is far from perpendicular to grating slits
            xmin = -((width * cos(beta))/2) + width/resolution - thick * (tan(alpha)-tan(beta));
        }
        
    }
    
    for(int n=-((rowsT-1)/2);n<=((rowsT-1)/2);n++) {
        for(ex=xmin; ex<xmax; ex +=width/resolution) { //copied from above; resolution = step resolution in x-axis. ex += height of window / steps = (40 nm / 1000)

          ////// THIS IS NOT IMPLEMENTED YET, ORIGINALLY IGOR PRO CODE THAT WAS IN HERE

          //// TODO: IMPLEMENT EFFECTS OF GRATING TILT ON INTERACTION
//            ph = 0
//            for(i = 0; i <=1; i  += 1)
//                beta *= -1
//                ex *= -1
//                ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                    /(2 * muonium_vel * Plancks) * thick/((width/2 + ex)
//                     * cos(alpha) * cos(beta))
//                if(sin(alpha + beta) != 0)
//                    ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                        /(2 * muonium_vel * Plancks * sin(alpha + beta))
//                         * log((width/2 + ex + thick * (tan(alpha) + tan(beta)))
//                        /(width/2 + ex))
//                else
//                    ph  += 2 * Pi * charge_ratio * e_charge^2 * Coulomb
//                        /(2 * muonium_vel * Plancks) * (width/2 + ex)
//                         * thick * cos(alpha)/cos(beta)
//                endif
//            endfor

          ////// END OF CODE THAT ISN'T IMPLEMENTED YET

          // ex is how far you are from the grating 'wall'
    
          // exnm is how far from the wall in nanometers.
          exnmleft = ex * 1.0e9;
          exnmright = (xmax - ex) * 1.0e9; 
          // fc is another electron thing; or the diffraction pattern? 2pi*n*x/period?
          fc = 2 * M_PI * n * ex/period; // so the first fc = 2  *  pi  *  -20  *  xmin / period, last fc = 2  *  pi  *  20  *  xmax / period. Looks like fc is a quantity proportional to distance from bottom/top of each 'window'/slit
            
          // abszloc is actually wrong; the code thinks the total z is 1m. It's actually closer to 2.8cm. Divide abszloc by 36.075 to get real value.
          double realzloc = abszloc / 36.075;

          // How much time has passed for a point in this system? x = vt, so t = x/v. x in this case is approximated by the z-location. We know v = vel.
          timeFreefall = (realzloc/ vel);

          // Both electrons and atoms will fall due to gravity. According to Dr. Daniel Kaplan's paper at arxiv.org/ftp/arxiv/papers/1308/1308.0878.pdf, the phase shift caused is 2 * pi * g * t^2 / d, where t is the time in free fall and d is the period of the gratings.
	if (accountGrav == 1)
		phGrav = (2 * M_PI * gravAccel * pow(timeFreefall, 2)) / period;  // phase shift due to gravity on particles

	else
		phGrav = 0;
         
          
          
	if (exnmleft == 0 | exnmright == 0)
		phM = 0;
 // phM is phase shift on Muonium/other neutral molecules due to Van der Waals effects through the gratings.
              	
	else
		phM = -C3 * Gthick / (difPlancks * nmvel * pow(exnmleft, 3)) -  -C3 * Gthick / (difPlancks * nmvel * pow(exnmright, 3));
                
	j=n + ((rowsT-1)/2); // j goes from 0 to 40 (right now rowsT = 41)
                
                
	if (ReTorImT == 1) // if it's the ReT array										
		ReTorImTar[j]  += cos(phM + fc + phGrav); // so fc and phM are both phase shifts; angles.

	else if (ReTorImT == 2) // if it's the ImT array
		ReTorImTar[j]  += sin(phM + fc + phGrav); // so fc and phM are both phase shifts; angles. 
  
          

            
        }
        
        
        
    }
    
    for (int i=0; i<rowsT; i++) { // must optimize all these for's that can probably be inside other for loops.
        
      ReTorImTar[i] = ReTorImTar[i]/resolution; // So is this some sort of normalization?
      
    }
    
    
    
}

