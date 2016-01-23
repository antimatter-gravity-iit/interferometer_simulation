// Gratings.c -> Holds the functions that determine the diffraction after each of the gratings in the system
// Mcomment - this file needs an overhaul.  It is basically incomprehensible, 
//            definitely immutable, and incredibly hard to debug.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "Gratings.h"
#include "BeamParams.h"
#include "Misc.h"
#include "PhaseShifts.h"
#include "complex.h"

double ( * gp0(double z,double r0,double el0, double w0, double Grat3x[], double Grat3I[],double energy, int rows, int col, int xpnts))
// get intensity profile 
{
    double xstart = -0.00020; // what is this? It's negative 200 microns.
    double xend = 0.00020;
    double pi = 3.14159265358979; // the constant irrational number pi.
    double w1;
    double jj;
    w1 = w(z,r0,el0,w0,energy); // width of beam incoming
    for(int i=0; i<rows; i++)
    {
        Grat3x[i]= xstart + (i) * ((xend-xstart)/(xpnts-1)); // current x-position at step i is put into a[i][0]
        jj = pow((Grat3x[i]/w1),2); //jj = (xpos/beamwidth)^2 
        Grat3I[i]=exp(-(pi * jj)); // a[i][1] is the intensity of the beam at the xposition at step i.
    }
}

double ( * gp1(double z12,double r1,double el1, double w1, double Grat3x[], double Grat3I[], double energy, int rows, int col, int elecOrAtom, double vel, int xpnts, double width, double abszloc, int accountGrav))
{
    // get intensity profile after one grating
    int rowsT =41;// rows of ReT and ImT array
    double xstart = -0.00020; // what is this? It's negative 200 microns.
    double xend = 0.00020;
    double pi = 3.14159265358979; // the constant irrational number pi.
    double period = 0.0000001;// period of grating - 100 nanometers.

    double wedgeangle = 0; // Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
    int useimagecharge = 0; // whether or not to consider image charge effects. 0 for False.
    double tilt =0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
    double cutoff = 0.000001; // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.
    double eta1 = .4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it.
    double eta2 = .4; //G2 open fraction; how open the second grating is.
    double coef; 
    double lim=5;
    double lambda = sqrt((1.5 * pow(10,-18))/(energy)); 
    // wavelength of what particles/waves we're working with; 
    double eta = width/period; // ratio of window 'height' to period of grating
    //double vel = pow(2 * energy * e_charge/e_mass,1/2); // electron velocity
    double alpha = wedgeangle * pi/180; // alpha and beta have been defined in almost every other function. Global variables? 
    double beta = tilt * pi; // defined in other functions too, same purpose.
    // Varnames could be better. w2, r2, el2? GSM width of beam, radius of GSM wavefront curvature, and GSM beam coherence width, respectively, after the first grating.
    double w2=w(z12,r1,el1,w1,energy); // width of beam between z1 and z2 (after grating 1)
    double r2 = v(z12,r1,el1,w1,energy); // radius of wavefront curvature between grating 1 and 2
    double el2 = el(z12, r1, el1, w1,energy); // beam coherence width
    int pos[41]={0};

    for (int i=0; i<rowsT; i++) // since rowsT is currently 41, pos[i] = -20 to 20.
    {
        pos[i]=i-((rowsT-1)/2);
    }

    double ReT[41]{0};
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
        for (int n=-lim; n<=lim; n++) 
        {
            for (int m=-lim; m<=lim; m++)
            {
                double dn =n-m; 
                // Varname could be better. n, m, dn, dm are all vague. I STILL don't know what they represent. Something to ask Dr. McMorran about?I know they're used to compute a part of the grating interaction, and I think they're distances from something, but not sure what. 
                double dm = (m + n)/2;

                if (useimagecharge==0) // if useimagecharge = 0, ignore image charge effects at G1. 
                {
                    coef = sinc(eta1 * pi * n)  *  (sinc(eta1 * pi * m) * pow((eta1), 2));
                    
                }

                else // if usechargeimage = 1, don't ignore image charge effects at G1.
                {
                    coef =  ReT[x2pnts(n, (int  * )pos)]  *  ReT[x2pnts(m,(int  * )pos)]  +  ImT[x2pnts(n,(int  * )pos)]  *  ImT[x2pnts(m,(int  * )pos)];
                }

                // lambda is the wavelength of our particles/waves
                coef = coef * exp(-pi * (dn * (lambda * z12/(pow(period * el2,2)))));
                // added isfinite macro in order to avoid inf values

                if (std::isfinite(coef)==0) // if coef is infinite, then:
                {
                    coef=0;
                }
              
                if (coef>=cutoff) // if coef ends up larger than cutoff value, add the values to the current a[i][1]'s intensities.
                {
                    
                    Grat3I[i] = Grat3I[i]  +   (coef * exp(-pi * pow(((Grat3x[i]-dm * lambda * z12/period)/w2),2) * cos(2 * pi * (dn/period) * (Grat3x[i]-dm * lambda * z12/period) * (1-z12/r2))));
                    // Since a[i][1] etc. is actually the ix array, and arrays essentially get passed by reference, this is modifying the ix array.
                    continue; // why have this here? Is there a point?
                }
            }
        }
    }
}

double ( * gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x, double Grat3x[], double Grat3I[], double energy, int rows, int col, int elecOrAtom, double vel, int xpnts, double width, double abszloc, int accountGrav))
{
    // get intensity profile after two grating
    int rowsT =41;// rows of ReT and ImT array
    double xstart = -0.00020; // what is this? It's negative 200 microns.
    double xend = 0.00020;
    double pi = 3.14159265358979; // the constant irrational number pi.
    double period = 0.0000001;// period of grating - 100 nanometers.
    double wedgeangle = 0; // Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
    int useimagecharge = 0; // whether or not to consider image charge effects. 0 for False.
    double tilt =0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.
    double cutoff = 0.000001; // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.
    double eta1 = .4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it.
    double eta2 = .4; //G2 open fraction; how open the second grating is.
    double lambda = sqrt((1.5 * pow(10,-18))/(energy)); // wavelength we're working with of particles/waves
    double res = 1000; // This is the resolution we want this graph at.
    double eta = width/period; // ratio of slit window 'height' to the period of the gratings
    // double vel = pow(2 * energy * e_charge/e_mass,1/2); electron velocity
    double alpha = wedgeangle * pi/180;
    double beta = tilt * pi; // 0 if beam is normal to gratings
    double theta = pi * mytheta/180;
    double d1=period; // period = period of gratings
    double d2=period;
    double z13 = z12  +  z23; // z distance between grating 1 and 3
    double phi = 0;
    double lim =5; // what is this?
    double _Complex coef;

    // THIS FUNCTION IS USING GSM MODEL FROM MCMORRAN, CRONIN 2008
    // IT MUST BE CHANGED IF YOU WANT TO USE OLDER MODEL FROM BREZGER 2003

    double dn = 0;
    double dm =0;
    double m=0;
    double n=0;
    int a5 =0;
    int b  =0;
    int c5 =0;
    int d5=0;
    
    double el3x = el(z13, r1x, el1x, w1x, energy);//G2z - G1z  +  zstart  +  0 * zres, r1, el1, w1; GSM coherence width in x-axis
    double w3x = w(z13,r1x,el1x,w1x, energy); // Beam width in x-axis
    double v3x = v(z13,r1x,el1x,w1x, energy); // Gaussian-Schell Model (GSM) radius of wavefront curvature in x-axis
    double el3y = el(z13,r1y,el1y, w1y, energy); // Coherence width in y-axis
    double w3y = w(z13,r1y,el1y,w1y, energy); // Beam width in y-axis
    double v3y = v(z13,r1y,el1y,w1y, energy); // radius of wavefront curvature in y-axis
    
    int pos[41]={0}; // array of 41 elements
    for (int i=0; i<rowsT; i++) // rowsT = rows of ReT and ImT arrays = 41 for now
    {
        pos[i]=i-((rowsT-1)/2); // so this goes from -20 to 20
    }
    
    double ReT[41]={0}; 
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
    phix = (double*) calloc(rows, sizeof(double)); // has same x-positions as Grat3x as shown below
    phiI = (double*) calloc(rows, sizeof(double)); // will affect Grat3I array;
    
    for (int i=0; i<rows; i++) {
      phix[i]= xstart + (i) * ((xend-xstart)/(xpnts-1)); // yep, just the xpositions so far.
    }
    
    for (int i=0; i<rows; i++)
    {
        for (int m1=-lim; m1<=lim; m1++) 
        {
            for (int m2=-lim; m2<=lim; m2++) 
            {
                for (int n1=-lim; n1<=lim; n1++) 
                {
                    for (int n2=-lim; n2<=lim; n2++) 
                    {
                        dn =n1-n2;
                        n = ((double)(n1 + n2))/2;
                        dm = m1-m2;
                        m = ((double)(m1 + m2))/2;
                        a5 = (x2pnts(m1, (int  * )pos));
                        b = (x2pnts(m2, (int  * )pos));
                        c5 = (x2pnts(n1, (int  * )pos));
                        d5 = (x2pnts(n2, (int  * )pos));
                        
                        if (useimagecharge==0) // 0 means ignore image charge effects, 1 means include image charge effects
                        {
                            coef = sinc(eta1 * pi * m1) +  0 * _Complex_I;
                            coef = coef * (sinc(eta1 * pi * m2)) +  0 * _Complex_I;
                        }
                        else // assumes G1 is identical to G2
                        {
                            coef = ReT[a5]  +  ImT[a5] * _Complex_I; // 
                            coef = coef  *  ( (ReT[b]-ImT[b] * _Complex_I) );
                        }
                        
                        coef = coef * (ReT[c5]  +  ImT[c5] * _Complex_I);
                        coef = coef * (ReT[d5]  +  ImT[d5] * _Complex_I);
                        // next factor responsible for the twist dependence of visibility
                        coef=coef * (exp(-pi * pow(((dn * sin(theta) * lambda * (z23))/(d2 * el3y)),2)));
                        coef=coef * (exp(-pi * pow((lambda * z23 * (dn * cos(theta) + dm * z13/z23)/(d1 * el3x)),2)));
                        // a[i][1] only has significant values if coef's real or imag parts are above the cutoff value. Otherwise a is returned as 0 intensity.
                        if (((__real__ coef)>=cutoff) || ((__imag__ coef)>=cutoff)) 
                        {
                            phi = dn * n * (1-z23/v3x) * pow((cos(theta)),2)  +  dn * n * (1-z23/v3y) * pow((sin(theta)),2)  +  dn * m * (1-z13/v3x) * cos(theta);
                            phi = phi  + (dm * n * (1-z13/v3x) * cos(theta)  +  dm * m * (z13/z23) * (1-z13/v3x));
                            phi = phi * (2 * pi * lambda * z23/(pow(d1,2)));
                            phi = phi - (2 * pi * dn * G2_x/d2);
                            phiI[i] = ((phi-(2 * pi * (phix[i])/d2) * (dn * cos(theta) * (1-z23/v3x)  +  dm * (1-z13/v3x))));
                            Grat3I[i] = Grat3I[i]  +  ((((__real__ coef) * cos(phiI[i]) - (__imag__ coef) * sin(phiI[i])) * exp(-pi * pow(((phix[i]-(lambda * z23/d1) * (n * cos(theta) + m * (z13/z23)))/w3x),2))));
                        }
                    }
                }
            }
        }
    }
    free(phix);
    free(phiI);
}
