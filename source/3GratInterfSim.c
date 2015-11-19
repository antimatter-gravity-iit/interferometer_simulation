//
//  2Darray_test2.c
//  
//
//  Created by Arthur Romero on 7/16/15.
//  Modified by Adam Denchfield afterwards. Oct. 20th, 2015.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string>

#include "complex.h"

// Mcomment: are these custom-made header files?  They have very undescriptive names.
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"

#include "SimplePlot.hh"
#include "Misc.hh"
#include "BeamParams.hh"
#include "Gratings.hh"
#include "PhaseShifts.hh"


// Mcomment: I'll also show you how to do section headers.

//global variables:
// Mcomment - there are a whooooooole lot of these.  How many are actually necessary?  I understand if many are - just checking.

// Mcomment: should you move to more convenient units?  This is *very* error-prone with so many zeros.
double e_charge = 0.00000000000000000016021765; // electric charge in Coulombs of electron (abs. value)
double e_mass = 0.00000000000000000000000000000091093819; // mass of an electron
// Mcomment: what is difPlancks? I see the explanation to the right, but I'm still unclear about the name.
double difPlancks = 0.000000000000658212; // hbar in mev * s; used for VdW calculations
double Plancks = 0.0000000000000000000000000000000006626068; // Planck's constant
double mu_mass = 0.000000000000000000000000000188353; // kg, 1.8835E-28 kg

// double mu_debroglie_wavelength = Plancks / (mu_mass * muonium_vel) // De Broglie Wavelength of the muonium atom
double mu_debroglie = 0.0000000005584; // the result of the above equation
double muonium_freq = 536890000000000000; // muonium frequency
double mu_lifetime = 0.0000022; // half life average decay time of a muon, in seconds
float C3 = 0.020453; // the VdW coefficient for hydrogen (assumed to be the same for muonium)

double Coulomb = 0.00000000898755179; // force; m^2/(Coulomb^-2)
double pi = 3.14159265358979; // the constant irrational number pi.
double const_e = 2.71828182845905; // the irrational constant e.

double cutoff = 0.000001; // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.

// Mcomment - I'm going to show you how to do block comments eventually.  There are many ways, and believe me, they make life better.
// Mcomment - also, we will use the standard terminal width for text.  More things I'll fix but show you. =)
// all the following helps compute the interference pattern throughout a two-grating interferometer
// izx is a slice of this interference in the x-z plane


double chargeratio =0.0; //strength of image charge (units of e, electron charge); values of 0.03, 0.05, or more can be had
// from McMorran's thesis, a value epsilon = ratio of permittivity of grating material to permittivity of free space
// and the image charge is  + e  *  (epsilon - 1) / (epsilon  +  1)
// for an ideal conductor image charge =  + e (electric charge), the energy due to this strength of image charge 1 nm from surface is U = -0.75eV

int useimagecharge = 0; // whether or not to consider image charge effects. 0 for False.

double eta1 = .4; //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it.
double eta2 = .4; //G2 open fraction; how open the second grating is.
double period = 0.0000001;// period of grating - 100 nanometers.

double r0 = -4.04;//initial radius of wavefront curvature; comes from initial beam
double el0 = 0.000001;// initial coherence width; 50e-9 can also be used. Depends on initial beam
double w0 = 0.00003;// initial beam width -- this is probably assumed, the muonium beam width. Can also be: 2e-6, 1e-6, depends on initial beam

double G1_z = 0.000001; // It being 1 micron high is arbitrary, pretty sure. Also same as thickness of gratings. Top down view of gratings:

// Mcomment - what is this?  I'm very confused - this is messy and unclear, and the diagram is not described.

//                                                ---------------------   -> at G2_z = 1
//
//                                                ---------------------   -> at G1_z = 0.000001
//                                                  (source going up)

double G2_z = 1; // assumed to be 1 meter away on z-axis.
double G2_x = 0.00000005; //50 nm. Initial lateral offset of G2.

double theta = 0; // could be 0.05 or more. This is the twist between gratings, in degrees

double thick = 0.000000014; // 14 nanometers. Not (real) thickness of gratings, most likely. Gratings are 1 micrometer thick.
double Gthick = 1000; // thickness of gratings; 1 micrometer = 1000 nm, this is in nm on purpose (see function ReTgenerator)

double wedgeangle = 0; // Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
// Going with the picture of the gratings above, this is if the horizontal grating has a rotation along the x-axis. i.e., twist around y-axis

double tilt =0; // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.

double res = 1000; // Resolution of the interaction in the gratings. Better varname?

double zstart = -0.1; // defining bounds of the grating structure? This definitely corresponds to a location, probably the bottom of the grating.
double zend = 2.1; // probably the end of the grating.
double xstart = -0.00020; // what is this? It's negative 200 microns.
double xend = 0.00020;
double ystart = -0.00011; // what is this? negative 110 microns?
double yend = 0.00011;

int elecOrAtom = 2; // whether or not we are dealing with an electron or a neutral atom/molecule, this is just the default, to go with atom (2)

int col = 2;// colomns of ix array // remnant of the old code, even though it's not used (it's an argument for the gp0, gp1, gp2 functions, but doesn't do anything)
int rowsT =41;// rows of ReT and ImT arrays; used to calculate phase shift



int main()
// the main function
// Mcomment - we know it's main - you labeled it as such. What does it do?  What does it call?  What does it need to run?  What should the user change?  Etc.
{
  double width = 0.000000050; // 50 nanometers is the default. It looks like 'height' of each slit in the gratings. 

    double energy = ((1.5 * pow(10,-18))/(pow(0.00000000001,2))) * (1);//((1.5 * pow(10,-18))/(pow(lambda,2))) * (1); //15000// energy
    

    double max; // what is this? The maximum intensity in each slice of the z-plane.
    // izx is a slice of the interference in the x-z plane

    int accountGrav; // ask the user if they want gravity accounted for.
    printf("Do you want gravity accounted for? [1 for yes]: ");
    scanf("%d", &accountGrav);

    // ask the user what kind of resolution they want for the simulation. This determines how many rows in the x and y planes.

    printf("Is this an electron beam [1] or atom beam [2]? ");
    scanf("%d", &elecOrAtom);

    int resolution;
    printf("How much resolution do you want in the simulation? [300-400 recommended]: ");
    scanf("%d", &resolution);

    double vel; // velocity of atom in question
    printf("What is the velocity [m/s] of the atoms/electrons? ");
    scanf("%lf", &vel);

    printf("What is the height of the slits in the gratings (nm)? ");
    scanf("%lf", &width);
    width = width / 1000000000;// turning nm into m.

    // initializes the two arrays to contain the intensities and xpositions of each intensity.
   
    double *Grat3I;
    double *Grat3x;
    Grat3I = (double*) calloc(resolution, sizeof(double));
    Grat3x = (double*) calloc(resolution, sizeof(double));

    
    int xpnts = resolution; // how many x, y, z points to consider. Resolution of computation in x, y, z directions.
    int ypnts = resolution;
    int zpnts = resolution;
    int rows = resolution;// rows of ix array

// ask the user if they want the total intensity simulation or just the final interference pattern
    int simchoice;
    printf("Do you want the total simulation [1] or the end interference pattern [2]? ");
    scanf("%d", &simchoice);

    int logchoice = 0; // if the user wants a total simulation, ask them if they want it to be scaled logarithmically so they can see where more of the particles go


   
    int izxnumels = xpnts * zpnts;
    int izxsize = izxnumels * sizeof(double);
    double *izx;
    izx = (double*) calloc(izxnumels, sizeof(double)); 
    
    double zres = (zend-zstart)/zpnts; // step resolution used in computation

    // Following three are used to calculate GSM values at the first grating. 
    // Mcomment: GSM?  Did you define this elsewhere?  Avoid acronyms when and where possible.  More keystrokes are better.
    double w1=w(G1_z,r0,el0,w0,energy); // okay, using a function w to initialize this variable w1... should change the names if possible. w1 = coherence width of beam. 
    double r1 = v(G1_z,r0,el0,w0,energy); // same as above. r1 defined in the function.
    double el1=el(G1_z,r0,el0,w0,energy); // same as above? All three different functions, but same parameters passed.
 
    //MAIN LOOP -- computes intensity profile at each step of optical axis
    // Mcomment - excellent comment here.
    // Mcomment - this is an example of embarassingly parallel computation.  Melanie will parallelize once cleaned up.
    // Mcomment - I'm not a fan of this indenting structure.  Really hard to follow.  Leave first { in-line?  Indent all comments to same level?   
    int zlocstart;
    if (simchoice == 1)
      {
        printf("Do you want to have a logarithmic scale? [1 for yes] [Makes small intensities more visible]: ");
        scanf("%d", &logchoice);

        zlocstart = 0;
      }
    else if (simchoice == 2)
      {
        zlocstart = zpnts - 1;
      }

    for (int i=(zlocstart); i<zpnts; i++) // i = 299 is just to get last row of the z
    {
      memset(Grat3x, 0, rows * sizeof(double)); // so each time the loop repeats, you reset the array's positions and intensities to zero. 
      memset(Grat3I, 0, rows * sizeof(double));
             
        double zloc = zstart  +  i * zres; // where you are with respect to z. 
        
  
        // if, else if, etc. is to determine where you are on zloc, and depending on where you are, you interact with different gratings. 

        if (zloc > G2_z) // if the location is above G2_z [which is 1 currently]
        {
            gp2(G2_z-G1_z, zloc-G2_z, theta, el1, w1, r1, el1, w1, r1, G2_x, Grat3x, Grat3I, energy, rows, col, elecOrAtom, vel, xpnts, width, zloc, accountGrav); // gp2 returns two values, has to, since q is an array of 2 doubles
                                                                      
            max = maximumvalue(Grat3I, rows); // the largest value of ix.

            ixgenerator(Grat3I, zloc, logchoice, rows); // ixgenerator returns two values, has to, since q1 is an array of 2 doubles
 
         
        }
        else                      // can probably be replaced, more efficiently, by elseifs. 
        {
            if (zloc > G1_z)
            {
                gp1(zloc - G1_z, r1, el1, w1,Grat3x, Grat3I,energy,rows,col, elecOrAtom, vel, xpnts, width, zloc, accountGrav); // if interacting with the first grating

                max = maximumvalue(Grat3I, rows); // max value of ix here.
                ixgenerator(Grat3I, zloc, logchoice, rows); // still ixgenerator is returning two values, since q1 is a pointer to an array of 2 doubles


            }
            else // simple GSM propagation until it hits the first grating
            {
                gp0(zloc, r0, el0, w0,Grat3x, Grat3I,energy,rows,col, xpnts); // if at the origin?

                max = maximumvalue(Grat3I, rows); // as above
                ixgenerator(Grat3I, zloc, logchoice, rows); // as above

            }
        }
        
        
   
        // REMEMBER: STILL INSIDE PREVIOUS FOR LOOP: for (i = 0; i < zpnts; i++ ) // Can optimize or make more efficient? Ranges or vectorize?
        // Mcomment - I'll look at efficiency soon - for now, if you have to mention this, it really indicates your code is unreadable.

        for (int j=0; j<rows; j++ ) // perhaps this for loop could be avoided with clever programming? Its purpose appears to just translate over values.
        {
            
            int f = rows * i + j; //400 * i + j; still keeping track of location.
            
            izx[f] = Grat3I[j]; // the f-th element of izx is set to be ix[j][1]; the intensity of the beam at the j-th point.
            

           // Grat3x[int(izxnumels - f)] = izxnumels - f; // xpos from 0 to 299.
           // Grat3I[int(izxnumels - f)] = izx[f]; // intensity for each expos 
            
            printf("the value of izx is:%0.15f \t %d \t %f \n", izx[f],f, max); // printing out izf, f, and the max value of ix.
            
        }
 
        // Mcomment - this isn't necessary to say.
        // There is some code that's approximately here in McMorran's thesis, but not here now. 

    
    }
 

    if (simchoice == 1)
      {
        free(Grat3x); // free the memory used by this array, since the simulation is over, izx has all the data.
        free(Grat3I); // same as above.

        SimplePlot::twoD("Intensity graph as particles diffract through gratings",izx,0.0,1.0,0.0,1.0,rows,rows); // using ROOT to plot izx
      }

    else if(simchoice == 2)
      {
            SimplePlot::graph("Intensity along the grating", Grat3x, Grat3I, rows);  // using ROOT to plot intensity vs. position at end of interferometer
            free(Grat3x); // free the memory used by this array, since the simulation is over.
            free(Grat3I); // same as above.
      }

    // free up the space used by the izx array.
    free(izx);

}

