/*
Current set of tasks.
- find new term for wave particle
- find unused variables
- flow chart of what everything does for drive
- create struct for waveparticle , charged particle
- insert main cmd line args for different predefined structs and user inputs.
- remove error testing and user input frames to replace with cmd line arguments.
*/

//***
//* Central Simulator
//* 
//*     Adam, put concise description here.
/* Icomment: This central program takes user input for:
			 electronic or atomic beam type,
			 inclusion or exclusion of gravitational acceleration,
			 resolution of output graph,
			 velocity of particles, 
			 pitch of gratings,
			 and intensity or final simulation paths,
			 and calculates simulation parameters in order to perform the simulation 
			 with the users desired arguments. 
*/
//*   
//* Code inspired by thesis by Dr. Benjamin McMorran
//* Put title of thesis and location here
//*
//* Collaborators:
//*     Arthur Romero, 16 July to 24 August 2015
//*         -- Original author
//*     Adam Denchfield, 20 October 2015 to present
//*         -- Additional parameters and functionality
//*     Melanie Cornelius (Dooley), 1 Nov 2015 to present
//*         -- Optimizations, standards, and readability edits
//*		Isaac Gewarges, January 2016 to present
//*			-- Optimizations, standards, and readability edits
//* 
//***

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cfloat>
#include <ctype.h>
#include <complex.h>
#include <string> 
#include <complex.h>
#include <limits.h> 		//fix for islimit error


#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"

#include "SimplePlot.h"          //Icomment: This program uses Root's TApplication.h and creates
						    	//  		a plot based on values passed into the function. 
#include "Misc.h"                //          contains functions of other miscellaneous tasks, 
							    //          such as taking care of boundary conditions, 
							    //          checking the max value in an array, etc.
#include "BeamParams.h"          //          This program contains the functions defining the GSM beam behavior
#include "Gratings.h"
#include "PhaseShifts.h"


/*
None of these variables are currently used.

//double mu_debroglie = 0.0000000005584;      // the result of the above equation  NOT USED IN MAIN PROGRAM
//double muonium_freq = 536890000000000000;   // muonium frequency
//double mu_lifetime = 0.0000022;             // half life average decay time of a muon, in seconds
//float C3 = 0.020453;                        // the VdW coefficient for hydrogen (assumed to be the same for muonium)


double Coulomb = 0.00000000898755179;       // force; m^2/(Coulomb^-2)
double pi = 3.14159265358979;               // the constant irrational number pi.
double const_e = 2.71828182845905;          // the irrational constant e.

double cutoff = 0.000001;                   // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.
*/

/* 
   from McMorran's thesis, a value epsilon = ratio of permittivity of grating material to permittivity of free space
   and the image charge is  + e  *  (epsilon - 1) / (epsilon  +  1)
   for an ideal conductor image charge =  + e (electric charge), the energy due to this strength of image charge 1 nm from surface is U = -0.75eV
*/
//int useimagecharge = 0;                     // whether or not to consider image charge effects. 0 for False. //not used in program.
/* not used
double eta1 = .4;                           //G1 open fraction; how open the first grating is. With .4 open, a little over than half the muonium should pass through it.
double eta2 = .4;                           //G2 open fraction; how open the second grating is.
double g_period = 0.0000001;                  // period of grating - 100 nanometers.
*/
double r0 = -4.04;                          //initial radius of wavefront curvature; comes from initial beam
double el0 = 0.000001;                      // initial coherence width; 50e-9 can also be used. Depends on initial beam
double w0 = 0.00003;                        // initial beam width -- this is probably assumed, the muonium beam width. Can also be: 2e-6, 1e-6, depends on initial beam

double G1_z = 0.000001;                     // It being 1 micron high is arbitrary, pretty sure. Also same as thickness of gratings. Top down view of gratings:

// Mcomment - what is this?  I'm very confused - this is messy and unclear, and the diagram is not described.

//                                                ---------------------   -> at G2_z = 1
//
//                                                ---------------------   -> at G1_z = 0.000001
//                                                  (source going up)

struct waveparticle {
double charge;
double charge_ratio;
float VdW_coef;
	
	
	
}waveparticle;

double G2_z = 1;                            // assumed to be 1 meter away on z-axis.
double G2_x = 0.00000005;                   //50 nm. Initial lateral offset of G2.

double theta = 0;                           // could be 0.05 or more. This is the twist between 1st and second gratings, in degrees. 2nd and 3rd grating are fixed to same rotational twist.

double thick = 0.000000014;                 // 14 nanometers. Not (real) thickness of gratings, most likely. Gratings are 1 micrometer thick.
double Gthick = 1000;                       // thickness of gratings; 1 micrometer = 1000 nm, this is in nm on purpose (see function ReTgenerator)

double wedgeangle = 0;                      // Grating wedge angle. The variable alpha below depends on this. This is a free parameter. Appears to be related to beam splitting.
// Going with the picture of the gratings above, this is if the horizontal grating has a rotation along the x-axis. i.e., twist around y-axis

double tilt =0;                             // A free parameter. Beta variable below depends on this. If beam is perp. to grating, then tilt (and thus Beta) are 0. This is the twist about the x-axis.

double res = 1000;                          // Resolution of the interaction in the gratings. Better varname?

double zstart = -0.1;                       // defining bounds of the grating structure? This definitely corresponds to a location, probably the bottom of the grating.
double zend = 2.1;                          // probably the end of the grating.
double xstart = -0.00020;                   // what is this? It's negative 200 microns.
double xend = 0.00020;
double ystart = -0.00011;                   // what is this? negative 110 microns?
double yend = 0.00011;


int col = 2;                                // colomns of ix array // remnant of the old code, even though it's not used (it's an argument for the gp0, gp1, gp2 functions, but doesn't do anything)
int rowsT =41;                              // rows of ReT and ImT arrays; used to calculate phase shift



int main(int argc, char *argv[]){ //added arguments
	// the main function
	/*arguments are currently in the order of:
argv[#]
	
	1. Account gravity? 1 = True, 0 = False.
	2. Electron beam or atom beam? Electron beam = 1, Atom beam = 2. //ELECTRON BEAM NOT MODELED!!!
	3. Resolution [300-400 recommended]. 
	4. Velocity of particles in m/s
	5. Pitch of gratings [in nm] 
	6. Total simulation [1]? or final interference pattern [2]?
	7. (if total simulation = 1), logscale[1] or normal scale[0]?
	
	*/ 
    
    // These are variables actually used.
    int elecOrAtom = 2;         // variable that takes user input and makes calculations on whether electron or atom beam is desired to simulate
    int resolution;             // determines how many x rows and y columns are displayed in the planes
    int xpnts;                  // x points of intensity 
    int ypnts;				    // y points of intensity
    int zpnts;					// z points of intensity
    int rows;					// max number of x points
    double vel; 				// velocity of particles
    double pitch;  				// distance from top of one grating to top of another
    double height;			    // max number of y points
    double max;                 // maximum intensity in each slice of the z-plane
    double *Grat3I;             // intensity array
    double *Grat3x;             // array of x position of intensity
    int simchoice;
    double *izx;
    int izxnumels;
    int izxsize;
    double zres;
    int zlocstart;
	int accountGrav;
    double energy = ((1.5 * pow(10,-18))/(pow(0.00000000001,2))) * (1); //energy of electron
	
	/*switch (atoi(argv[1])){
		case 0:
		accountGrav = atoi(argv[2]);
		elecOrAtom = atoi(argv[3]);
		xpnts= atoi(argv[4]);
		ypnts= xpnts;
		zpnts= xpnts;
		rows = xpnts;
		vel = atoi(argv[5]);
		pitch = atoi(argv[6]);
		height = (pitch / 2) / 1000000000;
		simchoice = atoi(argv[7]);

		//add statements for custom particle structure
		/*f(argv[2]!= NULL) accountGrav = atoi(argv[2]);
		else{
			printf("no gravity selection made. Cannot continue computation without this selection\n");
			exit(1);
		}
		if(argv[3]!= NULL) elecOrAtom = atoi(argv[3]);
		else{
			printf("no particle type selected. Cannot continue computation without this selection\n");
			exit(1);
		}
		if(argv[4]!= NULL){
			xpnts = atoi(argv[4]); //max of x vals to be calculated
			ypnts = atoi(argv[4]); //max of y vals to be calculated
			zpnts = atoi(argv[4]); //max of z vals to be calculated
			rows = atoi(argv[4]);  // rows of ix array
		}
		else{
			printf("no valid resolution entered, using default 300\n");
			xpnts = 300; //max of x vals to be calculated
			ypnts = 300; //max of y vals to be calculated
			zpnts = 300; //max of z vals to be calculated
			rows = 300;  // rows of ix array
		}
		if(argv[5]!= NULL){
			vel = atoi(argv[5]);
		}
		else{
			printf("no valid velocity entered, using default 6400 m/s\n");
			vel = 3600;
		}
		if(argv[6]!= NULL){
			pitch = atoi(argv[6]);
			height = (pitch / 2) / 1000000000;
		}
		else{
			printf("no valid pitch entered, using default 100nm\n");
			pitch = 0.0000001;
			height = (pitch / 2) / 1000000000;
		}
		if(argv[7]!= NULL){
			simchoice = atoi(argv[7]);
		}
		else{
			printf("no valid simulation choice entered, using default final interference pattern\n");
			simchoice = 2;
		}
		
		break;
		
		case 1:
		//charged photon particle to be simulated. 
		break;
		
		case 2:
		//x-ray photons to be simulated
		break;
		
		case 3:
		//hydrogen particle simulated
		break;
		
		case 4:
		//muonium particle simulated
		break;
		
//		default:
		// print incorrect argument passed
//		exit(1);
	}
	*/
	
    // Mcomment - we know it's main - you labeled it as such. What does it do?  What does it call?  What does it need to run?  What should the user change?  Etc.
//((1.5 * pow(10,-18))/(pow(lambda,2))) * (1); //15000// energy
	
	//Not used anymore due to cmd line arguments.
	accountGrav=atoi(argv[1]);
    // User input:
    /* Want gravity accounted for?
	printf("Do you want gravity accounted for? [1 for yes, 0 for no]: ");
    scanf("%d", &accountGrav);

	*/ // User input:
	// Beam type
	// 1 = electron beam, 2 = atom beam, default = 2
	/*
	//Not used anymore due to cmd line arguments.
	elecOrAtom = 2;
	printf("Is this an electron beam [1] or atom beam [2]? ");
	scanf("%d", &elecOrAtom);
	*/
	elecOrAtom =atoi(argv[2]);
	/*
    // User input:
	// resolution
    // Determines how many rows in the x and y planes.
    printf("How much resolution do you want in the simulation? [300-400 recommended]: ");
    scanf("%d", &resolution);
    xpnts = resolution; //max of x vals to be calculated
    ypnts = resolution; //max of y vals to be calculated
    zpnts = resolution; //max of z vals to be calculated
    rows = resolution;// rows of ix array
	*/
	resolution = atoi(argv[3]);
	xpnts = atoi(argv[3]);
	ypnts = atoi(argv[3]);
	zpnts = atoi(argv[3]);
    rows = resolution;
	
	/*
	// User input:
    // Velocity of particle or atom
	// Default velocity = 6400 m/s?
    printf("What is the velocity [m/s] of the atoms/electrons? (default 6400) : ");
    scanf("%lf", &vel);
	*/
	vel = atoi(argv[4]);
	/*
    //***
    //* User input:
    //* Pitch of grating (from top of one slit to top of next slit)
    //* Default is 100 nm 
    //***
    pitch = 0.0000001;  
    printf("Grating pitch (default is 100) (nm): ");
    scanf("%lf", &pitch);
    // Converting from nanometers to meters and from pitch to height
    */
	pitch = atoi(argv[5]);
	height = (pitch / 2) / 1000000000;
    
	/*
    //***
    //* User input:1
    //* Total intensity or just final intereference patern?
    //***
	printf("Do you want the total simulation [1] or the end interference pattern [2]? ");
    scanf("%d", &simchoice);
	*/
	simchoice = atoi(argv[6]);
	
    // Mcomment: you never ask the user for this.  Do you want to?
    
	
	int logchoice = 0; // if the user wants a total simulation, ask them if they want it to be scaled logarithmically so they can see where more of the particles go
	
    // initializing two arrays to contain the intensities and xpositions of each intensity -- Mcomment: I tried to fix this, but what is it actually saying?  It's the intensities of the x-positions and the intensities?  That doesn't make sense to me.
    Grat3I = (double*) calloc(resolution, sizeof(double));
    Grat3x = (double*) calloc(resolution, sizeof(double));
	
   
    // Mcomment: what is all of this?  Please describe.  Also, these do not have helpful
    // variable names - please fix.  It's very important for readability.  Just tell
    // me what you want them to be, and I'll make the change.
    izxnumels = xpnts * zpnts;
    izxsize = izxnumels * sizeof(double);
    izx = (double*) calloc(izxnumels, sizeof(double)); 
    zres = (zend-zstart)/zpnts; // step resolution used in computation

    // Following three are used to calculate GSM values at the first grating. (Gaussian-Schell Model)
    // Mcomment: GSM?  Did you define this elsewhere?  Avoid acronyms when and where possible.  More keystrokes are better.
    double w1=w(G1_z,r0,el0,w0,energy); 
    // okay, using a function w to initialize this variable w1... should change the names if possible. w1 = coherence width of beam. 
    double r1 = v(G1_z,r0,el0,w0,energy); 
    // same as above. r1 defined in the function.
    double el1=el(G1_z,r0,el0,w0,energy); 
    // same as above? All three different functions, but same parameters passed.
 
    // Mcomment - this is an example of embarassingly parallel computation.  Melanie will parallelize once cleaned up.
    // Mcomment - I cleaned up this indent structure.  Follow this indent structure in the future.
    //            The previous indent structure was outdated, unconventional, and unreadable.
    if (simchoice == 1) {
       //  printf("Do you want to have a logarithmic scale? [1 for yes, 2 for no] [Makes small intensities more visible]: ");
       // scanf("%d", &logchoice);
		logchoice = atoi(argv[7]);

        zlocstart = 0;
    }
    else if (simchoice == 2) {
        zlocstart = zpnts - 1;
    }

    for (int i=(zlocstart); i<zpnts; i++) {
        // i = 299 is just to get last row of the z if resolution is set to 300. 
        memset(Grat3x, 0, rows * sizeof(double)); 
        // so each time the loop repeats, you reset the array's positions and intensities to zero. 
        memset(Grat3I, 0, rows * sizeof(double));
             
        double zloc = zstart  +  i * zres; // where you are with respect to z. 
        
  
        // if, else if, etc. is to determine where you are on zloc, and depending on where you are, you interact with different gratings. 

        if (zloc > G2_z) { 
            // if the location is above G2_z [which is 1 currently]
        
            // Mcomment: This line is exactly and 100% why we MUST use better 
            // variable names.  This is an absolute mess.
            gp2(G2_z-G1_z, zloc-G2_z, theta, el1, w1, r1, el1, w1, r1, G2_x, Grat3x, Grat3I, energy, rows, col, elecOrAtom, vel, xpnts, height, zloc, accountGrav); 
            // gp2 returns two values, has to, since q is an array of 2 doubles
                                                                      
            max = maximumvalue(Grat3I, rows); 
            // the largest value of ix.

            ixgenerator(Grat3I, zloc, logchoice, rows); 
            // ixgenerator returns two values, has to, since q1 is an array of 2 doubles
 
         
        }
        else {
            // can probably be replaced, more efficiently, by elseifs. 
            if (zloc > G1_z) {
                gp1(zloc - G1_z, r1, el1, w1,Grat3x, Grat3I,energy,rows,col, elecOrAtom, vel, xpnts, height, zloc, accountGrav); 
                // if interacting with the first grating

                max = maximumvalue(Grat3I, rows); 
                // max value of ix here.
                ixgenerator(Grat3I, zloc, logchoice, rows); 
                // still ixgenerator is returning two values, since q1 is a pointer to an array of 2 doubles


            }
            else {
                // simple GSM propagation until it hits the first grating
                gp0(zloc, r0, el0, w0,Grat3x, Grat3I,energy,rows,col, xpnts); 
                // if at the origin?

                max = maximumvalue(Grat3I, rows); 
                // as above
                ixgenerator(Grat3I, zloc, logchoice, rows); 
                // as above

            }
        }   
        
        
   
        // REMEMBER: STILL INSIDE PREVIOUS FOR LOOP: for (i = 0; i < zpnts; i++ ) // Can optimize or make more efficient? Ranges or vectorize?
        // Mcomment - for now, if you have to mention this, it really indicates your code is unreadable.

        for (int j=0; j<rows; j++ ) {
        // perhaps this for loop could be avoided with clever programming? Its purpose appears to just translate over values.
            
            int f = rows * i + j; 
            //400 * i + j; still keeping track of location.
            
            izx[f] = Grat3I[j]; 
            // the f-th element of izx is set to be ix[j][1]; the intensity of the beam at the j-th point.
            

            // Grat3x[int(izxnumels - f)] = izxnumels - f; // xpos from 0 to 299.
            // Grat3I[int(izxnumels - f)] = izx[f]; // intensity for each expos 
            
printf("the value of izx is:%0.15f \t %d \t %f \n", izx[f],f, max); // debug print            
            // printing out izf, f, and the max value of ix.
        }
        // Mcomment - this isn't necessary to say.
        // There is some code that's approximately here in McMorran's thesis, but not here now. 
    }
    if (simchoice == 1) {
        free(Grat3x); 
        // free the memory used by this array, since the simulation is over, izx has all the data.
        free(Grat3I); 
        // same as above.

        SimplePlot::twoD("Intensity graph as particles diffract through gratings",izx,0.0,1.0,0.0,1.0,rows,rows); 
        // using ROOT to plot izx
    }
    else if(simchoice == 2) {
        SimplePlot::graph("Intensity along the last grating", Grat3x, Grat3I, rows);  
        // using ROOT to plot intensity vs. position at end of interferometer
        free(Grat3x); 
        // free the memory used by this array, since the simulation is over.
        free(Grat3I); 
        // same as above.
    }
    // free up the space used by the izx array.
    free(izx);
}
