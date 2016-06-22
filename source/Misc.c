// Misc.c - contains functions of other miscellaneous tasks, such as taking care of boundary conditions, checking the max value in an array, etc.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "Misc.h"
//#include "PhaseShifts.h"
//#include "Gratings.h"
//#include "BeamParams.h"

double maximumvalue(double arr[], int rows){    
    // finds the max value of the array
    double m =0;
    for (int i = 0; i < rows; i++ ){
        if (m<arr[i]){
            m = arr[i];// finds the maximum value of the array
        }
    }
    return m;
}


double sinc(double x){
    // this function avoids a division by x=0, sinc is defined to be sin(x)/x which will definied as 1 in such case
    double sinc_value=1;
    if (x!=0){
        sinc_value = (sin(x))/x;
    }
    return(sinc_value);
}


/* 
 * The function 'ixgenerator' is the x direction intensity calculator. It normalizes and compares intensities
 * to cutoff value, then determines which value to input to the array of x intensitites. Its arguments are:
 * 	intensity at the time of calculation (an initially empty array);
 * 	zlocation;
 * 	choice of scale for the plot; and
 * 	number of rows.
 * The function only modifies the array fed to it; it doesn't return any value. 
 */

double *ixgenerator(double a[], double zloc, int logchoice, int rows) 
{
    double cutoff = 0.000001; 
    // at what point does the intensity cut off and be treated as 0. Can also be 5e-5 like in McMorran thesis. Or 0.001.cutoff = 0.000001;
    double const_e = 2.71828182845905; 
    // the irrational constant e.
    double mu_lifetime = 0.0000022; 
    // half life average decay time of a muon, in seconds
    int rowsT =41;
    // rows of ReT and ImT arrays
    double max = maximumvalue(a, rows); 
    // maximum intensity value for a slice in z
    double min = cutoff;
    // the minimum value the intensity can be at before being set to 0; cutoff is a global variable
    // right now gratings are treated as being 0.5m apart. Divide zloc by 36.075 to get the proper distances, 1.4cm, between gratings.
    double realzloc = zloc / 36.075;

    for (int i=0; i<rows; i++) {
        /*
	 * HALF-LIFE DECAY - means normalizing the intensity at the point i to what the intensity should be after a certain time t.
	 * This will be approximated by substituting distance, zloc, for t; x = vt, v = 6300, so t = x/6300.
	 * In this case it's z instead of x. Since the electrons should become more and more as the antimuons decay, we actually want
	 * the electron presence to start out small, then go up.
	 */

        // Modeling electrons:
        // a[i][1] = a[i][1] * (1 - pow(const_e, (-1 * realzloc / 6300)/mu_lifetime)); // why is mu_lifetime included for an electron modeling array?

        // Modeling muonium: 
        a[i] = a[i] * (pow(const_e, (-1 * realzloc / 6300)/mu_lifetime));
        
        if (a[i]<min){ // could also have min = 0
            a[i]=0;
        }

        // a[i][1] = a[i][1]/max;//divides each element of the array by max, normalize it
        if (logchoice == 1 && a[i] > 0)
        {
            a[i] = log(a[i]/min); // it tries to make sure all elements are above 1, then takes the logarithm to scale the stuff better
        }
    }
}

int x2pnts(int value, int  * arr)
{
    int rowsT =41;
    // rows of ReT and ImT arrays
    for (int i = 0; i < rowsT; i++)
    {
        if (value ==  * (arr + i))
        {
            return(i);            
        }
    }
    printf("Error! n or m value does not match with any value of x2pnt array\n");
}
