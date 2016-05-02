//	SimplePlot.cc - simple plotting functions
//	if TEST is defined, also generates a simple main program to test it.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"

#include "SimplePlot.h"

TApplication *SimplePlot::app = 0;
int SimplePlot::nPlot = 0;

void SimplePlot::init()
{
	if(app != 0) return;

	int argc=1;
	char *argv[2]; argv[0] = (char *)"test"; argv[1] = 0;
	app = new TApplication("App",&argc,argv);
}

void SimplePlot::oneD(const char *title, const double y[], int nValues)
{
	if(nValues < 2) {
		fprintf(stderr,"ERROR SimplePlot::oneD: nValues < 2\n");
		return;
	}

	// create and fill x[] with indices into y[]
	double *x = new double[nValues];
	for(int i=0; i<nValues; ++i) x[i] = (double)i;

	graph(title,x,y,nValues);

	delete[] x;
}

void SimplePlot::oneD(const char *title, const float y[], int nValues)
{
	if(nValues < 2) {
		fprintf(stderr,"ERROR SimplePlot::oneD: nValues < 2\n");
		return;
	}

	// create and fill x[] with indices into y[]
	float *x = new float[nValues];
	for(int i=0; i<nValues; ++i) x[i] = (float)i;

	graph(title,x,y,nValues);

	delete[] x;
}

void SimplePlot::graph(const char *title, const double x[], const double y[],
								int nValues)
{
	init();

	// create a canvas with a unique name; it is drawn in its own window
	char name[64];
	sprintf(name,"Plot %d",++nPlot);
	TCanvas *c = new TCanvas(name,name, 400, 400);

	// create the graph (automatically goes into *c); this plots x vs y
	TGraph *g = new TGraph(nValues,x,y);
	g->SetTitle(title);
	g->GetXaxis()->SetTitle("Horizontal displacement x (meters)");
	g->GetYaxis()->SetTitle("Relative Intensity");
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->CenterTitle();
	g->Draw("APL"); // draw Axes, markers at Points, and Lines
	c->Update();

	// Run the Root event loop, until Quit ROOT is clicked
	fprintf(stderr,"Click on 'File/Quit ROOT' to proceed\n");
    app->Run(true);
	 //close the canvas and clean up
	c->Close();
	delete g;
	delete c;
}

void SimplePlot::graph(const char *title, const float x[], const float y[],
								int nValues)
{
	init();

	// create a canvas with a unique name; it is drawn in its own window
	char name[64];
	sprintf(name,"Plot %d",++nPlot);
	TCanvas *c = new TCanvas(name,name, 400, 400);

	// create the graph (automatically goes into *c); this plots x vs y
	TGraph *g = new TGraph(nValues,x,y);
	g->SetTitle(title);
	g->Draw("APL"); // draw Axes, markers at Points, and Lines
	c->Update();
	
	// Run the Root event loop, until Quit ROOT is clicked
	fprintf(stderr,"Click on 'File/Quit ROOT' to proceed\n");
	app->Run(true);
	// close the canvas and clean up
	c->Close();
	delete g;
	delete c;
}

void SimplePlot::twoD(const char *title, double value[],
		double Xmin, double Xmax, double Ymin, double Ymax,
		int NX, int NY, const char *options)
{
	init();

	// handle options
	if(options == 0) options = "colz,empty";
	std::string opt = std::string(",") + options + ",";
	int pos = opt.find(",empty,");
	if(pos < 0) pos = opt.find(",EMPTY,");
	if(pos >= 0) {
		opt.replace(pos,5,""); // remove ",empty"
		options = opt.c_str(); // (Root does not know it)
		for(int i=0; i<NX*NY; ++i)
            if(value[i] == 0.0) value[i] = 1.0E-300;
            //if(value[i] <= 0.05) value[i] = 0.0;//1.0E-300;
        
        gStyle->SetPalette(52,0);
        
	}
    
   // gStyle->SetPalette(52,0);


	// create a canvas with a unique name; it is drawn in its own window
	char name[64];
	sprintf(name,"Plot %d",++nPlot);
	TCanvas *c = new TCanvas(name,name, 400, 400);

	// create the plot (automatically goes into *c)
	TH2D *h = new TH2D(name, title, NX, Xmin, Xmax, NY, Ymin, Ymax);
	for(int ix=0; ix<NX; ++ix) { //  0   1    0     1   rows  rows
		for(int iy=0; iy<NY; ++iy) {
			h->SetBinContent(ix+1,iy+1,value[INDEX(ix,iy)]);
		}
	}
	double mean_x = h->GetMean(1);
	printf("Mean x = %.20f \n",mean_x);
	printf("Minimum measurement necessary to see gravitational effects: %.3fpm \n", fabs((mean_x)*1000000));       //convert um to nm

	h->Draw(options);
	h->GetXaxis()->SetTitle("Horizontal displacement x (micro meters)");
	h->GetYaxis()->SetTitle("Verticle displacement y (micro meters)");
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	c->Update();

	// Run the Root event loop, until Quit ROOT is clicked
	fprintf(stderr,"Click on 'File/Quit ROOT' to proceed\n");
	app->Run(true);

	// close the canvas and clean up
	c->Close();
	delete h;
	delete c;
}

void SimplePlot::setPalette(int ncolors, int colors[])
{
	gStyle->SetPalette(ncolors,colors);
}

#ifdef TEST


/**	main() is a simple test program for SimplePlot
 **/
int main(int argc, char *argv[])
{
	printf("oneD(double)\n");
	const int NVALUES=200;
	double values[NVALUES];
	for(int i=0; i<NVALUES; ++i)
		values[i] = sin(2.0*M_PI*i/50.0);
	SimplePlot::oneD("Test sin()",values,NVALUES);

	printf("oneD(float)\n");
	float floats[NVALUES];
	for(int i=0; i<NVALUES; ++i)
		floats[i] = cos(2.0*M_PI*i/50.0);
	SimplePlot::oneD("Test cos()",floats,NVALUES);

	printf("graph\n");
	float x[NVALUES], y[NVALUES];
	for(int i=0; i<NVALUES; ++i) {
		x[i] = cos(2.0*M_PI*i/31.0);
		y[i] = cos(2.0*M_PI*i/100.0);
	}
	SimplePlot::graph("Graph",x,y,NVALUES);

	printf("twoD\n");
	SimplePlot::selectGrayScale();
#define NX 100
#define NY 100
	double twod[NX*NY];
	for(int i=0; i<NX; ++i) {
		for(int j=0; j<NY; ++j) {
			twod[INDEX(i,j)] = 1E-20;
			if(j == 10) twod[INDEX(i,j)] = 0.1;
			if(j == 20) twod[INDEX(i,j)] = 0.2;
			if(j == 30) twod[INDEX(i,j)] = 0.3;
			if(j == 40) twod[INDEX(i,j)] = 0.4;
			if(j == 50) twod[INDEX(i,j)] = 0.5;
			if(j == 80) twod[INDEX(i,j)] = (double)i/100.0;
		}
	}
	SimplePlot::twoD("twoD",twod,0.0,1.0,0.0,1.0,NX,NY);

	printf("exiting\n");
	return 0;
}

#endif // TEST
