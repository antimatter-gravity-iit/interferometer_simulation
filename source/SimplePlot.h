/* 
 * SimplePlot.h
 * Copyright (C) 2016 Antimatter Gravity Interferometer Group, Illinois Institute of Technology (IIT). 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * 				*
 *
 * Code inspired by thesis by Dr. Benjamin McMorran
 * Electron Diffraction and Interferometry Using Nanostructures
 * http://gradworks.umi.com/33/52/3352633.html
 *
 * For a list of collaborators see README.md and CREDITS.
 *				
 *				*
 * 
 * DESCRIPTION:
 * Simple plotting functions using ROOT libraries.
 *
 */

#ifndef SIMPLEPLOT_H
#define SIMPLEPLOT_H

#include "TApplication.h"

/**	class SimplePlot contains functions to make Root plots.
 *
 *	Each function creates a Root plot in a new window.
 *	Each window blocks until "File/Quit ROOT" is clicked.
 *	The plots scale automatically.
 **/
class SimplePlot {
	static TApplication *app;
	static int nPlot;
	static void init();
public:
	/// oneD() plots a 1-d array.
	/// The plot is displayed in a new window.
	/// These functions block until "File/Quit ROOT" is clicked.
	/// The x axis is simply the integers [0,nValues).
	static void oneD(const char *title, const double y[], int nValues);
	static void oneD(const char *title, const float y[], int nValues);

	/// graph() plots a graph of (x[i],y{i]), with lines between adjacent
	/// points. x[] and y[] must be the same type (float or double), and
	/// have the same number of elements.
	/// The plot is displayed in a new window.
	/// These functions block until "File/Quit ROOT" is clicked.
	static void graph(const char *title, const double x[], const double y[],
								int nValues);
	static void graph(const char *title, const float x[], const float y[],
								int nValues);

	/// twoD() plots a 2-D array as a 2-d histogram. But to permit
	/// the caller to specify the size of the array, the actual array
	/// is 1-D, indexed by the INDEX(IX,IY) macro:
	///	#include "SimplePlot.h"
	///	#define NX 20
	///	#define NY 50
	///	double value[NX*NY];
	///	value[INDEX(ix,iy)] = ...
	///	NOTE: 0<=ix<NX, 0<=iy<Ny
	/// The plot is displayed in a new window.
	/// This function blocks until "File/Quit ROOT" is clicked.
	///
	/// options is a string containing a comma-separated set from the
	/// options given for THistPainter (space or no separataor not allowed).
	/// The most common options are:
	///	col	draw colored box for each bin, using the current palette
	///	colz	as col, displaying the palette
	///	surf2	draw surface plot using colors from the current palette
	///	cont4	draw contour plot using colors from the current palette
	///	scat	draw a scatter plot
	/// NEW OPTION:
	///	empty	draw 0.0 (empty) bins as 1.0E-300 (black in gray-scale)
	///		NOTE: This changes the 0.0 entries in the value array.
	/// if options== (the default)0, "colz,empty" is used.
	static void twoD(const char *title, double value[],
		double Xmin, double Xmax, double Ymin, double Ymax,
		int nx, int ny, const char *options=0);
#define INDEX(IX,IY) ((IX) + (NX)*(IY))

	/// Selects the palette -- see Root documentation for
	/// TColor::SetPalette()
	static void setPalette(int ncolors, int colors[]);

	/// select a gray-scale palette
	static void selectGrayScale() { setPalette(52,0); } // see Root docs
};

#endif
