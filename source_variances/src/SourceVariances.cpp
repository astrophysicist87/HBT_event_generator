#include <iostream>
#include <fstream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <random>
#include <algorithm>

#include "SourceVariances.h"
#include "estimate_error.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

void SourceVariances::initialize_all(
	ParameterReader * paraRdr_in,
	const vector<EventRecord> & allEvents_in )
{
	// Load parameters
	paraRdr			= paraRdr_in;

	// Copy in records of all events
	allEvents		= allEvents_in;
	total_N_events	= allEvents.size();
	number_of_completed_events
					= 0;

	//Set header info
	// - particle information
	particle_mass 	= paraRdr->getVal("mass");
	// - some mode options
	bin_mode 		= paraRdr->getVal("bin_mode");
	q_mode	 		= paraRdr->getVal("q_mode");
	// - bin parameters
	bin_epsilon		= paraRdr->getVal("bin_epsilon");
	//Define various grid sizes
	// - pair momenta points at which to evaluate correlation function
	n_KT_pts 		= paraRdr->getVal("n_KT_pts");
	KT_min 			= paraRdr->getVal("KTmin");
	KT_max 			= paraRdr->getVal("KTmax");
	n_Kphi_pts 		= paraRdr->getVal("n_Kphi_pts");
	Kphi_min 		= -M_PI;
	Kphi_max 		= M_PI;
	n_KL_pts 		= paraRdr->getVal("n_KL_pts");
	KL_min 			= paraRdr->getVal("KLmin");
	KL_max 			= paraRdr->getVal("KLmax");
	// - relative momentum points at which to evaluate

	n_KT_bins 		= n_KT_pts - 1;
	n_Kphi_bins 	= n_Kphi_pts - 1;
	n_KL_bins 		= n_KL_pts - 1;

	KT_pts 			= vector<double> (n_KT_pts);
	Kphi_pts 		= vector<double> (n_Kphi_pts);
	KL_pts 			= vector<double> (n_KL_pts);

	linspace(KT_pts, KT_min, KT_max);
	linspace(Kphi_pts, Kphi_min, Kphi_max);
	linspace(KL_pts, KL_min, KL_max);

	px_bin_width 	= bin_epsilon;
	py_bin_width 	= bin_epsilon;
	pz_bin_width 	= bin_epsilon;

	// need to know these for binning particle pairs efficiently
	KT_bin_width 	= KT_pts[1]-KT_pts[0];
	Kphi_bin_width 	= Kphi_pts[1]-Kphi_pts[0];
	KL_bin_width 	= KL_pts[1]-KL_pts[0];

	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	S 				= vector<double> (K_space_size);
	x_S				= vector<double> (K_space_size);
	x2_S			= vector<double> (K_space_size);
	y_S				= vector<double> (K_space_size);
	y2_S			= vector<double> (K_space_size);
	z_S				= vector<double> (K_space_size);
	z2_S			= vector<double> (K_space_size);
	xo_S			= vector<double> (K_space_size);
	xo2_S			= vector<double> (K_space_size);
	xs_S			= vector<double> (K_space_size);
	xs2_S			= vector<double> (K_space_size);
	xl_S			= vector<double> (K_space_size);
	xl2_S			= vector<double> (K_space_size);
	t_S				= vector<double> (K_space_size);
	t2_S			= vector<double> (K_space_size);
	x_t_S			= vector<double> (K_space_size);
	y_t_S			= vector<double> (K_space_size);
	z_t_S			= vector<double> (K_space_size);
	x_y_S			= vector<double> (K_space_size);
	x_z_S			= vector<double> (K_space_size);
	y_z_S			= vector<double> (K_space_size);
	xo_t_S			= vector<double> (K_space_size);
	xs_t_S			= vector<double> (K_space_size);
	xl_t_S			= vector<double> (K_space_size);
	xo_xs_S			= vector<double> (K_space_size);
	xo_xl_S			= vector<double> (K_space_size);
	xs_xl_S			= vector<double> (K_space_size);


	// Initializations finished
	// Check number of events and proceed if non-zero
	if ( allEvents.size() == 0 )
		return;
	else
		out << "allEvents.size() = " << allEvents.size() << ": doing this file!" << endl;

	Compute_radii_qmode_1();

	return;
}

SourceVariances::~SourceVariances()
{
	//clear everything

	return;
}


void SourceVariances::Update_records( const vector<EventRecord> & allEvents_in )
{

	// Copy in new records of all events
	// (erases old event information)
	allEvents		= allEvents_in;
	total_N_events	+= allEvents.size();

	// Check number of events and proceed if non-zero
	if ( allEvents.size() == 0 )
		return;
	else
		cout << "allEvents.size() = " << allEvents.size() << ": doing this file!" << endl;

	// Compute 1D radii from source variances
	Compute_radii_qmode_1();

	return;
}


//End of file
