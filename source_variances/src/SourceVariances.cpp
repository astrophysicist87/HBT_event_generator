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
#include <string>

#include "SourceVariances.h"
#include "estimate_error.h"
#include "Arsenal.h"
#include "Stopwatch.h"
#include "matrix.h"

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

	KT_bins 		= vector<double> (n_KT_bins);
	Kphi_bins 		= vector<double> (n_Kphi_bins);
	KL_bins 		= vector<double> (n_KL_bins);

	for (int iKT = 0; iKT < n_KT_bins; ++iKT)
		KT_bins[iKT] = 0.5*( KT_pts[iKT] + KT_pts[iKT+1] );
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
		Kphi_bins[iKphi] = 0.5*( Kphi_pts[iKphi] + Kphi_pts[iKphi+1] );
	for (int iKL = 0; iKL < n_KL_bins; ++iKL)
		KL_bins[iKL] = 0.5*( KL_pts[iKL] + KL_pts[iKL+1] );

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

	detHBT2D		= vector<double> (K_space_size);
	detHBT3D		= vector<double> (K_space_size);
	lambda_Correl	= vector<double> (K_space_size);
	R2				= vector<double> (K_space_size);
	R2o				= vector<double> (K_space_size);
	R2s				= vector<double> (K_space_size);
	R2l				= vector<double> (K_space_size);
	R2os			= vector<double> (K_space_size);
	R2ol			= vector<double> (K_space_size);
	R2sl			= vector<double> (K_space_size);

	// Initializations finished
	// Check number of events and proceed if non-zero
	if ( allEvents.size() == 0 )
		return;
	else
		out << "allEvents.size() = " << allEvents.size() << ": doing this file!" << endl;

	Compute_radii_qmode_3D();

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
	Compute_radii_qmode_3D();

	return;
}

void SourceVariances::Output_HBTradii( string outHBT_filename )
{
	int prec = 4;
	int extrawidth = 12;

	FILE * pFile = fopen ( outHBT_filename.c_str(), "w" );

	// Print header information
	//if ( q_mode == 0 )
	//{
		fprintf ( pFile, "# K_T      K_phi      K_L      ");
		fprintf ( pFile, "R2o      R2s      R2l      R2os      R2ol      ");
		fprintf ( pFile, "R2sl      det{2D}(R2ij)      det{3D}(R2ij)      R2o(err)      ");
		fprintf ( pFile, "R2s(err)      R2l(err)      R2os(err)      ");
		fprintf ( pFile, "R2ol(err)      R2sl(err)\n" );

		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							KT_bins[iKT], Kphi_bins[iKphi], KL_bins[iKL] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      %f      %f      ",
							R2o[idx], R2s[idx], R2l[idx],
							R2os[idx], R2ol[idx], R2sl[idx],
							detHBT2D[idx], detHBT3D[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f\n",
							0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );


			++idx;
		}

	//}
	/*else if ( q_mode == 1 )
	{
		fprintf ( pFile, "# K_T      K_phi      K_L      lambda      ");
		fprintf ( pFile, "R2      lambda(err)      R2(err)\n" );

		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							0.5*(KT_pts[iKT]+KT_pts[iKT+1]),
							0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]),
							0.5*(KL_pts[iKL]+KL_pts[iKL+1]) );
				fprintf (  pFile,  "%f      %f      %f      %f\n",
							lambda_Correl[idx], R2[idx],
							lambda_Correl_err[idx], R2_err[idx] );

			++idx;
		}

	}
	else
	{
		err << "source_variances(): q_mode = " << q_mode << " not supported!" << endl;
		exit(8);
	}*/

	return;
}


void SourceVariances::Output_source_moments( string outSM_filename, string coords )
{
	int prec = 4;
	int extrawidth = 12;

	FILE * pFile = fopen ( outSM_filename.c_str(), "w" );

	// Print header information
	if ( coords == "XYZ" )
	{

		fprintf ( pFile, "# K_T      K_phi      K_L      x      y      z      t      ");
		fprintf ( pFile, "x2      y2      z2      t2      ");
		fprintf ( pFile, "xy      xz      xt      yz      yt      zt      \n");
		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							KT_bins[iKT], Kphi_bins[iKphi], KL_bins[iKL] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							x_S[idx], y_S[idx], z_S[idx], t_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							x2_S[idx], y2_S[idx], z2_S[idx], t2_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      \n",
							x_y_S[idx], x_z_S[idx], x_t_S[idx],
							y_z_S[idx], y_t_S[idx], z_t_S[idx] );


			++idx;
		}

	}
	else if ( coords == "OSL" )
	{

		fprintf ( pFile, "# K_T      K_phi      K_L      xo      xs      xl      t      ");
		fprintf ( pFile, "xo2      xs2      xl2      t2      ");
		fprintf ( pFile, "xoxs      xoxl      xot      xsxl      xst      xlt      \n");
		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							KT_bins[iKT], Kphi_bins[iKphi], KL_bins[iKL] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							xo_S[idx], xs_S[idx], xl_S[idx], t_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							xo2_S[idx], xs2_S[idx], xl2_S[idx], t2_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      \n",
							xo_xs_S[idx], xo_xl_S[idx], xo_t_S[idx],
							xs_xl_S[idx], xs_t_S[idx], xl_t_S[idx] );


			++idx;
		}

	}
	else
	{
		err << "source_variances(): coords = " << coords << " not supported!" << endl;
		exit(8);
	}

	return;
}


void SourceVariances::Output_source_variances( string outSV_filename, string coords )
{
	int prec = 4;
	int extrawidth = 12;

	FILE * pFile = fopen ( outSV_filename.c_str(), "w" );

	// Print header information
	if ( coords == "XYZ" )
	{

		fprintf ( pFile, "# K_T      K_phi      K_L      ");
		fprintf ( pFile, "<x2>      <y2>      <z2>      <t2>      ");
		fprintf ( pFile, "<xy>      <xz>      <xt>      <yz>      <yt>      <zt>      \n");
		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							KT_bins[iKT], Kphi_bins[iKphi], KL_bins[iKL] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							x2_S[idx]-x_S[idx]*x_S[idx], y2_S[idx]-y_S[idx]*y_S[idx],
							z2_S[idx]-z_S[idx]*z_S[idx], t2_S[idx]-t_S[idx]*t_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      \n",
							x_y_S[idx]-x_S[idx]*y_S[idx], x_z_S[idx]-x_S[idx]*z_S[idx],
							x_t_S[idx]-x_S[idx]*t_S[idx], y_z_S[idx]-y_S[idx]*z_S[idx],
							y_t_S[idx]-y_S[idx]*t_S[idx], z_t_S[idx]-z_S[idx]*t_S[idx] );


			++idx;
		}

	}
	else if ( coords == "OSL" )
	{

		fprintf ( pFile, "# K_T      K_phi      K_L      ");
		fprintf ( pFile, "<xo2>      <xs2>      <xl2>      <t2>      ");
		fprintf ( pFile, "<xoxs>      <xoxl>      <xot>      <xsxl>      <xst>      <xlt>      \n");
		fprintf ( pFile, "#----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------" );
		fprintf ( pFile, "----------------------------------------\n" );

		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
				fprintf (  pFile,  "%f      %f      %f      ",
							KT_bins[iKT], Kphi_bins[iKphi], KL_bins[iKL] );
				fprintf (  pFile,  "%f      %f      %f      %f      ",
							xo2_S[idx]-xo_S[idx]*xo_S[idx], xs2_S[idx]-xs_S[idx]*xs_S[idx],
							xl2_S[idx]-xl_S[idx]*xl_S[idx], t2_S[idx]-t_S[idx]*t_S[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      \n",
							xo_xs_S[idx]-xo_S[idx]*xs_S[idx], xo_xl_S[idx]-xo_S[idx]*xl_S[idx],
							xo_t_S[idx]-xo_S[idx]*t_S[idx], xs_xl_S[idx]-xs_S[idx]*xl_S[idx],
							xs_t_S[idx]-xs_S[idx]*t_S[idx], xl_t_S[idx]-xl_S[idx]*t_S[idx] );


			++idx;
		}

	}
	else
	{
		err << "source_variances(): coords = " << coords << " not supported!" << endl;
		exit(8);
	}

	return;
}


//End of file
