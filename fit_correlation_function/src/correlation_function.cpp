#include <iostream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "correlation_function.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

void Correlation_function::initialize_all(
	ParameterReader * paraRdr_in,
	string filepath_in )
{
	err << "Starting initializations:" << endl;
	// Load parameters
	paraRdr = paraRdr_in;

	//Set header info
	// - particle information
	particle_mass 	= paraRdr->getVal("mass");
	// - some parameters
	bin_mode 		= paraRdr->getVal("bin_mode");
	q_mode 			= paraRdr->getVal("q_mode");
	// - pair momenta points at which to interpolate HBT results
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
	//   correlation function
	n_qo_pts 		= paraRdr->getVal("n_qo_pts");
	n_qs_pts 		= paraRdr->getVal("n_qs_pts");
	n_ql_pts 		= paraRdr->getVal("n_ql_pts");
	// - step size in q directions
	delta_qo 		= paraRdr->getVal("delta_qo");
	delta_qs 		= paraRdr->getVal("delta_qs");
	delta_ql 		= paraRdr->getVal("delta_ql");

	err << " --> all parameters read in." << endl;

	// - minimum value in each q direction
	init_qo 		= -0.5*double(n_qo_pts-1)*delta_qo;
	init_qs 		= -0.5*double(n_qs_pts-1)*delta_qs;
	init_ql 		= -0.5*double(n_ql_pts-1)*delta_ql;

	n_qo_bins 		= n_qo_pts - 1;
	n_qs_bins 		= n_qs_pts - 1;
	n_ql_bins 		= n_ql_pts - 1;

	n_KT_bins 		= n_KT_pts - 1;
	n_Kphi_bins 	= n_Kphi_pts - 1;
	n_KL_bins 		= n_KL_pts - 1;

	KT_pts 			= vector<double> (n_KT_pts);
	Kphi_pts 		= vector<double> (n_Kphi_pts);
	KL_pts 			= vector<double> (n_KL_pts);

	qo_pts 			= vector<double> (n_qo_pts);
	qs_pts 			= vector<double> (n_qs_pts);
	ql_pts 			= vector<double> (n_ql_pts);

	linspace(KT_pts, KT_min, KT_max);
	linspace(Kphi_pts, Kphi_min, Kphi_max);
	linspace(KL_pts, KL_min, KL_max);

	linspace(qo_pts, init_qo, -init_qo);
	linspace(qs_pts, init_qs, -init_qs);
	linspace(ql_pts, init_ql, -init_ql);

	KT_bin_width 	= KT_pts[1]-KT_pts[0];
	Kphi_bin_width 	= Kphi_pts[1]-Kphi_pts[0];
	KL_bin_width 	= KL_pts[1]-KL_pts[0];

	// Vectors for HBT radii and intercept parameters (values and errors)
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;
	const int q_space_size = ( q_mode == 0 ) ?
								n_qo_bins*n_qs_bins*n_ql_bins
								: n_Q_bins;

	lambda_Correl 		= vector<double> (K_space_size);
	R2					= vector<double> (K_space_size);
	R2_out				= vector<double> (K_space_size);
	R2_side				= vector<double> (K_space_size);
	R2_long				= vector<double> (K_space_size);
	R2_outside			= vector<double> (K_space_size);
	R2_outlong			= vector<double> (K_space_size);
	R2_sidelong			= vector<double> (K_space_size);

	lambda_Correl_err	= vector<double> (K_space_size);
	R2_err				= vector<double> (K_space_size);
	R2_out_err			= vector<double> (K_space_size);
	R2_side_err			= vector<double> (K_space_size);
	R2_long_err			= vector<double> (K_space_size);
	R2_outside_err		= vector<double> (K_space_size);
	R2_outlong_err		= vector<double> (K_space_size);
	R2_sidelong_err		= vector<double> (K_space_size);

	// For the correlation function (and related error) itself
	correlation_function 		= vector<double> (K_space_size*q_space_size);
	correlation_function_error 	= vector<double> (K_space_size*q_space_size);

	err << " --> all vectors initialized." << endl;

	// Read in correlation function
	Load_correlation_function( filepath_in );

	// Read in correlation function
	Fit_correlation_function();

	return;
}

Correlation_function::~Correlation_function()
{
	//clear everything

	return;
}


void Correlation_function::Load_correlation_function( string filepath )
{
	err << "  --> Loading the correlation function from " << filepath << endl;
	// For timebeing, just skip header lines
	string line;
	ifstream infile( filepath.c_str() );
	/*while ( getline(infile, line) )
	{
		if ( line.front() == '#'
				and line.at(1) == '-' )	// make this the last comment line for now
		{
			break;
		}
	}*/
	// ALT: just read in fixed number of header lines and be done with it
	getline(infile, line);
	getline(infile, line);

	// Load correlation function itself
	int idx = 0;
	double dummy = 0.0;
	if ( q_mode == 0 )
	{
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		for (int iqo = 0; iqo < n_qo_bins; iqo++)
		for (int iqs = 0; iqs < n_qs_bins; iqs++)
		for (int iql = 0; iql < n_ql_bins; iql++)
		{
			infile >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
					>> dummy >> dummy >> dummy
					>> correlation_function[idx]
					>> correlation_function_error[idx];

			++idx;
		}
	}
	else if ( q_mode == 0 )
	{
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		for (int iQ = 0; iQ < n_Q_bins; iQ++)
		{
			infile >> dummy >> dummy >> dummy >> dummy
					>> dummy >> dummy >> dummy
					>> correlation_function[idx]
					>> correlation_function_error[idx];

			++idx;
		}
	}
	else
	{
		err << "fit_correlation_function(): q_mode = " << q_mode << " not supported!" << endl;
		exit(8);
	}


	err << "  --> Finished loading the correlation function from " << filepath << endl;

	return;
}



void Correlation_function::Output_HBTradii( string outHBT_filename )
{
	int prec = 4;
	int extrawidth = 12;

	FILE * pFile = fopen ( outHBT_filename.c_str(), "w" );

	// Print header information
	if ( q_mode == 0 )
	{
		fprintf ( pFile, "# K_T      K_phi      K_L      lambda      ");
		fprintf ( pFile, "R2o      R2s      R2l      R2os      R2ol      ");
		fprintf ( pFile, "R2sl      lambda(err)      R2o(err)      ");
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
							0.5*(KT_pts[iKT]+KT_pts[iKT+1]),
							0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]),
							0.5*(KL_pts[iKL]+KL_pts[iKL+1]) );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      %f      ",
							lambda_Correl[idx],
							R2_out[idx], R2_side[idx], R2_long[idx],
							R2_outside[idx], R2_outlong[idx], R2_sidelong[idx] );
				fprintf (  pFile,  "%f      %f      %f      %f      %f      %f      %f\n",
							lambda_Correl_err[idx],
							R2_out_err[idx], R2_side_err[idx], R2_long_err[idx],
							R2_outside_err[idx], R2_outlong_err[idx], R2_sidelong_err[idx] );


			++idx;
		}

	}
	else if ( q_mode == 1 )
	{
		fprintf ( pFile, "# K_T      K_phi      K_L      lambda      ");
		fprintf ( pFile, "R2      lambda(err)      R2o(err)\n" );

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
							lambda_Correl[idx], R2_out[idx],
							lambda_Correl_err[idx], R2_out_err[idx] );

			++idx;
		}

	}
	else
	{
		err << "fit_correlation_function(): q_mode = " << q_mode << " not supported!" << endl;
		exit(8);
	}

	return;
}


//End of file
