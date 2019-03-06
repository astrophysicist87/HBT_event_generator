#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "ParameterReader.h"
#include "Arsenal.h"
#include "EventRecord.h"
#include "ParticleRecord.h"

using namespace std;

const complex<double> i(0.0, 1.0);
const double hbarC = 0.19733;	//GeV*fm

class Correlation_function
{
	private:
		ParameterReader * paraRdr;

		//header info
		string particle_name;
		double particle_mass;

		int bin_mode, q_mode, fit_mode;

		int n_Q_pts;
		int n_KT_pts, n_Kphi_pts, n_KL_pts;
		int n_qo_pts, n_qs_pts, n_ql_pts;
		int n_Kx_pts, n_Ky_pts, n_Kz_pts;
		int n_qx_pts, n_qy_pts, n_qz_pts;

		int n_Q_bins;
		int n_KT_bins, n_Kphi_bins, n_KL_bins;
		int n_qo_bins, n_qs_bins, n_ql_bins;
		int n_Kx_bins, n_Ky_bins, n_Kz_bins;
		int n_qx_bins, n_qy_bins, n_qz_bins;

		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;
		double Kx_min, Kx_max, Ky_min, Ky_max, Kz_min, Kz_max;
		double qx_min, qx_max, qy_min, qy_max, qz_min, qz_max;

		double init_Q, delta_Q;
		double init_qo, init_qs, init_ql;
		double delta_qo, delta_qs, delta_ql;

		double KT_bin_width, Kphi_bin_width, KL_bin_width;
		double Kx_bin_width, Ky_bin_width, Kz_bin_width;
		double qx_bin_width, qy_bin_width, qz_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> Q_pts;
		vector<double> KT_pts, Kphi_pts, KL_pts;
		vector<double> qo_pts, qs_pts, ql_pts;
		vector<double> Kx_pts, Ky_pts, Kz_pts;
		vector<double> qx_pts, qy_pts, qz_pts;
		
		//miscellaneous
		string path;
		ostream & out;
		ostream & err;

		vector<double> lambda_Correl, R2, R2_out, R2_side, R2_long,
						R2_outside, R2_outlong, R2_sidelong;
		vector<double> lambda_Correl_err, R2_err, R2_out_err, R2_side_err, R2_long_err,
						R2_outside_err, R2_outlong_err, R2_sidelong_err;

		vector<double> denominator, correlation_function, correlation_function_error;
		vector<bool> denominator_cell_was_filled;
		vector<complex<double> > numerator;


	public:

		// Constructors, destructors, and initializers
		Correlation_function( ParameterReader * paraRdr_in,
								const string filename_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( paraRdr_in, filename_in ); };


		void initialize_all(ParameterReader * paraRdr_in,
								const string filename_in);

		~Correlation_function();

		////////////////////
		// Library functions
		inline int indexerK(int iKT, int iKphi, int iKL)
		{
			return ( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL );
		}
		////////////////////
		inline int indexer_Q_K(int iKT, int iKphi, int iKL, int iQ)
		{
			return (
					( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL )
							* n_Q_bins + iQ
					);
		}
		////////////////////
		inline int indexer(int iKT, int iKphi, int iKL, int iqo, int iqs, int iql)
		{
			return (
					( ( ( ( iKT * n_Kphi_bins + iKphi )
								* n_KL_bins + iKL )
								* n_qo_bins + iqo )
								* n_qs_bins + iqs )
								* n_ql_bins + iql
					);
		}
		////////////////////


		// Input
		void Load_correlation_function( string filepath );

		// Output fit results
		void Output_HBTradii( string filepath );

		// Fit routines
		void Fit_correlation_function();
		void Fit_correlation_function_Q();
		void find_minimum_chisq_correlationfunction_full( int iKT, int iKphi, int iKL );
		void find_minimum_chisq_correlationfunction_Q( int iKT, int iKphi, int iKL );

};

#endif
