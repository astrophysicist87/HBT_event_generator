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

		int bin_mode;

		int n_pT_pts, n_pphi_pts, n_pY_pts;
		int n_KT_pts, n_Kphi_pts, n_KL_pts;
		int n_qo_pts, n_qs_pts, n_ql_pts;
		int n_px_pts, n_py_pts, n_pz_pts;
		int n_Kx_pts, n_Ky_pts, n_Kz_pts;
		int n_qx_pts, n_qy_pts, n_qz_pts;

		int n_pT_bins, n_pphi_bins, n_pY_bins;
		int n_KT_bins, n_Kphi_bins, n_KL_bins;
		int n_qo_bins, n_qs_bins, n_ql_bins;
		int n_px_bins, n_py_bins, n_pz_bins;
		int n_Kx_bins, n_Ky_bins, n_Kz_bins;
		int n_qx_bins, n_qy_bins, n_qz_bins;

		double pT_min, pT_max, pphi_min, pphi_max, pY_min, pY_max;
		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;
		double px_min, px_max, py_min, py_max, pz_min, pz_max;
		double Kx_min, Kx_max, Ky_min, Ky_max, Kz_min, Kz_max;
		double qx_min, qx_max, qy_min, qy_max, qz_min, qz_max;

		double init_qo, init_qs, init_ql;
		double delta_qo, delta_qs, delta_ql;

		double pT_bin_width, pphi_bin_width, pY_bin_width;
		double KT_bin_width, Kphi_bin_width, KL_bin_width;
		double px_bin_width, py_bin_width, pz_bin_width;
		double Kx_bin_width, Ky_bin_width, Kz_bin_width;
		double qx_bin_width, qy_bin_width, qz_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> pT_pts, pphi_pts, pY_pts;
		vector<double> KT_pts, Kphi_pts, KL_pts;
		vector<double> qo_pts, qs_pts, ql_pts;
		vector<double> px_pts, py_pts, pz_pts;
		vector<double> Kx_pts, Ky_pts, Kz_pts;
		vector<double> qx_pts, qy_pts, qz_pts;

		vector<double> dN_pTdpTdpphidpY;
		
		//miscellaneous
		string path;
		ostream & out;
		ostream & err;

		vector<double> lambda_Correl, R2_out, R2_side, R2_long,
						R2_outside, R2_outlong, R2_sidelong;
		vector<double> lambda_Correl_err, R2_out_err, R2_side_err, R2_long_err,
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
		inline int indexer(int ipT, int ipphi, int ipY)
		{
			return ( ( ipT * n_pphi_bins + ipphi )
							* n_pY_bins + ipY );
		}
		////////////////////
		inline int indexerK(int iKT, int iKphi, int iKL)
		{
			return ( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL );
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
		inline double get_q0(double m, double qo, double qs, double ql, double KT, double KL)
		{
			double xi2 = m*m + KT*KT + KL*KL + 0.25*(qo*qo + qs*qs + ql*ql);

			return ( sqrt(xi2 + qo*KT + ql*KL) - sqrt(xi2 - qo*KT - ql*KL) );
		}
		////////////////////


		// Input
		void Load_correlation_function( string filepath );

		// Output fit results
		void Output_HBTradii( string filepath );

		// Fit routines
		void Fit_correlation_function();
		void find_minimum_chisq_correlationfunction_full( int iKT, int iKphi, int iKL );
};

#endif
