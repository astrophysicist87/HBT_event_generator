#ifndef HBTEG_H
#define HBTEG_H

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

class HBT_event_generator
{
	private:
		ParameterReader * paraRdr;

		//header info
		string particle_name;
		double particle_mass;

		int bin_mode, total_N_events, number_of_completed_events;

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

		vector<double> denominator, correlation_function, correlation_function_error;
		vector<bool> denominator_cell_was_filled;
		vector<double> numerator, numerator2, denominator2, numerator_denominator;


	public:

		// Constructors, destructors, and initializers
		HBT_event_generator( ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( paraRdr_in, allEvents_in ); };


		void initialize_all(ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in);

		~HBT_event_generator();

		////////////////////
		// Library functions
		inline int bin_function( double datapoint, const vector<double> & points )
		{
			int result = (int)( ( datapoint - points[0] )
							* double( points.size()-1 )
							/ ( points[points.size()-1] - points[0] ) );

			// Assume uniform bin-widths for now
			return ( result );
		}
		////////////////////
		inline double bin_function( double p, double q, double bw, int mode )
		{
			return (
					/*double(mode==0) **/ double(abs(p-q) <= 0.5*bw) / bw
					// mode==1 too slow for now
					//+ double(mode==1) * exp( -(p-q)*(p-q) / (bw*bw) ) / ( sqrt(M_PI)*bw )
					);
		}
		////////////////////
		inline int indexerp(int ipT, int ipphi, int ipY)
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

		// For subsequent chunks of events
		void Update_records( const vector<EventRecord> & allEvents_in );


		// Functions to compute single-particle spectra
		//void Compute_spectra();
		// Sub-methods
		//void Compute_dN_pTdpTdpphidpY();
		//
		//void Compute_dN_pTdpTdpphi();
		//void Compute_dN_2pipTdpTdpY();
		//void Compute_dN_dpphidpY();
		//
		//void Compute_dN_2pipTdpT();
		//void Compute_dN_dpphi();
		//void Compute_dN_2pidpY();
		// total multiplicity
		//void Compute_N();


		// Separate parts of correlation function
		//void Compute_numerator(vector<complex<double> > & in_numerator);
		//void Compute_denominator(vector<double> & in_denominator);
		void Compute_numerator_and_denominator_with_errors(
							vector<double> & in_numerator,
							vector<double> & in_numerator2,
							vector<double> & in_denominator,
							vector<double> & in_denominator2,
							vector<double> & in_numerator_denominator
							);


		// Correlation function itself
		void Compute_correlation_function();


		// Input/output
		void Output_correlation_function();

};

#endif
