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
#include "gauss_quadrature.h"

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

		double bin_epsilon;

		int bin_mode, q_mode, method_mode, BE_mode;
		int total_N_events, number_of_completed_events;

		int n_KT_pts, n_Kphi_pts, n_KL_pts;
		int n_qo_pts, n_qs_pts, n_ql_pts;
		int n_Kx_pts, n_Ky_pts, n_Kz_pts;
		int n_qx_pts, n_qy_pts, n_qz_pts;
		int n_Q_pts, n_qRP_pts, n_thq_pts;

		int n_KT_bins, n_Kphi_bins, n_KL_bins;
		int n_qo_bins, n_qs_bins, n_ql_bins;
		int n_Kx_bins, n_Ky_bins, n_Kz_bins;
		int n_qx_bins, n_qy_bins, n_qz_bins;
		int n_Q_bins;

		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;
		double Kx_min, Kx_max, Ky_min, Ky_max, Kz_min, Kz_max;
		double qx_min, qx_max, qy_min, qy_max, qz_min, qz_max;

		double qo_min, qs_min, ql_min;
		double qo_max, qs_max, ql_max;
		double delta_qo, delta_qs, delta_ql;
		double Q_min, Q_max, delta_Q;

		double KT_bin_width, Kphi_bin_width, KL_bin_width;
		double px_bin_width, py_bin_width, pz_bin_width;
		double Kx_bin_width, Ky_bin_width, Kz_bin_width;
		double qx_bin_width, qy_bin_width, qz_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> KT_pts, Kphi_pts, KL_pts;
		vector<double> qo_pts, qs_pts, ql_pts;
		vector<double> Kx_pts, Ky_pts, Kz_pts;
		vector<double> qx_pts, qy_pts, qz_pts;
		vector<double> Q_pts, x_pts, x_wts, ttheta_q_pts, ttheta_q_wts;

		// results
		double n_pair_numerator, n_pair_denominator;
		vector<double> correlation_function, correlation_function_error;
		static vector<double> numerator, denominator, numerator2, denominator2, numerator_denominator;
		static vector<double> numPair, numPair2, denPair, denPair2;
		static vector<double> numerator_numPair, denominator_denPair;
		static vector<bool> denominator_cell_was_filled;
		static vector<int> numerator_bin_count, denominator_bin_count;
		
		// miscellaneous
		string path;
		ostream & out;
		ostream & err;


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
		inline int indexer_qmode_1(int iKT, int iKphi, int iKL, int iQ)
		{
			return (
					( ( iKT * n_Kphi_bins + iKphi )
								* n_KL_bins + iKL )
								* n_Q_bins + iQ 
					);
		}
		////////////////////
		inline int indexer_qmode_1(int iKT, int iKphi, int iKL, int iQ, int iqRP, int ithq)
		{
			return (
					( ( ( ( iKT * n_Kphi_bins + iKphi )
								* n_KL_bins + iKL )
								* n_Q_bins + iQ )
								* n_qRP_pts + iqRP )
								* n_thq_pts + ithq
					);
		}
		////////////////////
		inline int indexer_qmode_1(int iKT, int iKphi, int iKL, int iQ, int iqRP, int ithq, int iqs_root)
		{
			return (
					( ( ( ( ( iKT * n_Kphi_bins + iKphi )
								* n_KL_bins + iKL )
								* n_Q_bins + iQ )
								* n_qRP_pts + iqRP )
								* n_thq_pts + ithq )
								* 2 + iqs_root
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

		void get_random_angles(int n_mixed_events, vector<double> & random_angles);

		// Separate parts of correlation function
		void Compute_numerator_and_denominator();

		void Compute_numerator_and_denominator_momentum_space_only();
		void Compute_numerator_and_denominator_momentum_space_only_q_mode_3D();
		void Compute_numerator_and_denominator_momentum_space_only_q_mode_1D();

		void Compute_numerator_and_denominator_methodMode0();
		void Compute_numerator_and_denominator_methodMode0_q_mode_3D();
		void Compute_numerator_and_denominator_methodMode0_q_mode_1D();

		void Compute_numerator_and_denominator_methodMode1();
		void Compute_numerator_and_denominator_methodMode1_q_mode_3D();
		void Compute_numerator_and_denominator_methodMode1_q_mode_1D();

		void Compute_numerator_and_denominator_methodMode2();
		void Compute_numerator_and_denominator_methodMode2_q_mode_3D();
		void Compute_numerator_and_denominator_methodMode2_q_mode_1D();

		// Correlation function itself
		void Compute_correlation_function();

		// Input/output
		void Output_correlation_function( string filename );
		void Output_correlation_function_q_mode_3D( string filename );
		void Output_correlation_function_q_mode_1D( string filename );

};

#endif
