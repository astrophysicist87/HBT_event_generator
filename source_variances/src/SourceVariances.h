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

constexpr complex<double> i(0.0, 1.0);
constexpr double hbarC = 0.19733;	//GeV*fm

class SourceVariances
{
	private:
		ParameterReader * paraRdr;

		//header info
		string particle_name;
		double particle_mass;

		double bin_epsilon;

		int bin_mode, q_mode, method_mode, BE_mode;
		int total_N_events, number_of_completed_events;
		
		int n_mix_minimum;

		int n_KT_pts, n_Kphi_pts, n_KL_pts;
		int n_Kx_pts, n_Ky_pts, n_Kz_pts;

		int n_KT_bins, n_Kphi_bins, n_KL_bins;
		int n_Kx_bins, n_Ky_bins, n_Kz_bins;

		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;
		double Kx_min, Kx_max, Ky_min, Ky_max, Kz_min, Kz_max;

		double KT_bin_width, Kphi_bin_width, KL_bin_width;
		double px_bin_width, py_bin_width, pz_bin_width;
		double Kx_bin_width, Ky_bin_width, Kz_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> KT_pts, Kphi_pts, KL_pts;
		vector<double> Kx_pts, Ky_pts, Kz_pts;

		vector<double> S, x_S, x2_S, y_S, y2_S, z_S, z2_S, xo_S, xo2_S, xs_S, xs2_S, xl_S, xl2_S, t_S, t2_S;
		vector<double> x_t_S, y_t_S, z_t_S, x_y_S, x_z_S, y_z_S, xo_t_S, xs_t_S, xl_t_S, xo_xs_S, xo_xl_S, xs_xl_S;

		
		// miscellaneous
		string path;
		ostream & out;
		ostream & err;


	public:

		// Constructors, destructors, and initializers
		SourceVariances( ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( paraRdr_in, allEvents_in ); };


		void initialize_all(ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in);

		~SourceVariances();

		////////////////////
		// Library functions
		inline int indexerK(int iKT, int iKphi, int iKL)
		{
			return ( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL );
		}
		////////////////////

		// For subsequent chunks of events
		void Update_records( const vector<EventRecord> & allEvents_in );

		void Compute_radii_qmode_3D();
		void Average_source_moments();
		void Set_radii();

};

#endif
