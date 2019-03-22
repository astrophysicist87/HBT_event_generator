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

#include "HBT_event_generator.h"
#include "estimate_error.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

// Declare vectors statically to allow OpenMP parallelization
// (possible without this, but more cumbersome)
vector<double> HBT_event_generator::numerator,
			   HBT_event_generator::denominator;
vector<double> HBT_event_generator::numerator2,
			   HBT_event_generator::denominator2,
			   HBT_event_generator::numerator_denominator;
vector<double> HBT_event_generator::numPair,
			   HBT_event_generator::numPair2,
			   HBT_event_generator::denPair,
			   HBT_event_generator::denPair2;
vector<double> HBT_event_generator::numerator_numPair,
			   HBT_event_generator::denominator_denPair;
vector<bool>   HBT_event_generator::denominator_cell_was_filled;
vector<int>    HBT_event_generator::numerator_bin_count,
			   HBT_event_generator::denominator_bin_count;


void HBT_event_generator::initialize_all(
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
	q_mode 			= paraRdr->getVal("q_mode");
	method_mode 	= paraRdr->getVal("method_mode");
	BE_mode 		= paraRdr->getVal("BE_mode");
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
	//   correlation function
	n_qo_pts 		= paraRdr->getVal("n_qo_pts");
	n_qs_pts 		= paraRdr->getVal("n_qs_pts");
	n_ql_pts 		= paraRdr->getVal("n_ql_pts");
	n_Q_pts 		= paraRdr->getVal("n_Q_pts");
	n_qRP_pts		= paraRdr->getVal("n_qRP_pts");
	n_thq_pts		= paraRdr->getVal("n_thq_pts");
	// - step size in q directions
	delta_qo 		= paraRdr->getVal("delta_qo");
	delta_qs 		= paraRdr->getVal("delta_qs");
	delta_ql 		= paraRdr->getVal("delta_ql");
	delta_Q 		= paraRdr->getVal("delta_Q");
	// - minimum value in each q direction
	qo_min 			= -0.5*double(n_qo_pts-1)*delta_qo;
	qs_min 			= -0.5*double(n_qs_pts-1)*delta_qs;
	ql_min 			= -0.5*double(n_ql_pts-1)*delta_ql;
	Q_min 			= -0.5*double(n_Q_pts-1)*delta_Q;
	qo_max 			= -qo_min;
	qs_max 			= -qs_min;
	ql_max 			= -ql_min;
	Q_max 			= -Q_min;

	// - number of points to use when fleshing out correlation
	//   function in each direction
	//new_nqopts 		= ( n_qo_pts > 1 ) ? new_nqpts : 1;
	//new_nqspts 		= ( n_qs_pts > 1 ) ? new_nqpts : 1;
	//new_nqlpts 		= ( n_ql_pts > 1 ) ? new_nqpts : 1;

	n_qo_bins 		= n_qo_pts - 1;
	n_qs_bins 		= n_qs_pts - 1;
	n_ql_bins 		= n_ql_pts - 1;
	n_Q_bins 		= n_Q_pts - 1;

	n_KT_bins 		= n_KT_pts - 1;
	n_Kphi_bins 	= n_Kphi_pts - 1;
	n_KL_bins 		= n_KL_pts - 1;

	KT_pts 			= vector<double> (n_KT_pts);
	Kphi_pts 		= vector<double> (n_Kphi_pts);
	KL_pts 			= vector<double> (n_KL_pts);

	qo_pts 			= vector<double> (n_qo_pts);
	qs_pts 			= vector<double> (n_qs_pts);
	ql_pts 			= vector<double> (n_ql_pts);
	Q_pts 			= vector<double> (n_Q_pts);

	x_pts			= vector<double> (n_qRP_pts);
	x_wts			= vector<double> (n_qRP_pts);
	ttheta_q_pts	= vector<double> (n_thq_pts);
	ttheta_q_wts	= vector<double> (n_thq_pts);

	linspace(KT_pts, KT_min, KT_max);
	linspace(Kphi_pts, Kphi_min, Kphi_max);
	linspace(KL_pts, KL_min, KL_max);

	linspace(qo_pts, qo_min, qo_max);
	linspace(qs_pts, qs_min, qs_max);
	linspace(ql_pts, ql_min, ql_max);
	linspace(Q_pts, Q_min, Q_max);

	//gauss_quadrature(n_qRP_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);
	gauss_quadrature(n_qRP_pts, 2, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);	//chebyshev distribution works better
	gauss_quadrature(n_thq_pts, 1, 0.0, 0.0, -M_PI, M_PI, ttheta_q_pts, ttheta_q_wts);

	px_bin_width 	= bin_epsilon;
	py_bin_width 	= bin_epsilon;
	pz_bin_width 	= bin_epsilon;

	// need to know these for binning particle pairs efficiently
	KT_bin_width 	= KT_pts[1]-KT_pts[0];
	Kphi_bin_width 	= Kphi_pts[1]-Kphi_pts[0];
	KL_bin_width 	= KL_pts[1]-KL_pts[0];

	const int q_space_size = ( q_mode == 0 ) ?
								n_qo_bins*n_qs_bins*n_ql_bins :
								n_Q_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	n_pair_numerator 		= 0.0;
	n_pair_denominator 		= 0.0;

	// For the correlation function itself
	numerator 				= vector<double> (K_space_size*q_space_size);
	denominator 			= vector<double> (K_space_size*q_space_size);
	correlation_function 	= vector<double> (K_space_size*q_space_size);
	correlation_function_error 
							= vector<double> (K_space_size*q_space_size);

	numerator2 				= vector<double> (K_space_size*q_space_size);
	denominator2 			= vector<double> (K_space_size*q_space_size);
	numerator_denominator 	= vector<double> (K_space_size*q_space_size);

	denominator_cell_was_filled
							= vector<bool> (K_space_size*q_space_size, false);
	numerator_bin_count		= vector<int> (K_space_size*q_space_size);
	denominator_bin_count	= vector<int> (K_space_size*q_space_size);

	numPair 				= vector<double> (K_space_size);
	numPair2 				= vector<double> (K_space_size);
	denPair 				= vector<double> (K_space_size);
	denPair2 				= vector<double> (K_space_size);

	numerator_numPair 		= vector<double> (K_space_size*q_space_size);
	denominator_denPair 	= vector<double> (K_space_size*q_space_size);

	// Initializations finished
	// Check number of events and proceed if non-zero
	if ( allEvents.size() == 0 )
		return;
	else
		out << "allEvents.size() = " << allEvents.size() << ": doing this file!" << endl;

	// Compute numerator and denominator of correlation function
	Compute_numerator_and_denominator();

	return;
}

HBT_event_generator::~HBT_event_generator()
{
	//clear everything

	return;
}



void HBT_event_generator::Compute_numerator_and_denominator()
{
	if ( BE_mode == 0 )
	{
		switch ( method_mode )
		{
			case 0:
				Compute_numerator_and_denominator_methodMode0();
				break;
			case 1:
				Compute_numerator_and_denominator_methodMode1();
				break;
			case 2:
				Compute_numerator_and_denominator_methodMode2();
				break;
			default:
				err << "HBT_event_generator(): method_mode = "
					<< method_mode << " not supported!" << endl;
				exit(8);
				break;
		}
	}
	else if ( BE_mode == 1 )
	{
		Compute_numerator_and_denominator_momentum_space_only();
	}
	else
	{
		err << "HBT_event_generator(): BE_mode = "
			<< BE_mode << " not supported!" << endl;
		exit(8);
	}

	return;
}



void HBT_event_generator::Compute_correlation_function()
{
	bool verbose = false;

	double nev = (double)total_N_events;
	const double prefactor = nev / ( nev - 1.0 );

	const int q_space_size = ( q_mode == 0 ) ?
								n_qo_bins*n_qs_bins*n_ql_bins :
								n_Q_bins;

	const int iqCenter = (q_space_size-1)/2;

	if ( method_mode == 1 or BE_mode == 1 )
	{
		// Compute correlation function itself
		// (along with error estimates)
		int idx = 0, idxK = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
			for (int iq = 0; iq < q_space_size; iq++)
			{

				double CF_num = numerator[idx] / (numPair[idxK]+1.e-100);
				double CF_den = denominator[idx] / (denPair[idxK]+1.e-100);

				double R2 = CF_num / (CF_den+1.e-100);

				bool this_bin_is_safe = numPair[idxK] > 0
										and denPair[idxK] > 0
										and abs( CF_den ) > 1.e-25;

				//==============================
				//==== correlation function ====
				//==============================
				correlation_function[idx]
							= ( this_bin_is_safe )
								? 1.0 + R2 - static_cast<double>(BE_mode)
								: 1.0;

				//=========================
				//==== error estimates ====
				//=========================

				if ( this_bin_is_safe and nev > 1 )
				{
					double CF_num_err = 0.0, CF_den_err = 0.0;

					// numerator
					CF_num_err
						= estimate_ratio_error(
							numerator[idx], numPair[idxK],
							numerator2[idx], numPair2[idxK], numerator_numPair[idx],
							nev, verbose, err );

					// denominator
					CF_den_err
						= estimate_ratio_error(
							denominator[idx], denPair[idxK],
							denominator2[idx], denPair2[idxK], denominator_denPair[idx],
							nev, verbose, err );

					// equivalent version which is stable if CF_num = 0
					correlation_function_error[idx] =
						sqrt(
								( CF_num_err * CF_num_err / ( CF_den * CF_den+1.e-100 ) )
									+ ( CF_num * CF_num * CF_den_err * CF_den_err
										/ ( CF_den * CF_den * CF_den * CF_den+1.e-100 )
									   )
							);
				}
				else
				{
					correlation_function_error[idx]
								= 1.0e6;	// maximal uncertainty?
				}

				++idx;
			}

			++idxK;
		}
	}
	else if ( method_mode == 0 or method_mode == 2 )
	{
		// Compute correlation function itself
		// (along with error estimates)
		int idx = 0;
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		for (int iq = 0; iq < q_space_size; iq++)
		{
			double num = numerator[idx];
			double den = denominator[idx];
			double num2 = numerator2[idx];
			double den2 = denominator2[idx];
			double numden = numerator_denominator[idx];

			double R2 = num / (den+1.e-100); // total_N_events factors cancel

			bool methodMode0_is_safe = (method_mode == 0) and (denominator_cell_was_filled[idx]) and (nev > 1);
			bool methodMode2_is_safe = (method_mode == 2) and (abs(den)> 1.e-25) and (nev > 1);
			bool this_bin_is_safe = methodMode0_is_safe or methodMode2_is_safe;

			//==============================
			//==== correlation function ====
			//==============================
			correlation_function[idx]
						= ( this_bin_is_safe )
							? 1.0 + R2
							: 1.0;

			//=========================
			//==== error estimates ====
			//=========================

			bool do_not_get_error_in_center = false;
			bool in_center_cell = ( iq == iqCenter );

			if ( in_center_cell and do_not_get_error_in_center )
			{
				correlation_function_error[idx]
							= 0.0;		// no error at origin, by definition,
										// unless we calculate it
			}
			else if ( this_bin_is_safe )
			{
				correlation_function_error[idx]
					= estimate_ratio_error(
						num, den,
						num2, den2, numden,
						nev, verbose, err );
			}
			else
			{
				correlation_function_error[idx]
							= 1.0e6;	// maximal uncertainty?
			}


			++idx;
		}
	}

	return;
}


void HBT_event_generator::Output_correlation_function( string filename )
{
	switch(q_mode)
	{
		case 0:
			Output_correlation_function_q_mode_3D( filename );
			break;
		case 1:
			Output_correlation_function_q_mode_1D( filename );
			break;
		default:
			err << "Output_correlation_function(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}



void HBT_event_generator::Output_correlation_function_q_mode_3D( string filename )
{
	int prec = 6;
	int extrawidth = 6;

	ofstream ofs( filename.c_str() );

	std::ios oldState(nullptr);
	oldState.copyfmt(ofs);

	// Print header inforamtion
	ofs /*<< setfill('X') */<< setw(prec+extrawidth+2)
		<< left << "# K_T" << setw(prec+extrawidth)
		<< left << "K_phi" << setw(prec+extrawidth)
		<< left << "K_L" << setw(prec+extrawidth)
		<< left << "q_o" << setw(prec+extrawidth)
		<< left << "q_s" << setw(prec+extrawidth)
		<< left << "q_l" << setw(prec+16)
		<< left << "re(N)" << setw(prec+16)
		<< left << "im(N)" << setw(prec+16)
		<< left << "D" << setw(prec+36)
		<< left << "C" << endl;

	ofs << "# " << setfill('-') << setw(150) << " " << endl;

	ofs.copyfmt(oldState);

	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iqo = 0; iqo < n_qo_bins; iqo++)
	for (int iqs = 0; iqs < n_qs_bins; iqs++)
	for (int iql = 0; iql < n_ql_bins; iql++)
	{

		ofs /*<< setfill('X') */<< fixed << setprecision(prec) << "  "
			<< 0.5*(KT_pts[iKT]+KT_pts[iKT+1]) << setw(prec+extrawidth)
			<< 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << setw(prec+extrawidth)
			<< 0.5*(KL_pts[iKL]+KL_pts[iKL+1]) << setw(prec+extrawidth)
			<< 0.5*(qo_pts[iqo]+qo_pts[iqo+1]) << setw(prec+extrawidth)
			<< 0.5*(qs_pts[iqs]+qs_pts[iqs+1]) << setw(prec+extrawidth)
			<< 0.5*(ql_pts[iql]+ql_pts[iql+1]) << setw(prec+16)
			//<< real(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			//<< imag(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			//<< denominator[idx] / static_cast<double>(total_N_events) << setw(prec+36)
			<< scientific
			<< real(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< imag(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< denominator[idx] / static_cast<double>(total_N_events) << setw(prec+36)
			<< fixed
			<< setprecision(16) << correlation_function[idx] << setw(prec+36)
			<< setprecision(16) << correlation_function_error[idx] << endl;

		++idx;
	}

	ofs.close();

	return;
}


void HBT_event_generator::Output_correlation_function_q_mode_1D( string filename )
{
	int prec = 6;
	int extrawidth = 6;

	ofstream ofs( filename.c_str() );

	std::ios oldState(nullptr);
	oldState.copyfmt(ofs);

	// Print header inforamtion
	ofs /*<< setfill('X') */<< setw(prec+extrawidth+2)
		<< left << "# K_T" << setw(prec+extrawidth)
		<< left << "K_phi" << setw(prec+extrawidth)
		<< left << "K_L" << setw(prec+extrawidth)
		<< left << "Q" << setw(prec+extrawidth)
		<< left << "re(N)" << setw(prec+extrawidth)
		<< left << "im(N)" << setw(prec+extrawidth)
		<< left << "D" << setw(prec+3*extrawidth)
		<< left << "C" << setw(prec+3*extrawidth)
		<< left << "C(err)" << endl;

	ofs << "# " << setfill('-') << setw(130) << " " << endl;

	ofs.copyfmt(oldState);

	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iQ = 0; iQ < n_Q_bins; iQ++)
	{

		ofs /*<< setfill('X') */<< fixed << setprecision(prec) << "  "
			<< 0.5*(KT_pts[iKT]+KT_pts[iKT+1]) << setw(prec+extrawidth)
			<< 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << setw(prec+extrawidth)
			<< 0.5*(KL_pts[iKL]+KL_pts[iKL+1]) << setw(prec+extrawidth)
			<< 0.5*(Q_pts[iQ]+Q_pts[iQ+1]) << "   "//setw(prec+24)
			<< real(numerator[idx] / static_cast<double>(total_N_events)) << "   "//setw(prec+16)
			<< imag(numerator[idx] / static_cast<double>(total_N_events)) << "   "//setw(prec+16)
			<< denominator[idx] / static_cast<double>(total_N_events) << "   "//setw(prec+24)
			<< setprecision(16) << correlation_function[idx] << "   "//setw(prec+24)
			<< setprecision(16) << correlation_function_error[idx] << endl;

		++idx;
	}

	ofs.close();

	return;
}



void HBT_event_generator::Update_records( const vector<EventRecord> & allEvents_in )
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

	// Compute numerator and denominator of correlation function
	Compute_numerator_and_denominator();

	return;
}


void HBT_event_generator::Compute_numerator_and_denominator_methodMode0()
{
	switch(q_mode)
	{
		case 0:
			Compute_numerator_and_denominator_methodMode0_q_mode_3D();
			break;
		case 1:
			Compute_numerator_and_denominator_methodMode0_q_mode_1D();
			break;
		default:
			err << "Compute_numerator_and_denominator_methodMode0(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}

void HBT_event_generator::Compute_numerator_and_denominator_methodMode1()
{
	switch(q_mode)
	{
		case 0:
			Compute_numerator_and_denominator_methodMode1_q_mode_3D();
			break;
		case 1:
			Compute_numerator_and_denominator_methodMode1_q_mode_1D();
			break;
		default:
			err << "Compute_numerator_and_denominator_methodMode0(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}


void HBT_event_generator::Compute_numerator_and_denominator_methodMode2()
{
	switch(q_mode)
	{
		case 0:
			Compute_numerator_and_denominator_methodMode2_q_mode_3D();
			break;
		case 1:
			Compute_numerator_and_denominator_methodMode2_q_mode_1D();
			break;
		default:
			err << "Compute_numerator_and_denominator_methodMode0(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}


void HBT_event_generator::Compute_numerator_and_denominator_momentum_space_only()
{
	switch(q_mode)
	{
		case 0:
			Compute_numerator_and_denominator_momentum_space_only_q_mode_3D();
			break;
		case 1:
			Compute_numerator_and_denominator_momentum_space_only_q_mode_1D();
			break;
		default:
			err << "Compute_numerator_and_denominator_momentum_space_only(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}


void HBT_event_generator::get_random_angles(int n_mixed_events, vector<double> & random_angles)
{
	random_device rd;
	mt19937_64 mt(rd());
	uniform_real_distribution<double> distribution(-M_PI,M_PI);

	random_angles.clear();
	for (int rand_idx = 0; rand_idx < n_mixed_events; ++rand_idx)
		random_angles.push_back( distribution(mt) );

	return;
}




//End of file
