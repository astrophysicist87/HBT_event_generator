#include <iostream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "HBT_event_generator.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

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
	// - some parameters
	bin_mode 		= paraRdr->getVal("bin_mode");
	q_mode 			= paraRdr->getVal("q_mode");
	//Define various grid sizes
	// - SP momentum points at which to evaluate correlation function
	n_pT_pts 		= paraRdr->getVal("n_pT_pts");
	pT_min 			= paraRdr->getVal("pTmin");
	pT_max 			= paraRdr->getVal("pTmax");
	n_pphi_pts 		= paraRdr->getVal("n_pphi_pts");
	pphi_min 		= -M_PI;
	pphi_max 		= M_PI;
	n_pY_pts 		= paraRdr->getVal("n_pY_pts");
	pY_min 			= paraRdr->getVal("pYmin");
	pY_max 			= paraRdr->getVal("pYmax");
	//
	n_px_pts 		= paraRdr->getVal("n_px_pts");
	n_py_pts 		= paraRdr->getVal("n_py_pts");
	n_pz_pts 		= paraRdr->getVal("n_pz_pts");
	px_min 			= paraRdr->getVal("pxmin");
	px_max 			= paraRdr->getVal("pxmax");
	py_min 			= paraRdr->getVal("pymin");
	py_max 			= paraRdr->getVal("pymax");
	pz_min 			= paraRdr->getVal("pzmin");
	pz_max 			= paraRdr->getVal("pzmax");
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
	n_Q_pts 		= paraRdr->getVal("n_Q_pts");
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

	n_pT_bins 		= n_pT_pts  - 1;
	n_pphi_bins 	= n_pphi_pts  - 1;
	n_pY_bins 		= n_pY_pts  - 1;

	n_KT_bins 		= n_KT_pts - 1;
	n_Kphi_bins 	= n_Kphi_pts - 1;
	n_KL_bins 		= n_KL_pts - 1;

	pT_pts 			= vector<double> (n_pT_pts);
	pphi_pts 		= vector<double> (n_pphi_pts);
	pY_pts 			= vector<double> (n_pY_pts);
	px_pts 			= vector<double> (n_px_pts);
	py_pts 			= vector<double> (n_py_pts);
	pz_pts 			= vector<double> (n_pz_pts);

	KT_pts 			= vector<double> (n_KT_pts);
	Kphi_pts 		= vector<double> (n_Kphi_pts);
	KL_pts 			= vector<double> (n_KL_pts);

	qo_pts 			= vector<double> (n_qo_pts);
	qs_pts 			= vector<double> (n_qs_pts);
	ql_pts 			= vector<double> (n_ql_pts);
	Q_pts 			= vector<double> (n_Q_pts);

	dN_pTdpTdpphidpY = vector<double> (n_pT_bins*n_pphi_bins*n_pY_pts);

	linspace(pT_pts, pT_min, pT_max);
	linspace(pphi_pts, pphi_min, pphi_max);
	linspace(pY_pts, pY_min, pY_max);
	linspace(px_pts, px_min, px_max);
	linspace(py_pts, py_min, py_max);
	linspace(pz_pts, pz_min, pz_max);

	linspace(KT_pts, KT_min, KT_max);
	linspace(Kphi_pts, Kphi_min, Kphi_max);
	linspace(KL_pts, KL_min, KL_max);

	linspace(qo_pts, qo_min, qo_max);
	linspace(qs_pts, qs_min, qs_max);
	linspace(ql_pts, ql_min, ql_max);
	linspace(Q_pts, Q_min, Q_max);

	pT_bin_width 	= pT_pts[1]-pT_pts[0];
	pphi_bin_width 	= pphi_pts[1]-pphi_pts[0];
	pY_bin_width 	= pY_pts[1]-pY_pts[0];
	//px_bin_width 	= px_pts[1]-px_pts[0];
	//py_bin_width 	= py_pts[1]-py_pts[0];
	//pz_bin_width 	= pz_pts[1]-pz_pts[0];
	px_bin_width 	= paraRdr->getVal("n_px_pts");
	py_bin_width 	= paraRdr->getVal("n_py_pts");
	pz_bin_width 	= paraRdr->getVal("n_pz_pts");

	// assume uniform grid spacing for now
	KT_bin_width 	= KT_pts[1]-KT_pts[0];
	Kphi_bin_width 	= Kphi_pts[1]-Kphi_pts[0];
	KL_bin_width 	= KL_pts[1]-KL_pts[0];

	const int q_space_size = ( q_mode == 0 ) ?
								n_qo_bins*n_qs_bins*n_ql_bins :
								n_Q_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

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

	// Compute numerator and denominator of correlation function,
	// along with quantities needed to estimate error
	Compute_numerator_and_denominator_with_errors(
							numerator,
							numerator2,
							denominator,
							denominator2,
							numerator_denominator,
							numerator_bin_count,
							denominator_bin_count
					);

	return;
}

HBT_event_generator::~HBT_event_generator()
{
	//clear everything

	return;
}


void HBT_event_generator::Compute_correlation_function()
{
	switch(q_mode)
	{
		case 0:
			Compute_correlation_function_q_mode_3D();
			break;
		case 1:
			Compute_correlation_function_q_mode_1D();
			break;
		default:
			err << "Compute_correlation_function(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;
}


void HBT_event_generator::Compute_correlation_function_q_mode_3D()
{
	const double prefactor
					= static_cast<double>(total_N_events)
						/ static_cast<double>(total_N_events-1);
	const int iqoC = (n_qo_pts - 1) / 2;
	const int iqsC = (n_qs_pts - 1) / 2;
	const int iqlC = (n_ql_pts - 1) / 2;

	// Compute correlation function itself
	// (along with error estimates)
	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iqo = 0; iqo < n_qo_bins; iqo++)
	for (int iqs = 0; iqs < n_qs_bins; iqs++)
	for (int iql = 0; iql < n_ql_bins; iql++)
	{
		err << "Computing correlation function, loop: "
			<< iKT << "   " << iKphi << "   " << iKL << "   "
			<< iqo << "   " << iqs << "   " << iql << endl;

		double num = numerator[idx];
		double den = denominator[idx];
		double num2 = numerator2[idx];
		double den2 = denominator2[idx];
		double numden = numerator_denominator[idx];

		double R2 = num / den; // total_N_events factors cancel

		//==============================
		//==== correlation function ====
		//==============================
		correlation_function[idx]
					= ( denominator_cell_was_filled[idx] )
						? 1.0 + R2
						: 1.0;

		//=========================
		//==== error estimates ====
		//=========================

		bool do_not_get_error_in_center = false;
		bool in_center_cell = ( iqo == iqoC
								and iqs == iqsC
								and iql == iqlC );

		if ( in_center_cell and do_not_get_error_in_center )
		{
			correlation_function_error[idx]
						= 0.0;		// no error at origin, by definition,
									// unless we calculate it
		}
		else if ( denominator_cell_was_filled[idx] )
		{
bool verbose = (idx==1) or (idx==9);
verbose = false;

if (verbose) err << setprecision(8);
			// average of numerator in this q-K cell
			double EA = num / static_cast<double>(total_N_events);
if (verbose) err << "\t\t EA = " << EA << endl;
			// average of denominator in this q-K cell
			double EB = den / static_cast<double>(total_N_events);
if (verbose) err << "\t\t EB = " << EB << endl;
			// average of numerator^2 in this q-K cell
			double EA2 = num2 / static_cast<double>(total_N_events);
if (verbose) err << "\t\t EA2 = " << EA2 << endl;
			// average of denominator^2 in this q-K cell
			double EB2 = den2 / static_cast<double>(total_N_events);
if (verbose) err << "\t\t EB2 = " << EB2 << endl;
			// average of numerator*denominator in this q-K cell
			double EAB = numden / static_cast<double>(total_N_events);
if (verbose) err << "\t\t EAB = " << EAB << endl;
			// set variances and covariance
			double sigA2 = prefactor * (EA2 - EA*EA);
			double sigB2 = prefactor * (EB2 - EB*EB);
			double sigAB = prefactor * (EAB - EA*EB);
if (verbose) err << "\t\t sigA2 = " << sigA2 << endl;
if (verbose) err << "\t\t sigB2 = " << sigB2 << endl;
if (verbose) err << "\t\t sigAB = " << sigAB << endl;

			// want standard error, not variance itself
			sigA2 /= static_cast<double>(total_N_events);
			sigB2 /= static_cast<double>(total_N_events);
			sigAB /= static_cast<double>(total_N_events);
if (verbose) err << "\t\t sigA2 / total_N_events = " << sigA2 << endl;
if (verbose) err << "\t\t sigB2 / total_N_events = " << sigB2 << endl;
if (verbose) err << "\t\t sigAB / total_N_events = " << sigAB << endl;

			// set relative widths
			double cA = sigA2 / ( EA*EA+1.e-100 );
			double cB = sigB2 / ( EB*EB+1.e-100 );
			double cAB = sigAB / ( EA*EB+1.e-100 );
if (verbose) err << "\t\t cA = " << cA << endl;
if (verbose) err << "\t\t cB = " << cB << endl;
if (verbose) err << "\t\t cAB = " << cAB << endl;

			double disc = cA + cB - 2.0*cAB;
if (verbose) err << "\t\t disc = " << disc << endl;
if (verbose and disc < -1.e-6) err << "disc < 0!" << endl;

			// N.B.: R2 == EA / EB

			correlation_function_error[idx] =
				( disc < 0.0 )
				? 1.e-6
				: abs(R2) * sqrt( cA + cB - 2.0*cAB );
if (verbose) err << "\t\t Finally, check bin counts:" << endl;
if (verbose) err << "\t\t num. BC = " << numerator_bin_count[idx] << endl;
if (verbose) err << "\t\t den. BC = " << denominator_bin_count[idx] << endl;
if (verbose) err << "\t\t finished without difficulty" << endl;
		}
		/*else if ( iqo == iqoC
					and iqs == iqsC
					and iql == iqlC )
		{
			correlation_function_error[idx]
						= 0.0;		// no error at origin, by definition
		}*/
		else
		{
			correlation_function_error[idx]
						= 1.0e6;	// maximal uncertainty?
		}


		++idx;
	}

	return;
}



void HBT_event_generator::Compute_correlation_function_q_mode_1D()
{
	bool verbose = true;

	const double prefactor
					= static_cast<double>(total_N_events)
						/ static_cast<double>(total_N_events-1);
	const int iQC = (n_Q_pts - 1) / 2;

	// Compute correlation function itself
	// (along with error estimates)
	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iQ = 0; iQ < n_Q_bins; iQ++)
	{
		err << "Computing correlation function, loop: "
			<< iKT << "   " << iKphi << "   " << iKL << "   "
			<< iQ << endl;

		double num = numerator[idx];
		double den = denominator[idx];
		double num2 = numerator2[idx];
		double den2 = denominator2[idx];
		double numden = numerator_denominator[idx];

		double R2 = num / den; // total_N_events factors cancel

		//==============================
		//==== correlation function ====
		//==============================
		correlation_function[idx]
					= ( denominator_cell_was_filled[idx] )
						? 1.0 + R2
						: 1.0;

		//=========================
		//==== error estimates ====
		//=========================

		bool do_not_get_error_in_center = false;
		bool in_center_cell = ( iQ == iQC );

		if ( in_center_cell and do_not_get_error_in_center )
		{
			correlation_function_error[idx]
						= 0.0;		// no error at origin, by definition,
									// unless we calculate it
		}
		else if ( denominator_cell_was_filled[idx] )
		{

			double EA = num / static_cast<double>(total_N_events);
			double EB = den / static_cast<double>(total_N_events);

			double EA2 = num2 / static_cast<double>(total_N_events);
			double EB2 = den2 / static_cast<double>(total_N_events);
			double EAB = numden / static_cast<double>(total_N_events);

			double sigA2 = prefactor * (EA2 - EA*EA);
			double sigB2 = prefactor * (EB2 - EB*EB);
			double sigAB = prefactor * (EAB - EA*EB);

			sigA2 /= static_cast<double>(total_N_events);
			sigB2 /= static_cast<double>(total_N_events);
			sigAB /= static_cast<double>(total_N_events);

			double cA = sigA2 / ( EA*EA+1.e-100 );
			double cB = sigB2 / ( EB*EB+1.e-100 );
			double cAB = sigAB / ( EA*EB+1.e-100 );

			double disc = cA + cB - 2.0*cAB;
			if (verbose and disc < -1.e-6)
				err << "disc = " << disc << " < 0!" << endl;

			correlation_function_error[idx] =
				( disc < 0.0 )
				? 1.e-6
				: abs(R2) * sqrt( cA + cB - 2.0*cAB );
		}
		else
		{
			correlation_function_error[idx]
						= 1.0e6;	// maximal uncertainty?
		}

		++idx;
	}

	return;
}



void HBT_event_generator::Output_correlation_function()
{
	switch(q_mode)
	{
		case 0:
			Output_correlation_function_q_mode_3D();
			break;
		case 1:
			Output_correlation_function_q_mode_1D();
			break;
		default:
			err << "Output_correlation_function(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}



void HBT_event_generator::Output_correlation_function_q_mode_3D()
{
	int prec = 6;
	int extrawidth = 6;

	std::ios oldState(nullptr);
	oldState.copyfmt(out);

	// Print header inforamtion
	out /*<< setfill('X') */<< setw(prec+extrawidth+2)
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

	out << "# " << setfill('-') << setw(150) << " " << endl;

	out.copyfmt(oldState);

	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iqo = 0; iqo < n_qo_bins; iqo++)
	for (int iqs = 0; iqs < n_qs_bins; iqs++)
	for (int iql = 0; iql < n_ql_bins; iql++)
	{

		out /*<< setfill('X') */<< fixed << setprecision(prec) << "  "
			<< 0.5*(KT_pts[iKT]+KT_pts[iKT+1]) << setw(prec+extrawidth)
			<< 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << setw(prec+extrawidth)
			<< 0.5*(KL_pts[iKL]+KL_pts[iKL+1]) << setw(prec+extrawidth)
			<< 0.5*(qo_pts[iqo]+qo_pts[iqo+1]) << setw(prec+extrawidth)
			<< 0.5*(qs_pts[iqs]+qs_pts[iqs+1]) << setw(prec+extrawidth)
			<< 0.5*(ql_pts[iql]+ql_pts[iql+1]) << setw(prec+16)
			<< real(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< imag(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< denominator[idx] / static_cast<double>(total_N_events) << setw(prec+36)
			<< setprecision(16) << correlation_function[idx] << setw(prec+36)
			<< setprecision(16) << correlation_function_error[idx] << endl;

		++idx;
	}

	return;
}


void HBT_event_generator::Output_correlation_function_q_mode_1D()
{
	int prec = 6;
	int extrawidth = 6;

	std::ios oldState(nullptr);
	oldState.copyfmt(out);

	// Print header inforamtion
	out /*<< setfill('X') */<< setw(prec+extrawidth+2)
		<< left << "# K_T" << setw(prec+extrawidth)
		<< left << "K_phi" << setw(prec+extrawidth)
		<< left << "K_L" << setw(prec+extrawidth)
		<< left << "Q" << setw(prec+16)
		<< left << "re(N)" << setw(prec+16)
		<< left << "im(N)" << setw(prec+16)
		<< left << "D" << setw(prec+36)
		<< left << "C" << endl;

	out << "# " << setfill('-') << setw(150) << " " << endl;

	out.copyfmt(oldState);

	int idx = 0;
	for (int iKT = 0; iKT < n_KT_bins; iKT++)
	for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
	for (int iKL = 0; iKL < n_KL_bins; iKL++)
	for (int iQ = 0; iQ < n_Q_bins; iQ++)
	{

		out /*<< setfill('X') */<< fixed << setprecision(prec) << "  "
			<< 0.5*(KT_pts[iKT]+KT_pts[iKT+1]) << setw(prec+extrawidth)
			<< 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << setw(prec+extrawidth)
			<< 0.5*(KL_pts[iKL]+KL_pts[iKL+1]) << setw(prec+extrawidth)
			<< 0.5*(Q_pts[iQ]+Q_pts[iQ+1]) << setw(prec+16)
			<< real(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< imag(numerator[idx] / static_cast<double>(total_N_events)) << setw(prec+16)
			<< denominator[idx] / static_cast<double>(total_N_events) << setw(prec+24)
			<< setprecision(16) << correlation_function[idx] << setw(prec+24)
			<< setprecision(16) << correlation_function_error[idx] << endl;

		++idx;
	}

	return;
}



void HBT_event_generator::Update_records( const vector<EventRecord> & allEvents_in )
{

	// Copy in new records of all events
	// (erases old event information)
	allEvents		= allEvents_in;
	total_N_events	+= allEvents.size();

	// Compute numerator and denominator of correlation function,
	// along with quantities needed to estimate error
	Compute_numerator_and_denominator_with_errors(
							numerator,
							numerator2,
							denominator,
							denominator2,
							numerator_denominator,
							numerator_bin_count,
							denominator_bin_count
					);


	return;
}


void HBT_event_generator::Compute_numerator_and_denominator_with_errors(
							vector<double> & in_numerator,
							vector<double> & in_numerator2,
							vector<double> & in_denominator,
							vector<double> & in_denominator2,
							vector<double> & in_numerator_denominator,
							vector<int> & in_numerator_bin_count,
							vector<int> & in_denominator_bin_count
							)
{
	switch(q_mode)
	{
		case 0:
			Compute_numerator_and_denominator_with_errors_q_mode_3D(
							numerator, numerator2, denominator, denominator2,
							numerator_denominator, numerator_bin_count, denominator_bin_count
							);
			break;
		case 1:
			Compute_numerator_and_denominator_with_errors_q_mode_1D(
							numerator, numerator2, denominator, denominator2,
							numerator_denominator, numerator_bin_count, denominator_bin_count
							);
			break;
		default:
			err << "Compute_numerator_and_denominator_with_errors(): q_mode = "
				<< q_mode << " not supported!" << endl;
			exit(8);
			break;
	}

	return;	
}



//End of file
