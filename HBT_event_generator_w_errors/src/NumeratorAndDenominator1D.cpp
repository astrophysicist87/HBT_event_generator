#include <iostream>
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
#include "Arsenal.h"
#include "Stopwatch.h"


void HBT_event_generator::Compute_numerator_and_denominator_with_errors_q_mode_1D(
							vector<double> & in_numerator, 				vector<double> & in_numerator2,
							vector<double> & in_denominator, 			vector<double> & in_denominator2,
							vector<double> & in_numerator_denominator, 	vector<int>    & in_numerator_bin_count,
							vector<int>    & in_denominator_bin_count
							)
{
	const int iqoC = (n_qo_pts - 1) / 2;
	const int iqsC = (n_qs_pts - 1) / 2;
	const int iqlC = (n_ql_pts - 1) / 2;

	//int number_of_completed_events = 0;
	//err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	const int q_space_size = n_qo_bins*n_qs_bins*n_ql_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// check multiplicities
	//for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	//	err << "allEvents[" << iEvent << "].particles.size() = " << allEvents[iEvent].particles.size() << endl;

	// Sum over all events
	#pragma omp parallel for schedule(static) shared( in_numerator, in_numerator2,\
														in_denominator, in_denominator2,\
														in_numerator_denominator,\
														in_numerator_bin_count,\
														in_denominator_bin_count )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		//let these be fully six-dimensional
		vector<complex<double> > sum1(K_space_size*q_space_size);
		vector<double> sum2(K_space_size*q_space_size);
		vector<double> sum3(K_space_size*q_space_size);
		vector<double> sum4(K_space_size*q_space_size);
		vector<double> sum5(K_space_size*q_space_size);

		vector<int> private_ABC(K_space_size*q_space_size);
		vector<int> private_DBC(K_space_size*q_space_size);

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
			double KTmin = KT_pts[iKT];
			double KTmax = KT_pts[iKT+1];
			double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
			double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double Kx = KT * cKphi;
			double Ky = KT * sKphi;
			double Kz = KL;

			// Sum over particles
			for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
			{
				ParticleRecord pi = event.particles[iParticle];

				double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;


				/*#pragma omp critical
                {
                        //if (iqo == iqoC and iqs == iqsC and iql == iqlC)
                                err << "CHECKPOINT #0: "
                                        << Kx << "   " << Ky << "   " << Kz << "   "
                                        << pi.px << "   " << pi.py << "   " << pi.pz << "   "
                                        << Kx - pi.px << "   "
					<< Ky - pi.py << "   "
					<< Kz - pi.pz << endl;
                }*/


				bool num_bin_true = 	abs( Kx - pi.px ) <= 0.5*px_bin_width
									and abs( Ky - pi.py ) <= 0.5*py_bin_width
									and abs( Kz - pi.pz ) <= 0.5*pz_bin_width;

				if ( num_bin_true )
				{
					double num_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;

					for (int iqo = 0; iqo < n_qo_bins; iqo++)
					for (int iqs = 0; iqs < n_qs_bins; iqs++)
					for (int iql = 0; iql < n_ql_bins; iql++)
					{
						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
						double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
						double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
						double qx = qo * cKphi - qs * sKphi;
						double qy = qs * cKphi + qo * sKphi;
						double qz = ql;

						double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
						double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;
						double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
						double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);

						// rapidity cuts
						if ( impose_pair_rapidity_cuts
								and ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min )
								or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
							)
							continue;

						double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

						double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

						complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;

						sum1[index6D] += complex_num_term;
						sum2[index6D] += 1.0 / (num_bin_factor*num_bin_factor);
						private_ABC[index6D]++;

		/*#pragma omp critical
		{
			if (iqo == iqoC and iqs == iqsC and iql == iqlC)
				err << "CHECKPOINT #1: "
					<< KT << "   " << 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << "   " << KL << "   " 
					<< qo << "   " << qs << "   " << ql << "   " 
					<< q0 << "   " << index6D << "   "
					<< complex_num_term << "   " << 1.0 / (num_bin_factor*num_bin_factor) << "   "
					<< sum1[index6D] << "   " << sum2[index6D] << endl;
		}*/


					}	// end of q loops

				}		// end of if block

			}			// end of particle loop

		}				// end of K loops




		//=====================================
		//======== Doing denominator ==========
		//=====================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
			double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
			double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double Kx = KT * cKphi;
			double Ky = KT * sKphi;
			double Kz = KL;

			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
				double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;
				double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
				double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);

				//===================================
				// checking Qbar cell in denominator
				int index4D_pos = -1, index4D_neg = -1;	//make sure they get set
				{
					double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);
					double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
					if (Q2bar < -1.e-6)
					{
						err << "Compute_numerator_and_denominator_with_errors_q_mode_1D(warning): "
							<< "Q2bar = " << Q2bar << " < 0.0!" << endl;
						continue;
					}
					double Qbar = sqrt(abs(Q2bar));

					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					// N.B. - these vectors have different dimensions!!!
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					int iQ_pos = floor((Qbar - Q_min)/delta_Q);
					if ( iQ_pos < 0 or iQ_pos >= n_Q_bins )
						continue;
					index4D_pos = indexer_qmode_1(iKT, iKphi, iKL, iQ_pos);

					// note that both fail or neither fails, by symmetry
					Qbar *= -1.0;
					int iQ_neg = floor((Qbar - Q_min)/delta_Q);
					if ( iQ_neg < 0 or iQ_neg >= n_Q_bins )
						continue;
					index4D_neg = indexer_qmode_1(iKT, iKphi, iKL, iQ_neg);
				}
				//===================================

				// rapidity cuts
				if ( impose_pair_rapidity_cuts
						and ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min )
						or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
					)
					continue;

				// Sum over particles for pa bin
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pax) <= 0.5*px_bin_width
										and abs(pi.py-pay) <= 0.5*py_bin_width
										and abs(pi.pz-paz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;

						denominator_cell_was_filled[index4D_pos] = true;
						denominator_cell_was_filled[index4D_neg] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum3[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum3 particle loop

				// Sum over particles for pb bin
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pbx) <= 0.5*px_bin_width
										and abs(pi.py-pby) <= 0.5*py_bin_width
										and abs(pi.pz-pbz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;

						//denominator_cell_was_filled[index4D_pos] = true;
						//denominator_cell_was_filled[index4D_neg] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum4[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum4 particle loop

				// Sum over particles for pa and pb bins simultaneously
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pax) <= 0.5*px_bin_width
										and abs(pi.py-pay) <= 0.5*py_bin_width
										and abs(pi.pz-paz) <= 0.5*pz_bin_width
										and abs(pi.px-pbx) <= 0.5*px_bin_width
										and abs(pi.py-pby) <= 0.5*py_bin_width
										and abs(pi.pz-pbz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						// N.B. - use (bin volume)^2, here
						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width
							*px_bin_width*py_bin_width*pz_bin_width;

						//denominator_cell_was_filled[index4D_pos] = true;
						//denominator_cell_was_filled[index4D_neg] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum5[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum5 particle loop

			}			// end of q loops

		}				// end of K loops

		// Need this to avoid race conditions
		#pragma omp critical
		{
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				int idx = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
				double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);
				double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
				if (Q2bar < -1.e-6)
				{
					err << "Compute_numerator_and_denominator_with_errors_q_mode_1D(warning): "
						<< "Q2bar = " << Q2bar << " < 0.0!" << endl;
					continue;
				}
				double Qbar = sqrt(abs(Q2bar));

				double abs_sum1 = abs(sum1[idx]);
				double numerator_contribution_from_this_event
						= abs_sum1*abs_sum1 - sum2[idx];
				double denominator_contribution_from_this_event
						= sum3[idx]*sum4[idx] - sum5[idx];

				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// N.B. - these vectors have different dimensions!!!
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				int iQ_pos = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_pos < 0 or iQ_pos >= n_Q_bins )
					continue;
				int index4D_pos = indexer_qmode_1(iKT, iKphi, iKL, iQ_pos);

if (iqo == iqoC and iqs == iqsC and iql == iqlC)
	err << KT << "   " << 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]) << "   " << KL << "   " 
		<< qo << "   " << qs << "   " << ql << "   " 
		<< q0 << "   " << Qbar << "   " << iQ_pos << "   " << index4D_pos << "   "
		<< abs_sum1 << "   " << sum2[idx] << "   " << sum3[idx] << "   " << sum4[idx] << "   " << sum5[idx] << "   "
		<< numerator_contribution_from_this_event << "   "
		<< denominator_contribution_from_this_event << endl;

				// first moments
				in_numerator[index4D_pos]
					+= numerator_contribution_from_this_event;
				in_denominator[index4D_pos]
					+= denominator_contribution_from_this_event;

				// second moments
				in_numerator2[index4D_pos]
					+= numerator_contribution_from_this_event
						* numerator_contribution_from_this_event;
				in_denominator2[index4D_pos]
					+= denominator_contribution_from_this_event
						* denominator_contribution_from_this_event;
				in_numerator_denominator[index4D_pos]
					+= numerator_contribution_from_this_event
						* denominator_contribution_from_this_event;

				// track number of events where bin count was non-vanishing
				in_numerator_bin_count[index4D_pos] += int(private_ABC[idx] > 0);
				in_denominator_bin_count[index4D_pos] += int(private_DBC[idx] > 0);

				// don't count zero point twice
				//if ( Qbar < 1.e-6 )
				//	continue;

				Qbar *= -1.0;
				int iQ_neg = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_neg < 0 or iQ_neg >= n_Q_bins )
					continue;
				int index4D_neg = indexer_qmode_1(iKT, iKphi, iKL, iQ_neg);

				// first moments
				in_numerator[index4D_neg]
					+= numerator_contribution_from_this_event;
				in_denominator[index4D_neg]
					+= denominator_contribution_from_this_event;

				// second moments
				in_numerator2[index4D_neg]
					+= numerator_contribution_from_this_event
						* numerator_contribution_from_this_event;
				in_denominator2[index4D_neg]
					+= denominator_contribution_from_this_event
						* denominator_contribution_from_this_event;
				in_numerator_denominator[index4D_neg]
					+= numerator_contribution_from_this_event
						* denominator_contribution_from_this_event;

				// track number of events where bin count was non-vanishing
				in_numerator_bin_count[index4D_neg] += int(private_ABC[idx] > 0);
				in_denominator_bin_count[index4D_neg] += int(private_DBC[idx] > 0);

			}

			//err << "\t - finished " << ++number_of_completed_events << " of " << total_N_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	//err << "  * Finished!" << endl;

	return;
}


//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
void HBT_event_generator::Compute_numerator_and_denominator_with_errors_q_mode_1D_momentum_space_only(
							vector<double> & in_numerator, vector<double> & in_numerator2,
							vector<double> & in_numPair, vector<double> & in_numPair2,
							vector<double> & in_denominator, vector<double> & in_denominator2,
							vector<double> & in_denPair, vector<double> & in_denPair2,
							vector<double> & in_numerator_numPair, vector<double> & in_denominator_denPair
							)
{
	bool perform_random_rotation = false;
	bool perform_random_shuffle = false;

	int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	const int q_space_size = n_qo_bins*n_qs_bins*n_ql_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// Sum over all events
	#pragma omp parallel for schedule(static) shared( in_numerator, in_numerator2,\
														in_numPair, in_numPair2,\
														in_denominator, in_denominator2,\
														in_denPair, in_denPair2,\
														in_numerator_numPair, in_denominator_denPair )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_num(K_space_size*q_space_size, 0.0);
		vector<double> private_num2(K_space_size*q_space_size, 0.0);
		vector<double> private_numPair(K_space_size, 0.0);
		vector<double> private_numPair2(K_space_size, 0.0);
		vector<double> private_den(K_space_size*q_space_size, 0.0);
		vector<double> private_den2(K_space_size*q_space_size, 0.0);
		vector<double> private_denPair(K_space_size, 0.0);
		vector<double> private_denPair2(K_space_size, 0.0);

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			if (iParticle == jParticle)
				continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			//double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
			//double tj = pj.t, xj = pj.x, yj = pj.y, zj = pj.z;
			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);
			double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

			// rapidity cuts
			if ( impose_pair_rapidity_cuts
					and ( ( 2.0*Kz/(Ei+Ej) < Kz_over_K0_min )
					or ( 2.0*Kz/(Ei+Ej) > Kz_over_K0_max ) )
				)
				continue;

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);

			double qo = qx * cKphi + qy * sKphi;
			double qs = qy * cKphi - qx * sKphi;
			double ql = qz;

			// Get indices
			int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
			int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

			int qo_idx 	= floor((qo - qo_min) / delta_qo);
			int qs_idx 	= floor((qs - qs_min) / delta_qs);
			int ql_idx 	= floor((ql - ql_min) / delta_ql);

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
				continue;

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				continue;

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
				continue;

			if ( qo_idx < 0 or qo_idx >= n_qo_bins )
				continue;

			if ( qs_idx < 0 or qs_idx >= n_qs_bins )
				continue;

			if ( ql_idx < 0 or ql_idx >= n_ql_bins )
				continue;

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
			int index6D = indexer(KT_idx, Kphi_idx, KL_idx, qo_idx, qs_idx, ql_idx);

			//double arg =  q0 * (ti - tj)
			//			- qx * (xi - xj)
			//			- qy * (yi - yj)
			//			- qz * (zi - zj);

			//double num_term = cos(arg/hbarC);
			double num_term = 1.0;

			private_num[index6D] += num_term;
			private_numPair[index3D]++;

		}


		//=====================================
		//========= Doing denominator =========
		//=====================================

		// Randomly sample events to mix with
		const unsigned int n_mixing_events = allEvents.size()-1;
		//const unsigned int n_mixing_events = 100;

		vector<unsigned int> indices(allEvents.size());
		iota(indices.begin(), indices.end(), 0);

		if ( perform_random_shuffle )
			random_shuffle(indices.begin(), indices.end());
		vector<int> mixedEvents;
		for (int mix_idx = 0; mix_idx <= n_mixing_events; ++mix_idx)
			if ( indices[mix_idx] != iEvent
					and mixedEvents.size() < n_mixing_events )
			{
				mixedEvents.push_back( indices[mix_idx] );
			}

		vector<double> random_angles(n_mixing_events, 0.0), cos_rand_angles, sin_rand_angles;
		if ( perform_random_rotation )
			get_random_angles(n_mixing_events, random_angles);
		for (int mix_idx = 0; mix_idx < random_angles.size(); ++mix_idx)
		{
			cos_rand_angles.push_back( cos( random_angles[mix_idx] ) );
			sin_rand_angles.push_back( sin( random_angles[mix_idx] ) );
		}

		// Loop over mixed events
		for (int jEvent = 0; jEvent < mixedEvents.size(); ++jEvent)
		{
			EventRecord mixedEvent = allEvents[mixedEvents[jEvent]];

			double c_rand_phi = cos_rand_angles[jEvent],
					s_rand_phi = sin_rand_angles[jEvent];

			//c_rand_phi = 1.0, s_rand_phi = 0.0;

			// Sum over pairs of particles
			for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
			for (int jParticle = 0; jParticle < mixedEvent.particles.size(); ++jParticle)
			{

				ParticleRecord pi = event.particles[iParticle];
				ParticleRecord pj = mixedEvent.particles[jParticle];

				double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;

				// Randomly rotate the mixed event in transverse plane
				double Ej = pj.E,
						pjx = pj.px*c_rand_phi - pj.py*s_rand_phi,
						pjy = pj.px*s_rand_phi + pj.py*c_rand_phi,
						pjz = pj.pz;

				// New method of binning
				double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);
				double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

				// rapidity cuts
				if ( impose_pair_rapidity_cuts
						and ( ( 2.0*Kz/(Ei+Ej) < Kz_over_K0_min )
						or ( 2.0*Kz/(Ei+Ej) > Kz_over_K0_max ) )
					)
					continue;

				double KT = sqrt(Kx*Kx+Ky*Ky);
				double Kphi = atan2(Ky, Kx);
				double cKphi = cos(Kphi), sKphi = sin(Kphi);

				double qo = qx * cKphi + qy * sKphi;
				double qs = qy * cKphi - qx * sKphi;
				double ql = qz;
	
				// If pair survived cuts, get indices
				int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
				int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
				int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

				int qo_idx 	= floor((qo - qo_min) / delta_qo);
				int qs_idx 	= floor((qs - qs_min) / delta_qs);
				int ql_idx 	= floor((ql - ql_min) / delta_ql);

				// Momentum-space cuts
				if ( KT_idx < 0 or KT_idx >= n_KT_bins )
					continue;

				if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
					continue;

				if ( KL_idx < 0 or KL_idx >= n_KL_bins )
					continue;

				if ( qo_idx < 0 or qo_idx >= n_qo_bins )
					continue;

				if ( qs_idx < 0 or qs_idx >= n_qs_bins )
					continue;

				if ( ql_idx < 0 or ql_idx >= n_ql_bins )
					continue;

				int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, qo_idx, qs_idx, ql_idx);

				private_den[index6D]++;
				private_denPair[index3D]++;

			}
		}


		// Need this to avoid race conditions
		#pragma omp critical
		{
			int idx3D = 0, idx6D = 0;
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			{

				double private_numPair_val 	= private_numPair[idx3D];
				double private_denPair_val 	= private_denPair[idx3D];

				in_numPair[idx3D] 			+= private_numPair_val;
				in_numPair2[idx3D] 			+= private_numPair_val*private_numPair_val;
				in_denPair[idx3D] 			+= private_denPair_val;
				in_denPair2[idx3D] 			+= private_denPair_val*private_denPair_val;

				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);

				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
					double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
					double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
					double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

					// set Q2 value for this q-K cell
					double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
					if (Q2bar < -1.e-6)
					{
						err << "Compute_numerator_and_denominator_with_errors_q_mode_1D(warning): "
							<< "Q2bar = " << Q2bar << " < 0.0!" << endl;
						continue;
					}
					double Qbar = sqrt(abs(Q2bar));

					double private_num_val 			= private_num[idx6D];
					double private_den_val 			= private_den[idx6D];

					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					// N.B. - these vectors have different dimensions!!!
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					int iQ_pos = floor((Qbar - Q_min)/delta_Q);
					if ( iQ_pos < 0 or iQ_pos >= n_Q_bins )
						continue;
					int index4D_pos = indexer_qmode_1(iKT, iKphi, iKL, iQ_pos);

					// assume symmetry under Q --> -Q
					Qbar *= -1.0;
					int iQ_neg = floor((Qbar - Q_min)/delta_Q);
					if ( iQ_neg < 0 or iQ_neg >= n_Q_bins )
						continue;
					int index4D_neg = indexer_qmode_1(iKT, iKphi, iKL, iQ_neg);

					// first moments
					in_numerator[index4D_pos] 				+= private_num_val;
					in_denominator[index4D_pos] 			+= private_den_val;
					in_numerator[index4D_neg] 				+= private_num_val;
					in_denominator[index4D_neg] 			+= private_den_val;

					// second moments
					in_numerator2[index4D_pos] 				+= private_num_val * private_num_val;
					in_denominator2[index4D_pos] 			+= private_den_val * private_den_val;
					in_numerator2[index4D_neg] 				+= private_num_val * private_num_val;
					in_denominator2[index4D_neg] 			+= private_den_val * private_den_val;

					// for error normalizing by numbers of pairs
					in_numerator_numPair[index4D_pos] 		+= private_num_val*private_numPair_val;
					in_denominator_denPair[index4D_pos] 	+= private_den_val*private_denPair_val;
					in_numerator_numPair[index4D_neg] 		+= private_num_val*private_numPair_val;
					in_denominator_denPair[index4D_neg] 	+= private_den_val*private_denPair_val;

					++idx6D;
				}

				++idx3D;
			}

			err << "\t - finished " << ++number_of_completed_events << " of " << total_N_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	err << "  * Finished!" << endl;

	return;
}



