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


void HBT_event_generator::Compute_numerator_and_denominator_with_errors_q_mode_1D_wBA(
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
	out << "  * Computing numerator and denominator of correlation function with errors; qmode = 1D using bin-averaging" << endl;

	const int q_space_size = n_Q_bins;
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

			double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
			double tj = pj.t, xj = pj.x, yj = pj.y, zj = pj.z;
			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);

			// rapidity cuts
			if ( impose_pair_rapidity_cuts
					and ( ( 2.0*Kz/(Ei+Ej) < Kz_over_K0_min )
					or ( 2.0*Kz/(Ei+Ej) > Kz_over_K0_max ) )
				)
				continue;

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double KL = Kz;

			// Get indices
			int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
			int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	= floor((KL - KL_min)/KL_bin_width);

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
				continue;

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				continue;

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
				continue;

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

				double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

				// set Q2 value for this q-K cell
				double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
				if (Q2bar < -1.e-6)
				{
					err << "Compute_numerator_and_denominator_with_errors_q_mode_1D_wBA(warning): "
						<< "Q2bar = " << Q2bar << " < 0.0!" << endl;
					continue;
				}
				double Qbar = sqrt(abs(Q2bar));

				int iQ_pos = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_pos < 0 or iQ_pos >= n_Q_bins )
					continue;
				int index4D_pos = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, iQ_pos);

				// assume symmetry under Q --> -Q
				Qbar *= -1.0;
				int iQ_neg = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_neg < 0 or iQ_neg >= n_Q_bins )
					continue;
				int index4D_neg = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, iQ_neg);

				int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				double arg =  q0 * (ti - tj)
							- qx * (xi - xj)
							- qy * (yi - yj)
							- qz * (zi - zj);

				double num_term = cos(arg/hbarC);

				private_num[index4D_pos] += num_term;
				private_num[index4D_neg] += num_term;
				private_numPair[index3D]++;

			}
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

				// set Q2 value for this q-K cell
				double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
				if (Q2bar < -1.e-6)
				{
					err << "Compute_numerator_and_denominator_with_errors_q_mode_1D_wBA(warning): "
						<< "Q2bar = " << Q2bar << " < 0.0!" << endl;
					continue;
				}
				double Qbar = sqrt(abs(Q2bar));

				int iQ_pos = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_pos < 0 or iQ_pos >= n_Q_bins )
					continue;
				int index4D_pos = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, iQ_pos);

				// assume symmetry under Q --> -Q
				Qbar *= -1.0;
				int iQ_neg = floor((Qbar - Q_min)/delta_Q);
				if ( iQ_neg < 0 or iQ_neg >= n_Q_bins )
					continue;
				int index4D_neg = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, iQ_neg);

				int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, qo_idx, qs_idx, ql_idx);

				//private_den[index6D]++;
				private_den[index4D_pos]++;
				private_den[index4D_neg]++;
				private_denPair[index3D]++;

			}
		}

		//===================================
		//======== Storing results ==========
		//===================================

		// Need this to avoid race conditions
		#pragma omp critical
		{
			int idx3D = 0, idx6D = 0, idx4D = 0;
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

				for (int iQ = 0; iQ < n_Q_bins; iQ++)
				{
					double private_num_val 			= private_num[idx4D];
					double private_den_val 			= private_den[idx4D];

					// first moments
					in_numerator[idx4D] 			+= private_num_val;
					in_denominator[idx4D] 			+= private_den_val;

					// second moments
					in_numerator2[idx4D] 			+= private_num_val * private_num_val;
					in_denominator2[idx4D] 			+= private_den_val * private_den_val;

					// for error normalizing by numbers of pairs
					in_numerator_numPair[idx4D] 	+= private_num_val * private_numPair_val;
					in_denominator_denPair[idx4D] 	+= private_den_val * private_denPair_val;

					++idx4D;
				}

				++idx3D;
			}

			out << "\t - finished " << ++number_of_completed_events << " of " << total_N_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	out << "  * Finished!" << endl;

	return;
}



