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


void HBT_event_generator::Compute_numerator_and_denominator_with_errors_binPairsMode_qmode3D(
							vector<double> & in_numerator, vector<double> & in_numerator2,
							vector<double> & in_numPair, vector<double> & in_numPair2,
							vector<double> & in_denominator, vector<double> & in_denominator2,
							vector<double> & in_denPair, vector<double> & in_denPair2,
							vector<double> & in_numerator_numPair, vector<double> & in_denominator_denPair
							)
{
	bool perform_random_rotation = true;
	bool perform_random_shuffle = true;

	int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	const double KYmin = -0.5, KYmax = 0.5;
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

		vector<double> private_num(in_numerator.size(), 0.0);
		vector<double> private_num2(in_numerator2.size(), 0.0);
		vector<double> private_numPair(in_numPair.size(), 0.0);
		vector<double> private_numPair2(in_numPair2.size(), 0.0);
		vector<double> private_den(in_denominator.size(), 0.0);
		vector<double> private_den2(in_denominator2.size(), 0.0);
		vector<double> private_denPair(in_denPair.size(), 0.0);
		vector<double> private_denPair2(in_denPair2.size(), 0.0);

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
			double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

			// rapidity cuts
			if ( ( Kz / K0 < Kz_over_K0_min ) or ( Kz / K0 > Kz_over_K0_max ) )
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

			double arg =  q0 * (ti - tj)
						- qx * (xi - xj)
						- qy * (yi - yj)
						- qz * (zi - zj);

			double num_term = cos(arg/hbarC);

			private_num[index6D] += num_term;
			private_numPair[index3D]++;

		}


		//=====================================
		//========= Doing denominator =========
		//=====================================

		// Randomly sample events to mix with
		//const unsigned int n_mixing_events = allEvents.size()-1;
		const unsigned int n_mixing_events = 100;

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
				if ( ( Kz / K0 < Kz_over_K0_min ) or ( Kz / K0 > Kz_over_K0_max ) )
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

				if ( index6D >= n_KT_bins*n_Kphi_bins*n_KL_bins*n_qo_bins*n_qs_bins*n_ql_bins )
				{
					cout << "CHECK INDICES: "
						<< index6D << " of max. " << n_KT_bins*n_Kphi_bins*n_KL_bins*n_qo_bins*n_qs_bins*n_ql_bins-1 << "   "
						<< KT_idx << "   " << n_KT_bins << "   " << KT_max << "   " << KT << "   " << KT_min << "   "
						<< Kphi_idx << "   " << n_Kphi_bins << "   " << Kphi_max << "   " << Kphi << "   " << Kphi_min << "   "
						<< KL_idx << "   " << n_KL_bins << "   " << KL_max << "   " << Kz << "   " << KL_min << "   "
						<< qo_idx << "   " << n_qo_bins << "   " << qo_max << "   " << qo << "   " << qo_min << "   "
						<< qs_idx << "   " << n_qs_bins << "   " << qs_max << "   " << qs << "   " << qs_min << "   "
						<< ql_idx << "   " << n_ql_bins << "   " << ql_max << "   " << ql << "   " << ql_min << "   " << delta_ql << endl;
					exit(8);
				}

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

				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double private_num_val 			= private_num[idx6D];
					double private_den_val 			= private_den[idx6D];

					in_numerator[idx6D] 			+= private_num_val;
					in_denominator[idx6D] 			+= private_den_val;

					in_numerator2[idx6D] 			+= private_num_val*private_num_val;
					in_numerator_numPair[idx6D] 	+= private_num_val*private_numPair_val;
					in_denominator2[idx6D] 			+= private_den_val*private_den_val;
					in_denominator_denPair[idx6D] 	+= private_den_val*private_denPair_val;

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



//========================================================================
//========================================================================
//========================================================================
//========================================================================
//========================================================================
//========================================================================
//========================================================================

// this one uses Chun's method to get correlator vs. Q == sqrt(-(p1-p2)^2)
void HBT_event_generator::Compute_numerator_and_denominator_with_errors_binPairsMode_qmode1D(
							vector<double> & in_numerator, vector<double> & in_numerator2,
							vector<double> & in_numPair, vector<double> & in_numPair2,
							vector<double> & in_denominator, vector<double> & in_denominator2,
							vector<double> & in_denPair, vector<double> & in_denPair2,
							vector<double> & in_numerator_numPair, vector<double> & in_denominator_denPair
							)
{
	bool perform_random_rotation = true;
	bool perform_random_shuffle = false;

	int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	const double KYmin = -0.5, KYmax = 0.5;
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

		vector<double> private_num(in_numerator.size(), 0.0);
		vector<double> private_num2(in_numerator2.size(), 0.0);
		vector<double> private_numPair(in_numPair.size(), 0.0);
		vector<double> private_numPair2(in_numPair2.size(), 0.0);
		vector<double> private_den(in_denominator.size(), 0.0);
		vector<double> private_den2(in_denominator2.size(), 0.0);
		vector<double> private_denPair(in_denPair.size(), 0.0);
		vector<double> private_denPair2(in_denPair2.size(), 0.0);

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			// *distinct*
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
			double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

			// rapidity cuts
			if ( ( Kz / K0 < Kz_over_K0_min ) or ( Kz / K0 > Kz_over_K0_max ) )
				continue;

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);

			// Get indices
			int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
			int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
				continue;

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				continue;

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
				continue;

			// If we're in a valid K-bin, bin in Q too
			double Qbar 	= sqrt(abs(qx*qx + qy*qy + qz*qz - q0*q0));
			int Q_idx_pos 	= floor((Qbar - Q_min) / delta_Q);
			// two roots of Q^2==...
			Qbar 			*= -1.0;
			int Q_idx_neg 	= floor((Qbar - Q_min) / delta_Q);

			// Check if Q-bin is out of range
			if (   Q_idx_pos < 0
				or Q_idx_pos >= n_Q_bins
				or Q_idx_neg < 0
				or Q_idx_neg >= n_Q_bins )
				continue;

			// If not, set appropriate indices in vectors
			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
			int index4D_pos = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, Q_idx_pos);
			int index4D_neg = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, Q_idx_neg);

			// BE weight
			double arg =  q0 * (ti - tj)
						- qx * (xi - xj)
						- qy * (yi - yj)
						- qz * (zi - zj);

			double num_term = cos(arg/hbarC);

			private_num[index4D_pos] += num_term;
			private_num[index4D_neg] += num_term;
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

		// If we are randomly sampling mixed events...
		if ( perform_random_shuffle )
			random_shuffle(indices.begin(), indices.end());
		vector<int> mixedEvents;
		for (int mix_idx = 0; mix_idx <= n_mixing_events; ++mix_idx)
			if ( indices[mix_idx] != iEvent
					and mixedEvents.size() < n_mixing_events )
			{
				mixedEvents.push_back( indices[mix_idx] );
			}

		// If events are rotated by random angle
		// in transverse plane, set random angles
		// in advance
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
				if ( ( Kz / K0 < Kz_over_K0_min ) or ( Kz / K0 > Kz_over_K0_max ) )
					continue;

				double KT = sqrt(Kx*Kx+Ky*Ky);
				double Kphi = atan2(Ky, Kx);
				double cKphi = cos(Kphi), sKphi = sin(Kphi);

				// If pair survived cuts, get indices
				int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
				int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
				int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

				// Momentum-space cuts
				if ( KT_idx < 0 or KT_idx >= n_KT_bins )
					continue;

				if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
					continue;

				if ( KL_idx < 0 or KL_idx >= n_KL_bins )
					continue;

				// If we're in a valid K-bin, bin in Q too
				double Qbar 	= sqrt(abs(qx*qx + qy*qy + qz*qz - q0*q0));
				int Q_idx_pos 	= floor((Qbar - Q_min) / delta_Q);
				// two roots of Q^2==...
				Qbar 			*= -1.0;
				int Q_idx_neg 	= floor((Qbar - Q_min) / delta_Q);

				// Check if Q-bin is out of range
				if (   Q_idx_pos < 0
					or Q_idx_pos >= n_Q_bins
					or Q_idx_neg < 0
					or Q_idx_neg >= n_Q_bins )
					continue;

				// If not, set appropriate indices in vectors
				int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
				int index4D_pos = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, Q_idx_pos);
				int index4D_neg = indexer_qmode_1(KT_idx, Kphi_idx, KL_idx, Q_idx_neg);

				private_den[index4D_pos]++;
				private_den[index4D_neg]++;
				private_denPair[index3D]++;

			}
		}


		// Need this to avoid race conditions
		#pragma omp critical
		{
			int idx3D = 0, idx4D = 0;
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			{
				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);

				double private_numPair_val 	= private_numPair[idx3D];
				double private_denPair_val 	= private_denPair[idx3D];

				in_numPair[idx3D] 			+= private_numPair_val;
				in_numPair2[idx3D] 			+= private_numPair_val*private_numPair_val;
				in_denPair[idx3D] 			+= private_denPair_val;
				in_denPair2[idx3D] 			+= private_denPair_val*private_denPair_val;

				for (int iQ = 0; iQ < n_Q_bins; iQ++)
				{

					double private_num_val 			= private_num[idx4D];
					double private_den_val 			= private_den[idx4D];

					in_numerator[idx4D] 			+= private_num_val;
					in_denominator[idx4D] 			+= private_den_val;

					in_numerator2[idx4D] 			+= private_num_val*private_num_val;
					in_numerator_numPair[idx4D] 	+= private_num_val*private_numPair_val;
					in_denominator2[idx4D] 			+= private_den_val*private_den_val;
					in_denominator_denPair[idx4D] 	+= private_den_val*private_denPair_val;

					++idx4D;
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


