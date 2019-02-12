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
	err << "  * Entering Compute_numerator_and_denominator_with_errors_q_mode_1D()" << endl;

	//const int q_space_size = n_qo_bins*n_qs_bins*n_ql_bins;
	//const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

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
		vector<complex<double> > sum1(in_numerator.size());
		vector<double> sum2(in_numerator.size());
		vector<double> sum3(in_denominator.size());
		vector<double> sum4(in_denominator.size());
		vector<double> sum5(in_denominator.size());

		vector<int> private_ABC(in_numerator_bin_count.size());
		vector<int> private_DBC(in_denominator_bin_count.size());

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Loop over K-bins
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

				bool num_bin_true = 	abs( Kx - pi.px ) <= 0.5*px_bin_width
									and abs( Ky - pi.py ) <= 0.5*py_bin_width
									and abs( Kz - pi.pz ) <= 0.5*pz_bin_width;

				// If particle contributes to this K-bin
				if ( num_bin_true )
				{
					double num_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;

					// Loop over q-bins;
					// qs-integral performed with delta-function
					for (int iQ = 0; iQ < n_Q_bins; iQ++)
					for (int iqo = 0; iqo < n_qo_bins; iqo++)
					for (int iql = 0; iql < n_ql_bins; iql++)
					{
						double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
						double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
						double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

						int index4D = indexer_qmode_1(iKT, iKphi, iKL, iQ);
		
						const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
						const double xi1 = qo*KT+ql*KL;
						const double xi3 = Q0*Q0 - qo*qo - ql*ql;
		
						// Check if a solution even exists;
						// if not, we're in a (q,K)-bin which
						// doesn't contribute to this value of Q0!
						double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
						if ( disc < 0.0 )
						{
							//err << "Warning: no qs value solves this combination of Q0 and (q,K)!" << endl
							//	<< "\t KT      KL      qo      ql      Q0" << endl
							//	<< "\t " << KT << "      " << KL << "      "
							//	<< qo << "      " << ql << "      " << Q0 << endl;
							continue;
						}

						// Otherwise, set the positive root first
						double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );
		
						// Sum over +/- roots in q_s direction
						for (int i_qs_root = 0; i_qs_root <= 1; i_qs_root++)
						{

							double qx = qo * cKphi - qs0 * sKphi;
							double qy = qs0 * cKphi + qo * sKphi;
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

							double q0 = Ea - Eb;
							double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

							complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;

							sum1[index4D] += complex_num_term;
							sum2[index4D] += 1.0 / (num_bin_factor*num_bin_factor);
							private_ABC[index4D]++;

							// loop back and do negative root
							qs0 *= -1.0;

						}

					}	// end of q loops

				}		// end of if block

			}			// end of particle loop

		} 				// end of K loops




		//=====================================
		//======== Doing denominator ==========
		//=====================================

		// Loop over K-bins
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

			// Loop over q-bins;
			// qs-integral performed with delta-function
			for (int iQ = 0; iQ < n_Q_bins; iQ++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				int index4D = indexer_qmode_1(iKT, iKphi, iKL, iQ);
	
				const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
				const double xi1 = qo*KT+ql*KL;
				const double xi3 = Q0*Q0 - qo*qo - ql*ql;

				// Check if a solution even exists;
				// if not, we're in a (q,K)-bin which
				// doesn't contribute to this value of Q0!
				double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
				if ( disc < 0.0 )
				{
					//err << "Warning: no qs value solves this combination of Q0 and (q,K)!" << endl
					//	<< "\t KT      KL      qo      ql      Q0" << endl
					//	<< "\t " << KT << "      " << KL << "      "
					//	<< qo << "      " << ql << "      " << Q0 << endl;
					continue;
				}

				// Otherwise, set the positive root first
				double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );

				// Sum over +/- roots in q_s direction
				for (int i_qs_root = 0; i_qs_root <= 1; i_qs_root++)
				{

					double qx = qo * cKphi - qs0 * sKphi;
					double qy = qs0 * cKphi + qo * sKphi;
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

					double den_bin_factor =
						px_bin_width*py_bin_width*pz_bin_width;

					double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

					// Sum over particles for denominator
					for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
					{
						ParticleRecord pi = event.particles[iParticle];

						// Check if in pa-bin
						bool in_den_pa_bin = 	abs(pi.px-pax) <= 0.5*px_bin_width
											and abs(pi.py-pay) <= 0.5*py_bin_width
											and abs(pi.pz-paz) <= 0.5*pz_bin_width;

						// Check if in pb-bin
						bool in_den_pb_bin = 	abs(pi.px-pbx) <= 0.5*px_bin_width
											and abs(pi.py-pby) <= 0.5*py_bin_width
											and abs(pi.pz-pbz) <= 0.5*pz_bin_width;

						if ( in_den_pa_bin )
						{
							sum3[index4D] += den_term;
							private_DBC[index4D]++;
						}

						if ( in_den_pb_bin )
						{
							sum4[index4D] += den_term;
							private_DBC[index4D]++;
						}

						if ( in_den_pa_bin and in_den_pb_bin )
						{
							// N.B. - use (bin volume)^2, here
							sum5[index4D] += den_term*den_term;
							private_DBC[index4D]++;
						}

					}		// end of denominator particle loops

					// loop back and do negative root
					qs0 *= -1.0;

				}			// end of loop over qs-roots

			}				// end of q loops

		}					// end of K loops

		// Need this to avoid race conditions
		#pragma omp critical
		{
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			for (int iQ = 0; iQ < n_Q_bins; iQ++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
				double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
				const double xi1 = qo*KT+ql*KL;
				const double xi3 = Q0*Q0 - qo*qo - ql*ql;

				// Check if a solution even exists;
				// if not, we're in a (q,K)-bin which
				// doesn't contribute to this value of Q0!
				double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
				if ( disc < 0.0 )
				{
					//err << "Warning: no qs value solves this combination of Q0 and (q,K)!" << endl
					//	<< "\t KT      KL      qo      ql      Q0" << endl
					//	<< "\t " << KT << "      " << KL << "      "
					//	<< qo << "      " << ql << "      " << Q0 << endl;
					continue;
				}
				/*else
				{
					err << "Found this solution:" << endl
						<< "\t KT      KL      qo      ql      Q0      +/-qs0" << endl
						<< "\t " << KT << "      " << KL << "      "
						<< qo << "      " << ql << "      " << Q0 << "      +/-"
						<< sqrt( disc / ( 4.0*xi0 + xi3 ) ) << endl;
				}*/

				// Otherwise, set the |root|
				double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );
				// cut this point out of Riemann sum over qo and ql
				if (abs(qs0) < 1.e-6)
					continue;

				// weight factor from delta-function identities
				// to get the normalization right
				const double weight_num = abs( (4.0*xi0+xi3)*(4.0*xi0+xi3) - 4.0*xi1*xi1 );
				const double weight_den = qs0*( (4.0*xi0+xi3)*(4.0*xi0+xi3) + 4.0*xi1*xi1 + weight_num );
				const double weight_factor = (qs0 < 1.e-6) ? 0.0 : weight_num / weight_den;
				//err << "Check: weight_factor = " << qs0*weight_factor << endl;

				const int index4D = indexer_qmode_1(iKT, iKphi, iKL, iQ);

				double abs_sum1 = abs( sum1[index4D] );
				double numerator_contribution_from_this_event
						= weight_factor
							* ( abs_sum1*abs_sum1 - sum2[index4D] );
				double denominator_contribution_from_this_event
						= weight_factor
							* ( sum3[index4D]*sum4[index4D] - sum5[index4D] );

				// first moments
				in_numerator[index4D]
					+= numerator_contribution_from_this_event;
				in_denominator[index4D]
					+= denominator_contribution_from_this_event;

				// second moments
				in_numerator2[index4D]
					+= numerator_contribution_from_this_event
						* numerator_contribution_from_this_event;
				in_denominator2[index4D]
					+= denominator_contribution_from_this_event
						* denominator_contribution_from_this_event;
				in_numerator_denominator[index4D]
					+= numerator_contribution_from_this_event
						* denominator_contribution_from_this_event;

				// track number of events where bin count was non-vanishing
				in_numerator_bin_count[index4D] += int(private_ABC[index4D] > 0);
				in_denominator_bin_count[index4D] += int(private_DBC[index4D] > 0);

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



