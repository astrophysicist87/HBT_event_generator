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
	cout << "  * Computing numerator and denominator of correlation function with errors; qmode = 3D using bin-averaging" << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	const double KYmin = -0.5, KYmax = 0.5;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	if (true)
	{
		err << "Compute_numerator_and_denominator_with_errors_q_mode_1D_wBA() not yet fully validated!"
			<< "  Please try again later." << endl;
		exit(8);
	}

	// Sum over all events
	#pragma omp parallel for schedule(static) shared( in_numerator, in_numerator2,\
														in_denominator, in_denominator2,\
														in_numerator_denominator )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_num(2*in_numerator.size(), 0.0);
		vector<double> private_num2(2*in_numerator2.size(), 0.0);
		vector<double> private_den(2*in_denominator.size(), 0.0);
		vector<double> private_den2(2*in_denominator2.size(), 0.0);
		vector<double> private_num_den(2*in_numerator.size(), 0.0);

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

			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			bool num_bin_false = 	   abs( pix - pjx ) > px_bin_width
									or abs( piy - pjy ) > py_bin_width
									or abs( piz - pjz ) > pz_bin_width;

			if ( num_bin_false )
				continue;

			double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
			double tj = pj.t, xj = pj.x, yj = pj.y, zj = pj.z;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);

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

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);

			// Loop over q-bins;
			// qs-integral performed with delta-function
			for (int iQ = 0; iQ < n_Q_bins; iQ++)
			for (int iqRP = 0; iqRP < n_qRP_pts; iqRP++)	//using points, not bins!
			for (int ithq = 0; ithq < n_thq_pts; ithq++)	//using points, not bins!
			{
				double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
				if (abs(Q0)<1.e-20)
					Q0 = 1.e-20;
				double loc_alpha = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;

				double costthetaq = cos(ttheta_q_pts[ithq]);
				double sintthetaq = sin(ttheta_q_pts[ithq]);
				double qRP_min = 0.0;
				double num_loc = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;
				double den_loc = 4.0*(particle_mass*particle_mass + sintthetaq*sintthetaq*(KT*KT + KL*KL)) + Q0*Q0;
				double qRP_max = abs(Q0)*sqrt(num_loc/den_loc);
				double qRP_cen = 0.5*(qRP_max + qRP_min), qRP_hw = 0.5*(qRP_max - qRP_min);

				double qRP = qRP_cen + qRP_hw * x_pts[iqRP];
				double thetaq = ttheta_q_pts[ithq] + thetaK;	//shifted w.r.t. K
				double costhetaq = cos(thetaq);
				double sinthetaq = sin(thetaq);

				double qo = qRP * sinthetaq;
				double ql = qRP * costhetaq;

				const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
				const double xi1 = qo*KT+ql*KL;
				const double xi3 = Q0*Q0 - qo*qo - ql*ql;

				// Check if a solution even exists;
				// if not, we're in a (q,K)-bin which
				// doesn't contribute to this value of Q0!
				double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
				if ( disc < 0.0 )
					continue;

				// Otherwise, set the positive root first
				double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );
				//if (abs(qs0) < 1.e-6)
				//	continue;

				// Record +/- roots in q_s direction
				for (int i_qs_root = 0; i_qs_root <= 1; i_qs_root++)
				{
					int index7D = indexer_qmode_1(iKT, iKphi, iKL, iQ, iqRP, ithq, i_qs_root);

					double qx = qo * cKphi - qs0 * sKphi;
					double qy = qs0 * cKphi + qo * sKphi;
					double qz = ql;

					// loop back and do negative root
					qs0 *= -1.0;

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

					double q0  =  Ea - Eb;
					double arg =  q0 * (ti - tj)
								- qx * (xi - xj)
								- qy * (yi - yj)
								- qz * (zi - zj);

					double num_term = cos(arg/hbarC);

					private_num[index7D] += num_term;

				}
			}
		}


		//=====================================
		//========= Doing denominator =========
		//=====================================

		// Sum over pairs of particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			if ( iParticle == jParticle )
				continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);

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

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);

			// Loop over q-bins;
			// qs-integral performed with delta-function
			for (int iQ = 0; iQ < n_Q_bins; iQ++)
			for (int iqRP = 0; iqRP < n_qRP_pts; iqRP++)	//using points, not bins!
			for (int ithq = 0; ithq < n_thq_pts; ithq++)	//using points, not bins!
			{

				double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
				if (abs(Q0)<1.e-20)
					Q0 = 1.e-20;

				double loc_alpha = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;

				double costthetaq = cos(ttheta_q_pts[ithq]);
				double sintthetaq = sin(ttheta_q_pts[ithq]);
				double qRP_min = 0.0;
				double num_loc = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;
				double den_loc = 4.0*(particle_mass*particle_mass + sintthetaq*sintthetaq*(KT*KT + KL*KL)) + Q0*Q0;
				double qRP_max = abs(Q0)*sqrt(num_loc/den_loc);
				double qRP_cen = 0.5*(qRP_max + qRP_min), qRP_hw = 0.5*(qRP_max - qRP_min);

				double qRP = qRP_cen + qRP_hw * x_pts[iqRP];
				double thetaq = ttheta_q_pts[ithq] + thetaK;	//shifted w.r.t. K
				double costhetaq = cos(thetaq);
				double sinthetaq = sin(thetaq);

				double qo = qRP * sinthetaq;
				double ql = qRP * costhetaq;

				const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
				const double xi1 = qo*KT+ql*KL;
				const double xi3 = Q0*Q0 - qo*qo - ql*ql;

				// Check if a solution even exists;
				// if not, we're in a (q,K)-bin which
				// doesn't contribute to this value of Q0!
				double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
				if ( disc < 0.0 )
					continue;

				// Otherwise, set the positive root first
				double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );
				//if (abs(qs0) < 1.e-6)
				//	continue;

				// Record +/- roots in q_s direction
				for (int i_qs_root = 0; i_qs_root <= 1; i_qs_root++)
				{
					int index7D = indexer_qmode_1(iKT, iKphi, iKL, iQ, iqRP, ithq, i_qs_root);

					double qx = qo * cKphi - qs0 * sKphi;
					double qy = qs0 * cKphi + qo * sKphi;
					double qz = ql;

					// loop back and do negative root
					qs0 *= -1.0;

					// modified binning condition
					bool den_bin_false = 	   abs( pix - pjx - qx ) > px_bin_width
											or abs( piy - pjy - qy ) > py_bin_width
											or abs( piz - pjz - qz ) > pz_bin_width;

					if ( den_bin_false )
						continue;

					private_den[index7D]++;

				}	// end of loop over qs-roots
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
				/*for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double numerator_this_event = private_num[idx6D];
					double denominator_this_event = private_den[idx6D];

					// first moments
					in_numerator[idx6D]
								+= numerator_this_event;
					in_denominator[idx6D]
								+= denominator_this_event;

					// second moments
					in_numerator2[idx6D]
								+= numerator_this_event
									* numerator_this_event;
					in_denominator2[idx6D]
								+= denominator_this_event
									* denominator_this_event;
					in_numerator_denominator[idx6D]
								+= numerator_this_event
									* denominator_this_event;

					++idx6D;
				}*/
				for (int iQ = 0; iQ < n_Q_bins; iQ++)
				{
					const double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
					const double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
					const double thetaK = atan2(KT, KL);	// Usage: atan2(double y, double x)
													// y-->out, x-->long
					const double Kmag = sqrt(KT*KT+KL*KL);
					double Q0 = 0.5*(Q_pts[iQ]+Q_pts[iQ+1]);
					if (abs(Q0)<1.e-20)
						Q0 = 1.e-20;
					const double loc_alpha = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;

					// integrals over numerator and denominator
					double numerator_this_event = 0.0;
					double denominator_this_event = 0.0;

					const int index4D = indexer_qmode_1(iKT, iKphi, iKL, iQ);

					for (int ithq = 0; ithq < n_thq_pts; ithq++)	//using points, not bins!
					{
						const double thetaq = ttheta_q_pts[ithq] + thetaK;	//shifted w.r.t. K
						const double costthetaq = cos(ttheta_q_pts[ithq]);
						const double sintthetaq = sin(ttheta_q_pts[ithq]);
						const double costhetaq = cos(thetaq);
						const double sinthetaq = sin(thetaq);
						const double qRP_min = 0.0;
						double num_loc = 4.0*(particle_mass*particle_mass + KT*KT + KL*KL) + Q0*Q0;
						double den_loc = 4.0*(particle_mass*particle_mass + sintthetaq*sintthetaq*(KT*KT + KL*KL)) + Q0*Q0;
						double qRP_max = abs(Q0)*sqrt(num_loc/den_loc);
						const double qRP_cen = 0.5*(qRP_max + qRP_min), qRP_hw = 0.5*(qRP_max - qRP_min);

						for (int iqRP = 0; iqRP < n_qRP_pts; iqRP++)	//using points, not bins!
						{
							const double qRP = qRP_cen + qRP_hw * x_pts[iqRP];
							const double qRPwt = qRP_hw * x_wts[iqRP];
							const double qo = qRP * sinthetaq;
							const double ql = qRP * costhetaq;

							const double xi0 = particle_mass*particle_mass + KT*KT + KL*KL + 0.25*(qo*qo+ql*ql);
							const double xi1 = qo*KT+ql*KL;
							const double xi3 = Q0*Q0 - qo*qo - ql*ql;

							// Check if a solution even exists;
							// if not, we're in a (q,K)-bin which
							// doesn't contribute to this value of Q0!
							const double disc = 4.0*xi1*xi1 + 4.0*xi0*xi3 + xi3*xi3;
							if ( disc < 0.0 )
							{
								err << "Shouldn't have reached this point!" << endl;
								continue;
							}

							// Otherwise, set the |root|
							const double qs0 = sqrt( disc / ( 4.0*xi0 + xi3 ) );
							//if (abs(qs0) < 1.e-6)
							//	continue;

							// weight factor from delta-function identities
							// to get the normalization right
							const double weight_num = abs( (4.0*xi0+xi3)*(4.0*xi0+xi3) - 4.0*xi1*xi1 );
							const double weight_den = 1.e-100+qs0*( (4.0*xi0+xi3)*(4.0*xi0+xi3) + 4.0*xi1*xi1 + weight_num );
							const double weight_factor = /*(qs0 < 1.e-6) ? 0.0 :*/ weight_num / weight_den;
							const double integration_weight = qRP * qRPwt * ttheta_q_wts[ithq];

							// Sum over +/- roots in q_s direction
							for (int i_qs_root = 0; i_qs_root <= 1; i_qs_root++)
							{
								const int index7D = indexer_qmode_1(iKT, iKphi, iKL, iQ, iqRP, ithq, i_qs_root);

								numerator_this_event 	+= integration_weight * weight_factor * private_num[index7D];
								denominator_this_event 	+= integration_weight * weight_factor * private_den[index7D];

							}
						}
					}

					// input vectors have length of 4D space
					// first moments
					in_numerator[index4D]
						+= numerator_this_event;
					in_denominator[index4D]
						+= denominator_this_event;

					// second moments
					in_numerator2[index4D]
						+= numerator_this_event
							* numerator_this_event;
					in_denominator2[index4D]
						+= denominator_this_event
							* denominator_this_event;
					in_numerator_denominator[index4D]
						+= numerator_this_event
							* denominator_this_event;

				}

				++idx3D;
			}

			++number_of_completed_events;

			//cout << "\t - finished "
			//		<< number_of_completed_events + total_N_events - allEvents.size()
			//		<< " of " << total_N_events << endl;
		}

	}

	cout << "  * Finished " << total_N_events << " events so far!" << endl;

	return;
}



