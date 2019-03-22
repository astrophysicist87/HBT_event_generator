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



void HBT_event_generator::Compute_numerator_and_denominator_methodMode0_q_mode_1D()
{
	const int iqoC = (n_qo_pts - 1) / 2;
	const int iqsC = (n_qs_pts - 1) / 2;
	const int iqlC = (n_ql_pts - 1) / 2;

	err << "  * Entering Compute_numerator_and_denominator_with_errors_q_mode_1D()" << endl;

	const int q_space_size = n_Q_bins*n_qRP_pts*n_thq_pts;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		#pragma omp critical
		{
			cout << "\t -- working on event = " << event.eventID << endl;
		}

		//let these be fully six-dimensional (factor of 2 for +/- roots in qs)
		vector<complex<double> > sum1(2*q_space_size*K_space_size);
		vector<double> sum2(2*q_space_size*K_space_size);
		vector<double> sum3(2*q_space_size*K_space_size);
		vector<double> sum4(2*q_space_size*K_space_size);
		vector<double> sum5(2*q_space_size*K_space_size);

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
			double thetaK = atan2(KT, KL);	// Usage: atan2(double y, double x)
											// y-->out, x-->long
			double Kmag = sqrt(KT*KT+KL*KL);

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
					//double num_bin_factor =
					//		px_bin_width*py_bin_width*pz_bin_width;
					double num_bin_factor = 1.0;

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

						//int index6D = indexer_qmode_1(iKT, iKphi, iKL, iQ, iqRP, ithq);
		
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

							double q0 = Ea - Eb;
							double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

							complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;

							sum1[index7D] += complex_num_term;
							sum2[index7D] += 1.0 / (num_bin_factor*num_bin_factor);

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
			double thetaK = atan2(KT, KL);	// Usage: atan2(double y, double x)
											// y-->out, x-->long
			double Kmag = sqrt(KT*KT+KL*KL);

			// Loop over q-bins;
			// qs-integral performed with delta-function
			for (int iQ = 0; iQ < n_Q_bins; iQ++)
			for (int iqRP = 0; iqRP < n_qRP_pts; iqRP++)	//using points, not bins!
			for (int ithq = 0; ithq < n_thq_pts; ithq++)	//using points, not bins!
			{
				//int index6D = indexer_qmode_1(iKT, iKphi, iKL, iQ, iqRP, ithq);

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

					//double den_bin_factor =
					//	px_bin_width*py_bin_width*pz_bin_width;
					double den_bin_factor = 1.0;

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
							sum3[index7D] += den_term;

						if ( in_den_pb_bin )
							sum4[index7D] += den_term;

						if ( in_den_pa_bin and in_den_pb_bin )
							sum5[index7D] += den_term*den_term;	// N.B. - use (bin volume)^2, here

					}		// end of denominator particle loops

				}			// end of loop over qs-roots

			}				// end of q loops

		}					// end of K loops


		//===================================
		//======== Storing results ==========
		//===================================

		// Need this to avoid race conditions
		#pragma omp critical
		{
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
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
				double numerator_contribution_from_this_event = 0.0;
				double denominator_contribution_from_this_event = 0.0;

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

							// temporary sum* vectors have length of 7D space
							const double abs_sum1 = abs( sum1[index7D] );
							numerator_contribution_from_this_event
									+= integration_weight * weight_factor
										* ( abs_sum1*abs_sum1 - sum2[index7D] );
							denominator_contribution_from_this_event
									+= integration_weight * weight_factor
										* ( sum3[index7D]*sum4[index7D] - sum5[index7D] );

						}

/*err << "Check: weight_factor(etc.) = "
<< qo << "   " << qs0 << "   " << ql << "   " << Q0 << "   "
<< weight_factor << "   "
<< abs_sum1*abs_sum1 - sum2[index6D] << "   "
<< sum3[index6D]*sum4[index6D] - sum5[index6D]
<< endl;*/
					}
				}

				// input vectors have length of 4D space
				// first moments
				numerator[index4D]
					+= numerator_contribution_from_this_event;
				denominator[index4D]
					+= denominator_contribution_from_this_event;

				// second moments
				numerator2[index4D]
					+= numerator_contribution_from_this_event
						* numerator_contribution_from_this_event;
				denominator2[index4D]
					+= denominator_contribution_from_this_event
						* denominator_contribution_from_this_event;
				numerator_denominator[index4D]
					+= numerator_contribution_from_this_event
						* denominator_contribution_from_this_event;

			}

		}

	}

	err << "  * Exiting Compute_numerator_and_denominator_with_errors_q_mode_1D()" << endl;

	return;
}

//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================

void HBT_event_generator::Compute_numerator_and_denominator_momentum_space_only_q_mode_1D()
{
	bool perform_random_rotation = false;
	bool perform_random_shuffle = false;

	int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	//const int q_space_size = n_qo_bins*n_qs_bins*n_ql_bins;
	const int q_space_size = n_Q_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// Sum over all events
	#pragma omp parallel for schedule(static)
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

			// set Q2 value for this q-K cell
			double Q2bar = qo*qo + qs*qs + ql*ql - q0*q0;	// Q2bar>=0
			if (Q2bar < -1.e-6)
			{
				err << "Compute_numerator_and_denominator_with_errors_q_mode_1D(warning): "
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

			//double arg =  q0 * (ti - tj)
			//			- qx * (xi - xj)
			//			- qy * (yi - yj)
			//			- qz * (zi - zj);

			//double num_term = cos(arg/hbarC);
			double num_term = 1.0;

			//private_num[index6D] += num_term;
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
					err << "Compute_numerator_and_denominator_with_errors_q_mode_1D(warning): "
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

				numPair[idx3D] 			+= private_numPair_val;
				numPair2[idx3D] 			+= private_numPair_val*private_numPair_val;
				denPair[idx3D] 			+= private_denPair_val;
				denPair2[idx3D] 			+= private_denPair_val*private_denPair_val;

				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);

				for (int iQ = 0; iQ < n_Q_bins; iQ++)
				{
					double private_num_val 			= private_num[idx4D];
					double private_den_val 			= private_den[idx4D];

					// first moments
					numerator[idx4D] 			+= private_num_val;
					denominator[idx4D] 			+= private_den_val;

					// second moments
					numerator2[idx4D] 			+= private_num_val * private_num_val;
					denominator2[idx4D] 			+= private_den_val * private_den_val;

					// for error normalizing by numbers of pairs
					numerator_numPair[idx4D] 	+= private_num_val * private_numPair_val;
					denominator_denPair[idx4D] 	+= private_den_val * private_denPair_val;

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



