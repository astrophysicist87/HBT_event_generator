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



///*
void HBT_event_generator::Compute_numerator_and_denominator_methodMode2_q_mode_3D()
{
	// ignore q-dependence in denominator
	bool use_smoothness_approximation = false;
	const double SAfact = ( use_smoothness_approximation ) ? 0.0 : 1.0;

	bool perform_random_rotation = false;
	bool perform_random_shuffle = false;

	int number_of_completed_events = 0;
	cout << "  * Computing numerator and denominator of correlation function with errors; qmode = 3D using bin-averaging" << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	const double KYmin = -0.5, KYmax = 0.5;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_num(numerator.size(), 0.0);
		vector<double> private_num2(numerator2.size(), 0.0);
		vector<double> private_den(denominator.size(), 0.0);
		vector<double> private_den2(denominator2.size(), 0.0);
		vector<double> private_num_den(numerator.size(), 0.0);

		//===================================
		//======== Doing numerator ==========
		//===================================

		double num_pairs_this_event = static_cast<double>(event.particles.size() * (event.particles.size() - 1));

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = iParticle + 1; jParticle < event.particles.size(); ++jParticle)
		{
			//if (iParticle == jParticle)
			//	continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// check which cell we're in
			double cx = pix - pjx;
			double cy = piy - pjy;
			double cz = piz - pjz;

			bool num_bin_false = 	   cx*cx > px_bin_width*px_bin_width
									or cy*cy > py_bin_width*py_bin_width
									or cz*cz > pz_bin_width*pz_bin_width;

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
			{
				//err << "KT; shouldn't have made it here!" << endl;
				continue;
			}

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
			{
				//err << "Kphi; shouldn't have made it here!" << endl;
				continue;
			}

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
			{
				//err << "KL; shouldn't have made it here!" << endl;
				continue;
			}

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				bool in_center = ( iqo == (n_qo_bins-1)/2
						and iqs == (n_qs_bins-1)/2
						and iql == (n_ql_bins-1)/2 );
                                double overall_factor = 1.0;
                                if ( not in_center and
					not use_smoothness_approximation )
                                        overall_factor = 2.0;

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				double arg =  q0 * (ti - tj)
							- qx * (xi - xj)
							- qy * (yi - yj)
							- qz * (zi - zj);

				double num_term = cos(arg/hbarC);	// factor of 2 allows to collapse sum in half
				//						/ ( 1.0
				//							* px_bin_width
				//							* py_bin_width
				//							* pz_bin_width );

				private_num[index6D] += overall_factor * num_term;
				//						/ num_pairs_this_event;

			}
		}

		//=====================================
		//========= Doing denominator =========
		//=====================================

		// Sum over pairs of particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = iParticle + 1; jParticle < event.particles.size(); ++jParticle)
		{
			//if ( iParticle == jParticle )
			//	continue;

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

			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				//bool in_center = ( iqo == (n_qo_bins-1)/2 and iqs == (n_qs_bins-1)/2 and iql == (n_ql_bins-1)/2 );
				double overall_factor = 1.0;
				//if ( not use_smoothness_approximation )
				//{
				//	overall_factor = 0.5;
				//	if ( in_center )
				//		overall_factor = 1.0;
				//}

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				// check which cell we're in
				double cmx = pix - pjx - SAfact*qx;
				double cmy = piy - pjy - SAfact*qy;
				double cmz = piz - pjz - SAfact*qz;
				double cpx = pix - pjx + SAfact*qx;
				double cpy = piy - pjy + SAfact*qy;
				double cpz = piz - pjz + SAfact*qz;

				// modified binning condition
				bool this_pair_den_bin_false
						= 	   cmx*cmx > px_bin_width*px_bin_width
							or cmy*cmy > py_bin_width*py_bin_width
							or cmz*cmz > pz_bin_width*pz_bin_width;
				bool rev_pair_den_bin_false
						= 	   cpx*cpx > px_bin_width*px_bin_width
							or cpy*cpy > py_bin_width*py_bin_width
							or cpz*cpz > pz_bin_width*pz_bin_width;

				if ( this_pair_den_bin_false and rev_pair_den_bin_false )
					continue;

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				//private_den[index6D]++;
				private_den[index6D] += overall_factor;

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
				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double numerator_this_event = private_num[idx6D];
					double denominator_this_event = private_den[idx6D];

					// first moments
					numerator[idx6D]
								+= numerator_this_event;
					denominator[idx6D]
								+= denominator_this_event;

					// second moments
					numerator2[idx6D]
								+= numerator_this_event
									* numerator_this_event;
					denominator2[idx6D]
								+= denominator_this_event
									* denominator_this_event;
					numerator_denominator[idx6D]
								+= numerator_this_event
									* denominator_this_event;

					++idx6D;
				}

				++idx3D;
			}

			++number_of_completed_events;

			cout << "\t - finished "
					<< number_of_completed_events + total_N_events - allEvents.size()
					<< " of " << total_N_events << endl;
		}

	}

	cout << "  * Finished " << total_N_events << " events so far!" << endl;

	return;
}
//*/


/*
void HBT_event_generator::Compute_numerator_and_denominator_methodMode2_q_mode_3D()
{
	//int number_of_completed_events = 0;
	//err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	constexpr bool oneDim_slices = true;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	const int q_space_size = n_qo_bins*n_qs_bins*n_ql_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;
	const int sums_lengths = q_space_size * K_space_size * n_KT_pts_per_bin * n_Kphi_pts_per_bin * n_KL_pts_per_bin;

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<complex<double> > sum1(sums_lengths);
		vector<double> sum2(sums_lengths);
		vector<double> sum3(sums_lengths);
		vector<double> sum4(sums_lengths);
		vector<double> sum5(sums_lengths);

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KTmin = KT_pts[iKT];
			double KTmax = KT_pts[iKT+1];
			double Kphimin = Kphi_pts[iKT];
			double Kphimax = Kphi_pts[iKT+1];
			double KLmin = KL_pts[iKT];
			double KLmax = KL_pts[iKT+1];

			double KT_hw = 0.5*( KTmax - KTmin ), KT_cen = 0.5*( KTmax + KTmin );
			double Kphi_hw = 0.5*( Kphimax - Kphimin ), Kphi_cen = 0.5*( Kphimax + Kphimin );
			double KL_hw = 0.5*( KLmax - KLmin ), KL_cen = 0.5*( KLmax + KLmin );

			for (int ijKT 	= 0; ijKT 	< n_KT_pts_per_bin; 	ijKT++)
			for (int ijKphi = 0; ijKphi < n_Kphi_pts_per_bin; 	ijKphi++)
			for (int ijKL 	= 0; ijKL 	< n_KL_pts_per_bin; 	ijKL++)
			{

				double KT = KT_hw * xKT_pts[ijKT] + KT_cen;
				double Kphi = Kphi_hw * xKphi_pts[ijKphi] + Kphi_cen;
				double KL = KL_hw * xKL_pts[ijKL] + KL_cen;
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

					if ( num_bin_true )
					{
						double num_bin_factor =
								px_bin_width*py_bin_width*pz_bin_width;
						//num_bin_factor = 1.0;

						for (int iqo = 0; iqo < n_qo_bins; iqo++)
						for (int iqs = 0; iqs < n_qs_bins; iqs++)
						for (int iql = 0; iql < n_ql_bins; iql++)
						{
							// allows to do only slices
							if ( oneDim_slices )
							{
								bool iqo_not_center = ( iqo != (n_qo_bins-1)/2 );
								bool iqs_not_center = ( iqs != (n_qs_bins-1)/2 );
								bool iql_not_center = ( iql != (n_ql_bins-1)/2 );

								// if we're not on an axis slice, skip this q-bin
								if (
										( iqo_not_center and iqs_not_center )
										or ( iqo_not_center and iql_not_center )
										or ( iqs_not_center and iql_not_center )
									)
									continue;
							}

							int index9D = indexer(iKT, iKphi, iKL, ijKT, ijKphi, ijKL, iqo, iqs, iql);

							double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
							double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
							double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
							double qx = qo * cKphi - qs * sKphi;
							double qy = qs * cKphi + qo * sKphi;
							double qz = ql;

							// rapidity cuts
							if ( impose_pair_rapidity_cuts )
							{
								double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
								double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;
								double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
								double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);

								bool too_far_backward = ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min );
								bool too_far_forward = ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max );

								if ( too_far_backward or too_far_forward ) continue;
							}

							double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

							double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

							complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;
							//err << "CHECKCOMPLEX " << arg/hbarC << "   " << i*arg/hbarC << "   " << exp(i*arg/hbarC) << endl;

							sum1[index9D] += complex_num_term;
							sum2[index9D] += 1.0 / (num_bin_factor*num_bin_factor);

						}	// end of q loops

					}		// end of if block

				}			// end of particle loop

			}				// end of loop over bin integration points

		}					// end of K loops




		//=====================================
		//======== Doing denominator ==========
		//=====================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KTmin = KT_pts[iKT];
			double KTmax = KT_pts[iKT+1];
			double Kphimin = Kphi_pts[iKT];
			double Kphimax = Kphi_pts[iKT+1];
			double KLmin = KL_pts[iKT];
			double KLmax = KL_pts[iKT+1];

			double KT_hw = 0.5*( KTmax - KTmin ), KT_cen = 0.5*( KTmax + KTmin );
			double Kphi_hw = 0.5*( Kphimax - Kphimin ), Kphi_cen = 0.5*( Kphimax + Kphimin );
			double KL_hw = 0.5*( KLmax - KLmin ), KL_cen = 0.5*( KLmax + KLmin );


			for (int ijKT 	= 0; ijKT 	< n_KT_pts_per_bin; 	ijKT++)
			for (int ijKphi = 0; ijKphi < n_Kphi_pts_per_bin; 	ijKphi++)
			for (int ijKL 	= 0; ijKL 	< n_KL_pts_per_bin; 	ijKL++)
			{

				double KT = KT_hw * xKT_pts[ijKT] + KT_cen;
				double Kphi = Kphi_hw * xKphi_pts[ijKphi] + Kphi_cen;
				double KL = KL_hw * xKL_pts[ijKL] + KL_cen;
				double cKphi = cos(Kphi), sKphi = sin(Kphi);
				double Kx = KT * cKphi;
				double Ky = KT * sKphi;
				double Kz = KL;

				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{

					// allows to do only slices
					if ( oneDim_slices )
					{
						bool iqo_not_center = ( iqo != (n_qo_bins-1)/2 );
						bool iqs_not_center = ( iqs != (n_qs_bins-1)/2 );
						bool iql_not_center = ( iql != (n_ql_bins-1)/2 );

						// if we're not on an axis slice, skip this q-bin
						if (
								( iqo_not_center and iqs_not_center )
								or ( iqo_not_center and iql_not_center )
								or ( iqs_not_center and iql_not_center )
							)
							continue;
					}

					double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
					double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
					double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
					double qx = qo * cKphi - qs * sKphi;
					double qy = qs * cKphi + qo * sKphi;
					double qz = ql;

					double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
					double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;

					// rapidity cuts
					if ( impose_pair_rapidity_cuts )
					{
						double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
						double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);

						bool too_far_backward = ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min );
						bool too_far_forward = ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max );

						if ( too_far_backward or too_far_forward ) continue;
					}

					// Sum over particles for pa bin
					for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
					{
						ParticleRecord pi = event.particles[iParticle];

						bool den_bin_true = 	abs(pi.px-pax) <= 0.5*px_bin_width
											and abs(pi.py-pay) <= 0.5*py_bin_width
											and abs(pi.pz-paz) <= 0.5*pz_bin_width;

						if ( den_bin_true )
						{

							int index9D = indexer(iKT, iKphi, iKL, ijKT, ijKphi, ijKL, iqo, iqs, iql);

							double den_bin_factor =
								px_bin_width*py_bin_width*pz_bin_width;
							//den_bin_factor = 1.0;

							denominator_cell_was_filled[index9D] = true;

							double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

							sum3[index9D] += den_term;

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

							int index9D = indexer(iKT, iKphi, iKL, ijKT, ijKphi, ijKL, iqo, iqs, iql);

							double den_bin_factor =
								px_bin_width*py_bin_width*pz_bin_width;
							//den_bin_factor = 1.0;

							denominator_cell_was_filled[index9D] = true;

							double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

							sum4[index9D] += den_term;

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

							int index9D = indexer(iKT, iKphi, iKL, ijKT, ijKphi, ijKL, iqo, iqs, iql);

							// N.B. - use (bin volume)^2, here
							double den_bin_factor =
								px_bin_width*py_bin_width*pz_bin_width
								*px_bin_width*py_bin_width*pz_bin_width;
							//den_bin_factor = 1.0;

							denominator_cell_was_filled[index9D] = true;

							double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

							sum5[index9D] += den_term;

						}	// end of if block

					}		// end of sum5 particle loop

				}			// end of q loops

			}				// end of loop over bin integration points

		}					// end of K loops

		// Need this to avoid race conditions
		#pragma omp critical
		{
			const int iqoC = (n_qo_pts - 1) / 2;
			const int iqsC = (n_qs_pts - 1) / 2;
			const int iqlC = (n_ql_pts - 1) / 2;

			int idx = 0;
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				double KTmin = KT_pts[iKT];
				double KTmax = KT_pts[iKT+1];
				double Kphimin = Kphi_pts[iKT];
				double Kphimax = Kphi_pts[iKT+1];
				double KLmin = KL_pts[iKT];
				double KLmax = KL_pts[iKT+1];

				double KT_hw = 0.5*( KTmax - KTmin ), Kphi_hw = 0.5*( Kphimax - Kphimin ), KL_hw = 0.5*( KLmax - KLmin );
				double KT_cen = 0.5*( KTmax + KTmin );

				double numerator_this_event = 0.0;
				double denominator_this_event = 0.0;

				// integrate over this K-bin
				for (int ijKT = 0; ijKT < n_KT_pts_per_bin; ijKT++)
				for (int ijKphi = 0; ijKphi < n_Kphi_pts_per_bin; ijKphi++)
				for (int ijKL = 0; ijKL < n_KL_pts_per_bin; ijKL++)
				{
					int index9D = indexer(iKT, iKphi, iKL, ijKT, ijKphi, ijKL, iqo, iqs, iql);

					double KT = KT_cen + KT_hw * xKT_pts[ijKT];
					double bin_weight = KT * xKT_wts[ijKT] * xKphi_wts[ijKphi] * xKL_wts[ijKL];	// note jacobian

					numerator_this_event += bin_weight * ( norm(sum1[index9D]) - sum2[index9D] );	// norm returns |z|^2
					denominator_this_event += bin_weight * ( sum3[index9D]*sum4[index9D] - sum5[index9D] );
				}

				// first moments
				numerator[idx] += numerator_this_event;
				denominator[idx] += denominator_this_event;

				// second moments
				numerator2[idx]
					+= numerator_this_event
						* numerator_this_event;
				denominator2[idx]
					+= denominator_this_event
						* denominator_this_event;
				numerator_denominator[idx]
					+= numerator_this_event
						* denominator_this_event;

				++idx;
			}

			//err << "\t - finished " << ++number_of_completed_events << " of " << total_N_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	//err << "  * Finished!" << endl;

	return;
}
*/






