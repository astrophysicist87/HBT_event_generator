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



/*
void HBT_event_generator::Compute_numerator_and_denominator_methodMode2_q_mode_3D()
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

				double num_term = cos(arg/hbarC);
				//						/ ( 1.0
				//							* px_bin_width
				//							* py_bin_width
				//							* pz_bin_width );

				private_num[index6D] += num_term;
				//						/ num_pairs_this_event;

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

				// modified binning condition
				bool den_bin_false = 	   abs( pix - pjx - qx ) > px_bin_width
										or abs( piy - pjy - qy ) > py_bin_width
										or abs( piz - pjz - qz ) > pz_bin_width;

				if ( den_bin_false )
					continue;

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				private_den[index6D]++;

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
*/







void HBT_event_generator::Compute_numerator_and_denominator_methodMode2_q_mode_3D()
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

	const int number_of_threads = omp_get_num_threads();
	vector<vector<double> > num_per_thread( number_of_threads, vector<double>( numerator.size(), 0.0 ) );
	vector<vector<double> > num2_per_thread( number_of_threads, vector<double>( numerator2.size(), 0.0 ) );
	vector<vector<double> > den_per_thread( number_of_threads, vector<double>( denominator.size(), 0.0 ) );
	vector<vector<double> > den2_per_thread( number_of_threads, vector<double>( denominator2.size(), 0.0 ) );
	vector<vector<double> > num_den_per_thread( number_of_threads, vector<double>( numerator.size(), 0.0 ) );

	// Define a reduction for std::vector<double> objects
	//#pragma omp declare reduction( vec_double_plus : std::vector<double> : \
	//			std::transform( omp_out.begin(), omp_out.end(), \
	//			omp_in.begin(), omp_out.begin(), std::plus<double>() ) ) \
	//			initializer(omp_priv = omp_orig)

	// Sum over all events
	#pragma omp parallel for schedule(static)
	//			reduction( vec_double_plus : numerator, denominator, \
	//			numerator2, denominator2, numerator_denominator )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		const int current_thread = omp_get_thread_num();

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

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				double arg =  q0 * (ti - tj)
							- qx * (xi - xj)
							- qy * (yi - yj)
							- qz * (zi - zj);

				double num_term = cos(arg/hbarC);

				private_num[index6D] += num_term;

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

				// modified binning condition
				bool den_bin_false = 	   abs( pix - pjx - qx ) > px_bin_width
										or abs( piy - pjy - qy ) > py_bin_width
										or abs( piz - pjz - qz ) > pz_bin_width;

				if ( den_bin_false )
					continue;

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				private_den[index6D]++;

			}
		}

		///*
		// Need this to avoid race conditions
		//#pragma omp critical
		//{
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
					//numerator[idx6D]
					num_per_thread[current_thread][idx6D]
								+= numerator_this_event;
					//denominator[idx6D]
					den_per_thread[current_thread][idx6D]
								+= denominator_this_event;

					// second moments
					//numerator2[idx6D]
					num2_per_thread[current_thread][idx6D]
								+= numerator_this_event
									* numerator_this_event;
					//denominator2[idx6D]
					den2_per_thread[current_thread][idx6D]
								+= denominator_this_event
									* denominator_this_event;
					//numerator_denominator[idx6D]
					num_den_per_thread[current_thread][idx6D]
								+= numerator_this_event
									* denominator_this_event;

					++idx6D;
				}

				++idx3D;
			}


		if (verbose)
		{
			#pragma omp critical
			{
				//++number_of_completed_events;

				cout << "\t - finished "
						<< ++number_of_completed_events
							+ total_N_events
							- allEvents.size()
						<< " of " << total_N_events << endl;
			}
		}
		//}*/

	}

	// Do my own make-shift reduction here
	const int idx6D = 0;
	for (int iThread = 0; iThread < number_of_threads; ++iThread)
	{
		numerator[idx6D] += num_per_thread[iThread][idx6D];
		denominator[idx6D] += den_per_thread[iThread][idx6D];
		numerator2[idx6D] += num2_per_thread[iThread][idx6D];
		denominator2[idx6D] += den2_per_thread[iThread][idx6D];
		numerator_denominator[idx6D] += num_den_per_thread[iThread][idx6D];
		idx6D++;
	}

	cout << "  * Finished " << total_N_events << " events so far!" << endl;

	return;
}



