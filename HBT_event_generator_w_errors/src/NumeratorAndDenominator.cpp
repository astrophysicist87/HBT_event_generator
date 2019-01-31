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


void HBT_event_generator::Compute_numerator_and_denominator_with_errors_q_mode_3D(
							vector<double> & in_numerator,
							vector<double> & in_numerator2,
							vector<double> & in_denominator,
							vector<double> & in_denominator2,
							vector<double> & in_numerator_denominator,
							vector<int> & in_numerator_bin_count,
							vector<int> & in_denominator_bin_count
							)
{
	//int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// check multiplicities
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
		err << "allEvents[" << iEvent << "].particles.size() = " << allEvents[iEvent].particles.size() << endl;

	// Sum over all events
	#pragma omp parallel for schedule(static) shared( in_numerator, in_numerator2,\
														in_denominator, in_denominator2,\
														in_numerator_denominator,\
														in_numerator_bin_count,\
														in_denominator_bin_count )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

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
						if ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min ) or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
							continue;

						double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

						double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

						complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;

						sum1[index6D] += complex_num_term;
						sum2[index6D] += 1.0 / (num_bin_factor*num_bin_factor);
						private_ABC[index6D]++;

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
				//double paT = sqrt(pax*pax+pay*pay), pbT = sqrt(pbx*pbx+pby*pby);
				//double pa_phi = atan2(pay,pax), pb_phi = atan2(pby,pbx);

				// rapidity cuts
				if ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min ) or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
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

						denominator_cell_was_filled[index6D] = true;

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

						denominator_cell_was_filled[index6D] = true;

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

						denominator_cell_was_filled[index6D] = true;

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
				double abs_sum1 = abs(sum1[idx]);
				double numerator_contribution_from_this_event
						= abs_sum1*abs_sum1 - sum2[idx];
				double denominator_contribution_from_this_event
						= sum3[idx]*sum4[idx] - sum5[idx];

				// first moments
				in_numerator[idx]
					+= numerator_contribution_from_this_event;
				in_denominator[idx]
					+= denominator_contribution_from_this_event;

				// second moments
				in_numerator2[idx]
					+= numerator_contribution_from_this_event
						* numerator_contribution_from_this_event;
				in_denominator2[idx]
					+= denominator_contribution_from_this_event
						* denominator_contribution_from_this_event;
				in_numerator_denominator[idx]
					+= numerator_contribution_from_this_event
						* denominator_contribution_from_this_event;

				// track number of events where bin count was non-vanishing
				in_numerator_bin_count[idx] += int(private_ABC[idx] > 0);
				in_denominator_bin_count[idx] += int(private_DBC[idx] > 0);
				//in_numerator_bin_count[idx] += private_ABC[idx];
				//in_denominator_bin_count[idx] += private_DBC[idx];

				++idx;
			}

			err << "\t - finished " << ++number_of_completed_events << " of " << total_N_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	err << "  * Finished!" << endl;

	return;
}


