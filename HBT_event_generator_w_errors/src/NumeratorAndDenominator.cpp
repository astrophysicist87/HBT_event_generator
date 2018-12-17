#include <iostream>
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


void HBT_event_generator::Compute_numerator_and_denominator_with_errors(
							vector<double> & in_numerator,
							vector<double> & in_numerator2,
							vector<double> & in_denominator,
							vector<double> & in_denominator2,
							vector<double> & in_numerator_denominator
							)
{
	int number_of_completed_events = 0;
	err << "  * Computing numerator and denominator of correlation function with errors" << endl;

	// Sum over all events
	#pragma omp parallel for schedule(static) shared( in_numerator, in_numerator2,\
														in_denominator, in_denominator2,\
														in_numerator_denominator )
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_A(in_numerator.size());
		vector<double> private_B(in_denominator.size());

		vector<double> private_A2(in_numerator2.size());
		vector<double> private_B2(in_denominator2.size());
		vector<double> private_AB(in_numerator_denominator.size());


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

			// Get particle mass
			double m = particle_mass;

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

				//===================================
				//======== Doing numerator ==========
				//===================================

				bool num_bin_true = abs( Kx - pi.px ) <= 0.5*px_bin_width
									and abs( Ky - pi.py ) <= 0.5*py_bin_width
									and abs( Kz - pi.pz ) <= 0.5*pz_bin_width
									and abs( Kx - pj.px ) <= 0.5*px_bin_width
									and abs( Ky - pj.py ) <= 0.5*py_bin_width
									and abs( Kz - pj.pz ) <= 0.5*pz_bin_width;

				if ( num_bin_true )
				{
					double num_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width
							*px_bin_width*py_bin_width*pz_bin_width;

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

						double q0 = get_q0(m, qo, qs, ql, KT, KL);

						double arg =  q0 * (ti - tj)
									- qx * (xi - xj)
									- qy * (yi - yj)
									- qz * (zi - zj);

						//complex<double> num_term = exp(i*arg/hbarC) * num_bin_factor;
						double num_term = cos(arg/hbarC) * num_bin_factor;

						private_A[index6D] += num_term;

					}
				}

				//=====================================
				//======== Doing denominator ==========
				//=====================================
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

					double p1x = Kx + 0.5 * qx, p1y = Ky + 0.5 * qy, p1z = Kz + 0.5 * qz;
					double p2x = Kx - 0.5 * qx, p2y = Ky - 0.5 * qy, p2z = Kz - 0.5 * qz;

					bool den_bin_true = abs(pi.px-p1x) <= 0.5*px_bin_width
										and abs(pi.py-p1y) <= 0.5*py_bin_width
										and abs(pi.pz-p1z) <= 0.5*pz_bin_width
										and abs(pj.px-p2x) <= 0.5*px_bin_width
										and abs(pj.py-p2y) <= 0.5*py_bin_width
										and abs(pj.pz-p2z) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width
							*px_bin_width*py_bin_width*pz_bin_width;

						denominator_cell_was_filled[index6D] = true;

						double den_term = den_bin_factor;	//no phase factor in denominator

						private_B[index6D] += den_term;

					}

				}

			}

		}

		// Need this to avoid race conditions
		#pragma omp critical
		{
			int idx = 0;
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				double private_A_val = private_A[idx];
				double private_B_val = private_B[idx];

				in_numerator[idx] += private_A_val;
				in_denominator[idx] += private_B_val;

				in_numerator2[idx] += private_A_val*private_A_val;
				in_denominator2[idx] += private_B_val*private_B_val;
				in_numerator_denominator[idx] += private_A_val*private_B_val;

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


