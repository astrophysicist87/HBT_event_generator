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


void HBT_event_generator::Compute_numerator( vector<complex<double> > & in_numerator )
{
	// Sum over all events
	#pragma omp parallel for schedule(static) shared(in_numerator)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<complex<double> > private_numerator(in_numerator.size());

		//err << "Processing numerator of event#" << iEvent << ":" << endl;

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


			// Get particle mass (can also read this in)
			//double m = sqrt(pi.E*pi.E - pi.px*pi.px - pi.py*pi.py - pi.pz*pi.pz);
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

				if ( 	   abs(pi.px-Kx) > 0.5*px_bin_width
						or abs(pi.py-Ky) > 0.5*py_bin_width
						or abs(pi.pz-Kz) > 0.5*pz_bin_width
						or abs(pj.px-Kx) > 0.5*px_bin_width
						or abs(pj.py-Ky) > 0.5*py_bin_width
						or abs(pj.pz-Kz) > 0.5*pz_bin_width
					) continue;

				double bin_factor =
						  bin_function( Kx, pi.px, px_bin_width, bin_mode )
						* bin_function( Ky, pi.py, py_bin_width, bin_mode )
						* bin_function( Kz, pi.pz, pz_bin_width, bin_mode )
						* bin_function( Kx, pj.px, px_bin_width, bin_mode )
						* bin_function( Ky, pj.py, py_bin_width, bin_mode )
						* bin_function( Kz, pj.pz, pz_bin_width, bin_mode );

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

					double q0 = get_q0(m, qo, qs, ql, KT, KL);

					double arg =  q0 * (ti - tj)
								- qx * (xi - xj)
								- qy * (yi - yj)
								- qz * (zi - zj);

					double term = exp(i*arg/hbarC) * bin_factor;
					int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

					private_numerator[index6D] += term;
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
				in_numerator[idx] += private_numerator[idx];
				++idx;
			}
		}

	}

	return;
}


