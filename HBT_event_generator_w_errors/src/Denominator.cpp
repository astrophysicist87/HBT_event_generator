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


void HBT_event_generator::Compute_denominator( vector<double> & in_denominator )
{
	// Sum over all events
	#pragma omp parallel for schedule(static) shared(in_denominator)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_denominator(in_denominator.size());

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			if (iParticle == jParticle)
				continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
				double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
				double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				double cKphi = cos(Kphi), sKphi = sin(Kphi);
				double Kx = KT * cKphi;
				double Ky = KT * sKphi;
				double Kz = KL;
				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				double p1x = Kx + 0.5 * qx, p1y = Ky + 0.5 * qy, p1z = Kz + 0.5 * qz;
				double p2x = Kx - 0.5 * qx, p2y = Ky - 0.5 * qy, p2z = Kz - 0.5 * qz;

				
				private_denominator[indexer(iKT, iKphi, iKL, iqo, iqs, iql)] += 
						  double(abs(p1x-pi.px) <= 0.5*px_bin_width)
						* double(abs(p1y-pi.py) <= 0.5*py_bin_width)
						* double(abs(p1z-pi.pz) <= 0.5*pz_bin_width)
						* double(abs(p2x-pj.px) <= 0.5*px_bin_width)
						* double(abs(p2y-pj.py) <= 0.5*py_bin_width)
						* double(abs(p2z-pj.pz) <= 0.5*pz_bin_width)
						/ ( px_bin_width * py_bin_width * pz_bin_width
							* px_bin_width * py_bin_width * pz_bin_width );


				denominator_cell_was_filled[indexer(iKT, iKphi, iKL, iqo, iqs, iql)]
						= ( abs(pi.px-p1x) > 0.5*px_bin_width
							or abs(pi.py-p1y) > 0.5*py_bin_width
							or abs(pi.pz-p1z) > 0.5*pz_bin_width
							or abs(pj.px-p2x) > 0.5*px_bin_width
							or abs(pj.py-p2y) > 0.5*py_bin_width
							or abs(pj.pz-p2z) > 0.5*pz_bin_width );

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
				in_denominator[idx] += private_denominator[idx];
				++idx;
			}
		}

	}

	return;
}


