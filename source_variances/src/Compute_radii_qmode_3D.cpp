#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstdlib>

#include "SourceVariances.h"
#include "matrix.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

void SourceVariances::Compute_radii_qmode_3D()
{
	// Sum over all events
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{

			ParticleRecord p = event.particles[iParticle];

			double t = p.t, x = p.x, y = p.y, z = p.z;
			double E = p.E, px = p.px, py = p.py, pz = p.pz;

			// New method of binning
			double K0 = E, Kx = px, Ky = py, Kz = pz;

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
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

/*cout << "particles which passed: "
		<< iEvent << "   " << iParticle << "   "
		<< E << "   " << px << "   " << py << "   " << pz << "   "
		<< t << "   " << x << "   " << y << "   " << z << endl;*/

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);

			//check this
			//double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double cKphi = cos(Kphi_bins[Kphi_idx]), sKphi = sin(Kphi_bins[Kphi_idx]);
			double xo = x * cKphi + y * sKphi;
			double xs = y * cKphi - x * sKphi;
			double xl = z;

			S[index3D] 		 += 1.0;
			x_S[index3D]	 += x;
			x2_S[index3D]	 += x*x;
			y_S[index3D]	 += y;
			y2_S[index3D]	 += y*y;
			z_S[index3D]	 += z;
			z2_S[index3D]	 += z*z;
			xo_S[index3D] 	 += xo;
			xo2_S[index3D] 	 += xo*xo;
			xs_S[index3D] 	 += xs;
			xs2_S[index3D] 	 += xs*xs;
			xl_S[index3D] 	 += xl;
			xl2_S[index3D] 	 += xl*xl;
			t_S[index3D] 	 += t;
			t2_S[index3D] 	 += t*t;
			x_t_S[index3D] 	 += x*t;
			y_t_S[index3D] 	 += y*t;
			z_t_S[index3D] 	 += z*t;
			x_y_S[index3D] 	 += x*y;
			x_z_S[index3D] 	 += x*z;
			y_z_S[index3D] 	 += y*z;
			xo_t_S[index3D]  += xo*t;
			xs_t_S[index3D]  += xs*t;
			xl_t_S[index3D]  += xl*t;
			xo_xs_S[index3D] += xo*xs;
			xo_xl_S[index3D] += xo*xl;
			xs_xl_S[index3D] += xs*xl;

		}

	}

}


void SourceVariances::Average_source_moments()
{
	int index3D = 0;
	for (int iKT = 0; iKT < n_KT_bins; ++iKT)
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
	for (int iKL = 0; iKL < n_KL_bins; ++iKL)
	{
		double norm = S[index3D];

		x_S[index3D]	 /= norm;
		x2_S[index3D]	 /= norm;
		y_S[index3D]	 /= norm;
		y2_S[index3D]	 /= norm;
		z_S[index3D]	 /= norm;
		z2_S[index3D]	 /= norm;
		xo_S[index3D] 	 /= norm;
		xo2_S[index3D] 	 /= norm;
		xs_S[index3D] 	 /= norm;
		xs2_S[index3D] 	 /= norm;
		xl_S[index3D] 	 /= norm;
		xl2_S[index3D] 	 /= norm;
		t_S[index3D] 	 /= norm;
		t2_S[index3D] 	 /= norm;
		x_t_S[index3D] 	 /= norm;
		y_t_S[index3D] 	 /= norm;
		z_t_S[index3D] 	 /= norm;
		x_y_S[index3D] 	 /= norm;
		x_z_S[index3D] 	 /= norm;
		y_z_S[index3D] 	 /= norm;
		xo_t_S[index3D]  /= norm;
		xs_t_S[index3D]  /= norm;
		xl_t_S[index3D]  /= norm;
		xo_xs_S[index3D] /= norm;
		xo_xl_S[index3D] /= norm;
		xs_xl_S[index3D] /= norm;

		++index3D;

	}	
}


void SourceVariances::Set_radii()
{

	int idx3D = 0;
	double mass = particle_mass;

	// Get indices
	for (int iKT = 0; iKT < n_KT_bins; ++iKT)
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
	for (int iKL = 0; iKL < n_KL_bins; ++iKL)
	{

		double KT_loc = KT_bins[iKT];
		double KL_loc = KL_bins[iKL];

		double bT 	= KT_loc
						/ sqrt(mass*mass
							+ KT_loc*KT_loc
							+ KL_loc*KL_loc);
		double bL 	= KL_loc
						/ sqrt(mass*mass
							+ KT_loc*KT_loc
							+ KL_loc*KL_loc);

		double xo 	 = xo_S[idx3D];
		double xo2 	 = xo2_S[idx3D];
		double xs 	 = xs_S[idx3D];
		double xs2 	 = xs2_S[idx3D];
		double xl 	 = xl_S[idx3D];
		double xl2 	 = xl2_S[idx3D];
		double t 	 = t_S[idx3D];
		double t2 	 = t2_S[idx3D];
		double xo_t  = xo_t_S[idx3D];
		double xs_t  = xs_t_S[idx3D];
		double xl_t  = xl_t_S[idx3D];
		double xo_xs = xo_xs_S[idx3D];
		double xo_xl = xo_xl_S[idx3D];
		double xs_xl = xs_xl_S[idx3D];

		R2o[idx3D] = ( xo2 - xo*xo )
						-2.0*bT*( xo_t - xo*t )
						+bT*bT*( t2 - t*t );

		R2s[idx3D] = xs2 - xs*xs;

		R2l[idx3D] = ( xl2 - xl*xl )
						-2.0*bL*( xl_t - xl*t )
						+bL*bL*( t2 - t*t );

		R2os[idx3D] = ( xo_xs - xo*xs )
						-bT*( xs_t - xs*t );

		R2ol[idx3D] = ( xo_xl - xo*xl )
						-bT*( xl_t - xl*t )
						-bL*( xo_t - xo*t )
						+bT*bL*( t2 - t*t );

		R2sl[idx3D] = ( xs_xl - xs*xl )
						-bL*( xs_t - xs*t );

		double HBT2D[4] = { R2o[idx3D], R2os[idx3D],
							R2os[idx3D], R2s[idx3D] };
		double HBT3D[9] = { R2o[idx3D],  R2os[idx3D], R2ol[idx3D],
							R2os[idx3D], R2s[idx3D],  R2sl[idx3D],
							R2ol[idx3D], R2sl[idx3D], R2l[idx3D] };

		//2D and 3D "HBT volumes"
		detHBT2D[idx3D] = det2D( HBT2D );
		detHBT3D[idx3D] = det3D( HBT3D );

		++idx3D;

	}

}


//End of file
