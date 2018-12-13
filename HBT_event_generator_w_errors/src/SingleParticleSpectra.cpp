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



void HBT_event_generator::Compute_spectra()
{
	Compute_dN_pTdpTdpphidpY();

	Compute_dN_2pipTdpTdpY();
	Compute_dN_pTdpTdpphi();
	Compute_dN_dpphidpY();

	Compute_dN_2pipTdpT();
	Compute_dN_dpphi();
	Compute_dN_2pidpY();

	Compute_N();

	return;
}



void HBT_event_generator::Compute_dN_pTdpTdpphidpY()
{
	// Sum over all events
	int NEvents = allEvents.size();
	for (int iEvent = 0; iEvent < NEvents; ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];

			int ipT 	= bin_function(p.pT, pT_pts);
			int ipphi 	= bin_function(p.pphi, pphi_pts);
			int ipY 	= bin_function(p.pY, pY_pts);

			if ( ipT < 0 or ipT >= n_pT_bins
					or ipphi < 0 or ipphi >= n_pphi_bins
					or ipY < 0 or ipY >= n_pY_bins
				) 	{
						//err << "Warning: discard particle with (pT, pphi, pY) = " << p.pT << "   " << p.pphi << "   " << p.pY << endl;
						continue;
					}

			dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)]+=1.0;
		}
	}

	// Average over events
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);

		dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)]
				/= ( pT_bin_center * double(NEvents)
						* pT_bin_width * pphi_bin_width * pY_bin_width );
	}

	return;
}


/////////////
void HBT_event_generator::Compute_dN_2pipTdpTdpY()
{
	vector<double> dN_2pipTdpTdpY(n_pT_bins*n_pY_bins);

	int idx2D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		dN_2pipTdpTdpY[idx2D] = 0.0;

		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
			dN_2pipTdpTdpY[idx2D] += dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];

		dN_2pipTdpTdpY[idx2D++] *= pphi_bin_width / ( 2.0 * M_PI );
	}

	return;
}



/////////////
void HBT_event_generator::Compute_dN_pTdpTdpphi()
{
	vector<double> dN_pTdpTdpphi(n_pT_bins*n_pphi_bins);

	int idx2D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	{
		dN_pTdpTdpphi[idx2D] = 0.0;

		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
			dN_pTdpTdpphi[idx2D] += dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];

		dN_pTdpTdpphi[idx2D++] *= pY_bin_width;
	}

	return;
}

void HBT_event_generator::Compute_dN_dpphidpY()
{
	vector<double> dN_dpphidpY(n_pT_bins*n_pY_bins);

	int idx2D = 0;
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		dN_dpphidpY[idx2D] = 0.0;

		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		{
			double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
			dN_dpphidpY[idx2D] += pT_bin_center
									* dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];
		}

		dN_dpphidpY[idx2D++] *= pT_bin_width;
	}

	return;
}

void HBT_event_generator::Compute_dN_2pipTdpT()
{
	vector<double> dN_2pipTdpT(n_pT_bins);

	int idx1D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	{
		dN_2pipTdpT[idx1D] = 0.0;

		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);

		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
			dN_2pipTdpT[idx1D] += dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];

		dN_2pipTdpT[idx1D++] *= pY_bin_width * pphi_bin_width / ( 2.0 * M_PI );
	}

	//for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	//	out << 0.5*(pT_pts[ipT]+pT_pts[ipT+1]) << "   " << dN_2pipTdpT[ipT] << endl;

	return;
}

void HBT_event_generator::Compute_dN_dpphi()
{
	vector<double> dN_dpphi(n_pphi_bins);

	int idx1D = 0;
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	{
		dN_dpphi[idx1D] = 0.0;

		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
			dN_dpphi[idx1D] += pT_bin_center
								* dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];
		}

		dN_dpphi[idx1D++] *= pT_bin_width * pY_bin_width;
	}

	return;
}

void HBT_event_generator::Compute_dN_2pidpY()
{
	vector<double> dN_2pidpY(n_pY_bins);

	int idx1D = 0;
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		dN_2pidpY[idx1D] = 0.0;

		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		{
			double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
			dN_2pidpY[idx1D] += pT_bin_center
								* dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];
		}

		dN_2pidpY[idx1D++] *= pT_bin_width * pphi_bin_width / ( 2.0 * M_PI );
	}

	return;
}



void HBT_event_generator::Compute_N()
{
	double N = 0.0;
	int NEvents = allEvents.size();

	int idx3D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
		N += pT_bin_center * dN_pTdpTdpphidpY[idx3D++]
				* pT_bin_width * pphi_bin_width * pY_bin_width;
	}

	out << "Event-averaged multiplicity <N> = " << N << endl;
	out << "Total multiplicity N = " << NEvents*N << endl;

	return;
}



