#ifndef RANDOM_EVENTS_H
#define RANDOM_EVENTS_H

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>

#include "./EventRecord.h"
#include "./ParticleRecord.h"
#include "./ParameterReader.h"

using namespace std;

void generate_events(vector<EventRecord> & allEvents, ParameterReader * paraRdr)
{
	//cout << "Using random number generator for toy model calculation!" << endl;

	allEvents.clear();

	double mass 	= paraRdr->getVal("mass");
	double RNG_R 	= paraRdr->getVal("RNG_R");
	double RNG_a 	= paraRdr->getVal("RNG_a");

	int RNG_Nev 	= paraRdr->getVal("RNG_Nev");
	int RNG_mult 	= paraRdr->getVal("RNG_mult");
	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	bool RNG_seed	= (bool)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed )
		generator = default_random_engine (seed);

	normal_distribution<double> distribution(0.0,RNG_R/sqrt(2.0));

	// this toy function uses the model of
	// Zhang, Wiedemann, Slotta, and Heinz (1997)
	for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
	{
		EventRecord event;

		for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
		{
			
			double tP = 0.0;	// cf. paper
			double xP = RNG_xDir ? distribution(generator) : 0.0;
			double yP = RNG_yDir ? distribution(generator) : 0.0;
			double zP = RNG_zDir ? distribution(generator) : 0.0;

			double px = RNG_a * xP;
			double py = RNG_a * yP;
			double pz = RNG_a * zP;
			double Ep = sqrt( mass*mass + px*px + py*py + pz*pz );

			ParticleRecord particle;
			particle.eventID 	= iEvent;
			particle.particleID = iParticle;
			particle.E 			= Ep;
			particle.px 		= px;
			particle.py 		= py;
			particle.pz 		= pz;
			particle.t 			= tP;
			particle.x 			= xP;
			particle.y 			= yP;
			particle.z 			= zP;

			event.particles.push_back( particle );

		} 

		allEvents.push_back( event );

	}

	return;
}


void generate_events_v2(vector<EventRecord> & allEvents, ParameterReader * paraRdr)
{
	//cout << "Using random number generator for toy model calculation!" << endl;

	allEvents.clear();

	double mass 	= paraRdr->getVal("mass");
	double RNG_R 	= paraRdr->getVal("RNG_R");
	double RNG_a 	= paraRdr->getVal("RNG_a");

	int RNG_Nev 	= paraRdr->getVal("RNG_Nev");
	int RNG_mult 	= paraRdr->getVal("RNG_mult");
	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	bool RNG_seed	= (bool)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed )
		generator = default_random_engine (seed);

	normal_distribution<double> distribution(0.0,RNG_R/sqrt(2.0));

	int model_to_use = 2;
	// model_to_use values:
	//		0 - Zhang, Wiedemann, Slotta, and Heinz (1997) or slight variant
	//		1 - my own toy emission function
	//		2 - another of my own toy emission functions

	if ( model_to_use == 0 )
	{

		normal_distribution<double> distribution( 0.0, RNG_R/sqrt(2.0) );

		for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
		{
			EventRecord event;

			for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
			{
			
				double tP = 0.0;	// cf. paper
				double xP = RNG_xDir ? distribution(generator) : 0.0;
				double yP = RNG_yDir ? distribution(generator) : 0.0;
				double zP = RNG_zDir ? distribution(generator) : 0.0;

				double px = RNG_a * xP;
				double py = RNG_a * yP;
				double pz = RNG_a * zP;
				double Ep = sqrt( mass*mass + px*px + py*py + pz*pz );

				ParticleRecord particle;
				particle.eventID 	= iEvent;
				particle.particleID = iParticle;
				particle.E 			= Ep;
				particle.px 		= px;
				particle.py 		= py;
				particle.pz 		= pz;
				particle.t 			= tP;
				particle.x 			= xP;
				particle.y 			= yP;
				particle.z 			= zP;

				event.particles.push_back( particle );

			} 

			allEvents.push_back( event );

		}
	}
	else if ( model_to_use == 1 and RNG_xDir
				and RNG_yDir and RNG_zDir )								// 3D only for now
	{

		const double TFO = 0.12;										// FO temp of 120 MeV

		normal_distribution<double> distribution( 0.0, 1.0/sqrt(2.0) );	// rescale by variable source radius
		uniform_real_distribution<double> KPhi_distribution( 0.0, 1.0 );
		exponential_distribution<double> KT_distribution( TFO );		// set width from FO temp

		// set some emission function parameters
		const double t0 = 7.5, z0 = 0.0, R0 = 3.0, Del_t = 2.0, Del_z = 2.5;

		// my toy emission function
		for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
		{
			EventRecord event;

			for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
			{
				
				// sample KT, KL first
				double KT = KT_distribution(generator);
				double KL = 0.0;	// for simplicity

				double Rscale = R0 / sqrt(1.0 + 0.1*sqrt(mass*mass+KT*KT)/TFO);
				//Rscale = R0;

				double tP = t0 + Del_t * distribution(generator);
				double xP = Rscale * distribution(generator);
				double yP = Rscale * distribution(generator);
				double zP = z0 + Del_z * distribution(generator);

				double phiT = atan2(yP, xP);

				//double Kphi = 2.0*M_PI*KPhi_distribution(generator);
				double Kphi = phiT;

				double px = KT * cos(Kphi);
				double py = KT * sin(Kphi);
				double pz = 0.0;

				double Ep = sqrt( mass*mass + px*px + py*py + pz*pz );

				ParticleRecord particle;
				particle.eventID 	= iEvent;
				particle.particleID = iParticle;
				particle.E 			= Ep;
				particle.px 		= px;
				particle.py 		= py;
				particle.pz 		= pz;
				particle.t 			= tP;
				particle.x 			= xP;
				particle.y 			= yP;
				particle.z 			= zP;

				event.particles.push_back( particle );

			} 

			allEvents.push_back( event );

		}
	}
	else if ( model_to_use == 2 )
	{

		// test sensitivity to bin-width
		double bin_epsilon = paraRdr->getVal("bin_epsilon");
		double KLmin = paraRdr->getVal("KLmin");
		double KLmax = paraRdr->getVal("KLmax");
		double KTmin = paraRdr->getVal("KTmin");
		double KTmax = paraRdr->getVal("KTmax");

		const double TFO = 0.12;										// FO temp of 120 MeV

		normal_distribution<double> distribution( 0.0, 1.0/sqrt(2.0) );	// rescale by variable source radius
		//uniform_real_distribution<double> KPhi_distribution( 0.0, 2.0*M_PI );
		//exponential_distribution<double> KT_distribution( 1.0 / TFO );		// set width from FO temp
		//uniform_real_distribution<double> KT_distribution( KTmin, KTmax );
		//weibull_distribution<double> weird_KT_distribution( 2.0, sqrt(2.0)*TFO );			// Weibull corresponds to distribution
																				// of magnitude of two normal R.V.s
		//uniform_real_distribution<double> KL_distribution( KLmin, KLmax );

		// set some emission function parameters
		const double t0 = 7.5, z0 = 0.0, R0 = 3.0, Del_t = 2.0, Del_z = 2.5;

		// my toy emission function
		for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
		{
			EventRecord event;

			for (int iParticle = 0; iParticle < RNG_mult; ++iParticle)
			{
				
				// sample KT, KL first
				// just use uniform K distribution for simplicity
				//double KT = KT_distribution(generator);
				double KT = 0.0;
				//double KL = KL_distribution(generator);
				//double KL = 0.0;
				//double KL = KLmax * distribution(generator);
				//double KL = bin_epsilon * distribution(generator);

				double Rscale = R0 / sqrt(1.0 + 0.1*sqrt(mass*mass+KT*KT)/TFO);
				Rscale = R0;

				//double tP = t0 + Del_t * distribution(generator);
				double tP = 0.0;
				double xP = RNG_xDir ? Rscale * distribution(generator) : 0.0;
				double yP = RNG_yDir ? Rscale * distribution(generator) : 0.0;
				double zP = RNG_zDir ? Rscale * distribution(generator) : 0.0;

				//double phiT = atan2(yP, xP);

				//double Kphi = KPhi_distribution(generator);
				//double Kphi = phiT;
				//double Kphi = 0.0;

				/*do
				{
					KT = sqrt(2.0) * TFO * distribution(generator);
				} while ( KT < 0.0 );*/

				//KT = weird_KT_distribution(generator);

				//double px = KT * cos(Kphi);
				double px = TFO * distribution(generator);
				//double py = KT * sin(Kphi);
				double py = TFO * distribution(generator);
				//double pz = KL;
				double pz = TFO * distribution(generator);
				//pz = 0.0;

				double Ep = sqrt( mass*mass + px*px + py*py + pz*pz );

				ParticleRecord particle;
				particle.eventID 	= iEvent;
				particle.particleID = iParticle;
				particle.E 			= Ep;
				particle.px 		= px;
				particle.py 		= py;
				particle.pz 		= pz;
				particle.t 			= tP;
				particle.x 			= xP;
				particle.y 			= yP;
				particle.z 			= zP;
//cout << "particle check: " << Ep << "   " << px << "   " << py << "   " << pz << "   " << tP << "   " << xP << "   " << yP << "   " << zP << endl;

				event.particles.push_back( particle );

			} 

			allEvents.push_back( event );

		}
	}
	else
	{
		cout << "CRASH!!!" << endl;
		exit(10);
	}

	return;
}

#endif
