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
	cout << "Using random number generator for toy model calculation!" << endl;

	allEvents.clear();

	double mass 	= paraRdr->getVal("mass");
	double RNG_R 	= paraRdr->getVal("RNG_R");
	double RNG_a 	= paraRdr->getVal("RNG_a");

	int RNG_Nev 	= paraRdr->getVal("RNG_Nev");
	int RNG_mult 	= paraRdr->getVal("RNG_mult");
	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	//default_random_engine generator;
	normal_distribution<double> distribution(0.0,RNG_R/sqrt(2.0));

	//allEvents.resize( RNG_Nev );

	// this toy function uses the model of
	// Zhang, Wiedemann, Slotta, and Heinz (1997)
	for (int iEvent = 0; iEvent < RNG_Nev; ++iEvent)
	{
		EventRecord event;
		//event.particles.resize( RNG_mult );

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
			//event.particles[iParticle] = particle;

		}

		allEvents.push_back( event );
		//allEvents[iEvent] = event;

	}

	return;
}




#endif
