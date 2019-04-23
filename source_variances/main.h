#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "src/ParameterReader.h"
#include "src/ensemble.h"

using namespace std;

//this is just to give this file a reason to exist for the moment...
const double plumbergtest = 0.;
const bool CONVERT_MM_TO_FM = true;	// needs to be true if running on Pythia output, false for Vishnu output
const double MmPerFm = ( CONVERT_MM_TO_FM ) ? 1.e-12 : 1.0;	//mm-to-fm conversion

vector<EventMultiplicity> ensemble_multiplicites;

void display_intro(int logo_mode)
{
	cout << endl
			<< "              HBT event generator              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 10/2018" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo( logo_mode ); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
}

string get_filename( string path, string chosen_particle_name, string extension )
{
	ostringstream filename_stream;
	filename_stream << path << "HBT_"
					<< chosen_particle_name << chosen_particle_name
					<< "SVradii." << extension;
	return ( filename_stream.str() );
}



#endif
