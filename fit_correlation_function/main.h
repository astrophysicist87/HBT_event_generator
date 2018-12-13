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

using namespace std;

//this is just to give this file a reason to exist for the moment...
const double plumbergtest = 0.;
const double MmPerFm = 1.e-12;	//mm/fm

void print_header(int logostyle)
{
	cout << endl
			<< "              Fit correlation function              " << endl
			<< endl
			<< "  Ver 1.0   ----- (c) Christopher Plumberg, 11/2018" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo(logostyle); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;

	return;
}

// function to read in catalogue containing name
// of correlation function file to be processed
string read_file_catalogue(string catalogue_name)
{
	ifstream catalogue_in(catalogue_name.c_str());

	string filename;
	getline(catalogue_in, filename);

	return ( filename );
}


// returns directory component of
// full file pathname
string dirname(string pathname)
{
    string buf;					// Have a buffer string
    stringstream ss(pathname);	// Insert the string into a stream
    vector<string> tokens;		// Create vector to hold our words

	// create vector containing each
	// 'chunk' of full path
	while(getline(ss, buf, '/'))
		tokens.push_back(buf);

	// get rid of the filename itself
	tokens.pop_back();

	// concatenate everything else back
	// into directory name
	string dir = "";
	for( vector<string>::iterator it = tokens.begin();
			it != tokens.end();
			++it)
		dir += (*it) + '/';

	return ( dir );
}

// returns basename component of
// full file pathname
string basename(string pathname)
{
    string buf;					// Have a buffer string
    stringstream ss(pathname);	// Insert the string into a stream
    vector<string> tokens;		// Create vector to hold our words

	while(getline(ss, buf, '/'))
		tokens.push_back(buf);

	return ( tokens.back() );
}

// returns filename with new extension
string change_file_extension( string old_filename,
								string old_file_extension,
								string new_file_extension )
{
	int OF_len = old_filename.length();
	int OFE_len = old_file_extension.length();
	return ( old_filename.erase( OF_len - OFE_len )
				+ new_file_extension );
}

#endif
