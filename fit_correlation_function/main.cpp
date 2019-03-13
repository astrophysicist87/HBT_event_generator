#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <sys/time.h>
#include <fenv.h>

#include "src/Stopwatch.h"
#include "src/correlation_function.h"
#include "src/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
	//===============
	// Display intro
	print_header(2);   

	//===================================
	// Read-in free parameters
	ParameterReader * paraRdr
		= new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	// Read-in HBT particle information
	string particle_info_filename
		= read_file_catalogue("./particle_catalogue.dat");
	paraRdr->readFromFile(particle_info_filename);

	// Read-in command-line arguments
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();

	// Read in name of file(s) to process
	string input_filename
		= read_file_catalogue("./catalogue.dat");

	//===================================
	// Main calculation starts here
	//===================================
	// Start timing
	Stopwatch sw;
	sw.Start();

	// Throws exception if NaNs are encountered
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	// Set-up output files
	string path
		= dirname( input_filename );	// make sure this directory exists
	string base_filename
		= basename( input_filename );
	ostringstream out_filename_stream, err_filename_stream;
	out_filename_stream << path << "fit_"
						<< change_file_extension(
								base_filename,
								"dat", "out" );
	err_filename_stream << path << "fit_"
						<< change_file_extension(
								base_filename,
								"dat", "err" );

	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());

	cout << "fit_correlation_function(): set-up completed." << endl;

	// Define correlation function object
	// (loading and fitting correlation function
	//  performed automatically)
	Correlation_function
		HBT_event_ensemble( paraRdr, input_filename,
							outmain, errmain );
	
	string HBTradii_filename
		= "./results/HBTradii.dat";
	HBT_event_ensemble.Output_HBTradii( HBTradii_filename );


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;

	// Wrap it up!
	return (0);
}

//End of file
