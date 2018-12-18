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
#include "src/HBT_event_generator.h"
#include "src/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Display intro
	cout << endl
			<< "              HBT event generator              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 10/2018" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo(2); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
   

	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	// Read-in particle and ensemble information
	vector<string> particle_info_filename, ensemble_info_filename;
	read_file_catalogue("./particle_catalogue.dat", particle_info_filename);
	paraRdr->readFromFile(particle_info_filename[0]);

	// Read-in command-line arguments
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();

	// Start timing
	Stopwatch sw;
	sw.Start();

	// Throws exception if NaNs are encountered
	feenableexcept(FE_INVALID | FE_OVERFLOW);


	// Specify files containing all position-momentum information
	// from which to construct HBT correlation function
	vector<string> all_file_names;
	read_file_catalogue("./catalogue.dat", all_file_names);


	// Process multiplicity and ensemble information
	vector<string> ensemble_info;
	read_file_catalogue("./ensemble_catalogue.dat", ensemble_info);

	// allows to give files appropriate names
	string collision_system_info = ensemble_info[0];
	string target_name, projectile_name, beam_energy;
	int Nevents;
	double centrality_minimum, centrality_maximum;
	istringstream iss(collision_system_info);
	iss >> target_name >> projectile_name >> beam_energy >> centrality_minimum >> centrality_maximum >> Nevents;
	
	// select only those events falling into specificed centrality range
	string multiplicity_filename = ensemble_info[1];
	get_events_in_centrality_class(
				multiplicity_filename, ensemble_multiplicites,
				centrality_minimum, centrality_maximum );

	cout << "run_HBT_event_generator(): "
			<< "Using " << ensemble_multiplicites.size()
			<< " events in centrality class "
			<< centrality_minimum << "-"
			<< centrality_maximum << "%!" << endl;


	// Set-up output files
	string path = "./results/";	// make sure this directory exists
	string chosen_particle_name = "pi";
	//string collision_system = "pp_13TeV";
	//string collision_system = "pPb_5.02TeV";
	string collision_system = "PbPb_2.760TeV";
	ostringstream out_filename_stream, err_filename_stream;
	/*out_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF_" << collision_system << "_100000events.dat";
	err_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF_" << collision_system << "_100000events.err";*/
	out_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.dat";
	err_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.err";
	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());


	// Proceed with HBT calculations
	string mode = "read stream";

	if ( mode == "default" or mode == "read stream" )
	{

		// Vector to hold all event information
		vector<EventRecord> allEvents;


		// Read in the first file
		int iFile = 0;
		cout << "Processing " << all_file_names[iFile] << "..." << endl;

		get_all_events(all_file_names[iFile], allEvents, paraRdr);


		// Create HBT_event_generator object here
		// note: numerator and denominator computed automatically
		HBT_event_generator
			HBT_event_ensemble( paraRdr, allEvents,
								outmain, errmain );


		// Loop over the rest of the files
		for (iFile = 1; iFile < all_file_names.size(); ++iFile)
		{

			cout << "Processing " << all_file_names[iFile] << "..." << endl;

			// Read in the next file
			get_all_events(all_file_names[iFile], allEvents, paraRdr);


			// - for each file, update numerator and denominator
			HBT_event_ensemble.Update_records( allEvents );

		}

		// Compute correlation function itself (after
		// all events have been read in)
		HBT_event_ensemble.Compute_correlation_function();


		// Output results
		HBT_event_ensemble.Output_correlation_function();

	}
	else if ( mode == "read all" )
	{

		// Vector to hold all event information
		vector<EventRecord> allEvents;


		// Read in the files
		get_all_events(all_file_names, allEvents, paraRdr);


		// Create HBT_event_generator object from allEvents
		HBT_event_generator
			HBT_event_ensemble( paraRdr, allEvents,
								outmain, errmain );


		// Compute correlation function itself
		HBT_event_ensemble.Compute_correlation_function();


		// Output correlation function
		HBT_event_ensemble.Output_correlation_function();

	}


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}

//End of file
