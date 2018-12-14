#! /usr/bin/env bash

#=====================================
# Header info
HOME_DIRECTORY=~/HBT_event_generator
# Pythia
PYTHIA_DIRECTORY=~/pythia8235/examples
ALT_PYTHIA_DIRECTORY=/scratch/blixen/plumberg
PYTHIA_RESULTS_DIRECTORY=$ALT_PYTHIA_DIRECTORY/results
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HOME_DIRECTORY/fit_correlation_function

# Make sure results directory exists
RESULTS_DIRECTORY=$HOME_DIRECTORY/results
if [ ! -d "$RESULTS_DIRECTORY" ]; then
	mkdir $RESULTS_DIRECTORY
fi

#===================
# Main calculation
#===================
runPythia=true

projectile="Pb"
target="Pb"
beamEnergy="2760.0"
Nevents="100000"

#=====================================
# Run Pythia (if desired)
(

	cd $PYTHIA_DIRECTORY

	if [ $runPythia ]
	then

		# make sure results directory exists
		if [ ! -d "$PYTHIA_RESULTS_DIRECTORY" ]
		then
			    mkdir $PYTHIA_RESULTS_DIRECTORY
		fi

		# time and run
		./run_mainHIC.sh $projectile $target $beamEnergy $Nevents "$PYTHIA_RESULTS_DIRECTORY/"
		#cp $PYTHIA_RESULTS_DIRECTORY/* $RESULTS_DIRECTORY/
	fi

	# Get the filenames which need to be processed
	recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents"_S_x_p_filenames.dat"
	recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents"_total_N_filename.dat"
	for line in `cat $recordOfOutputFilenames_Sxp`
	do
		readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	done
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_EVENT_GEN_DIRECTORY/particle_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_FITCF_DIRECTORY/particle_catalogue.dat

)

#=====================================
# Run HBT_event_generator
(

	cd $HBT_EVENT_GEN_DIRECTORY

	# using OpenMP (leave a couple cores free)
	export OMP_NUM_THREADS=10

	if [ ! -d "./results" ]; then
	        mkdir results
	fi

	# time and run
	nohup time ./run.e 1>> HBT_event_generator.out 2>> HBT_event_generator.err
	cp HBT_event_generator.[oe]* ./results/* $RESULTS_DIRECTORY/

	readlink -f ./results/HBT_pipiCF.dat > $HBT_FITCF_DIRECTORY/catalogue.dat

)

#=====================================
# Run fit_correlation_function
(

	cd $HBT_FITCF_DIRECTORY

	if [ ! -d "./results" ]
	then
		mkdir results
	fi

	# time and run
	nohup time ./run.e 1> fit_correlation_function.out 2> fit_correlation_function.err
	cp fit_correlation_function.[oe]* ./results/* $RESULTS_DIRECTORY/

)

zip -r $HOME_DIRECTORY/$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents"_results.zip" $RESULTS_DIRECTORY

# Clean-up HBT directories (but not Pythia results directory!!!)
rm -rf $HBT_EVENT_GEN_DIRECTORY/HBT_event_generator.[oe]* $HBT_EVENT_GEN_DIRECTORY/results\
       $HBT_FITCF_DIRECTORY/fit_correlation_function.[oe]* $HBT_FITCF_DIRECTORY/results

echo 'Finished everything!'


# End of file
