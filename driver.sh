#! /usr/bin/env bash

#=====================================
# Header info
HOME_DIRECTORY=~/HBT_event_generator
# Pythia
PYTHIA_DIRECTORY=~/pythia8235/examples
ALT_PYTHIA_DIRECTORY=/scratch/blixen/plumberg
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=./HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=./fit_correlation_function

# Make sure results directory exists
RESULTS_DIRECTORY=$HOME_DIRECTORY/results
if [ ! -d "$RESULTS_DIRECTORY" ]; then
	mkdir $RESULTS_DIRECTORY
fi

#===================
# Main calculation
#===================
runPythia=true

#=====================================
# Run Pythia (if desired)
(

	cd $PYTHIA_DIRECTORY

	projectile="Pb"
	target="Pb"
	beamEnergy="2760.0"
	Nevents="100"

	if [ $runPythia ]
	then
		# time and run
		./run_mainHIC.sh $projectile $target $beamEnergy $Nevents
		cp ./results/* $RESULTS_DIRECTORY/
	fi

	readlink -f ./results/$projectile$target"_"$beamEnergy"GeV"_*.dat > $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	readlink -f ./results/HBT_particle.dat > $HBT_EVENT_GEN_DIRECTORY/particle_catalogue.dat
	readlink -f ./results/HBT_particle.dat > $HBT_FITCF_DIRECTORY/particle_catalogue.dat

)

#=====================================
# Run HBT_event_generator
(

	cd $HBT_EVENT_GEN_DIRECTORY

	# using OpenMP (leave a couple cores free)
	export OMP_NUM_THREADS=10

	# time and run
	nohup time ./HBT_event_generator.e 1>> HBT_event_generator.out 2>> HBT_event_generator.err
	cp ./results/* $RESULTS_DIRECTORY/

	readlink -f ./results/HBT_pipiCF_*.dat > $HBT_EVENT_GEN_DIRECTORY/catalogue.dat

)

#=====================================
# Run fit_correlation_function
(

	cd $HBT_FITCF_DIRECTORY

	# time and run
	nohup time ./run.e 1> run.out 2> run.err
	cp ./results/* $RESULTS_DIRECTORY/

)

zip -r ./$projectile$target"_"`echo $beamEnergy`GeV_results.zip ./$RESULTS_DIRECTORY

echo 'Finished everything!'


# End of file
