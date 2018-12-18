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

# Clean out previous results
RESULTS_DIRECTORY=$HOME_DIRECTORY/results
if [ -d "$RESULTS_DIRECTORY" ]; then
	rm -rf $RESULTS_DIRECTORY
fi

mkdir $RESULTS_DIRECTORY

#===================
# Main calculation
#===================
runPythia=true

projectile="Pb"
target="Pb"
beamEnergy="2760.0"
Nevents="1000"
centralityCutString="0-10%"
centralityCut=(`echo $centralityCutString | sed 's/-/ /g' | sed 's/%//g'`)
lowerLimit=${centralityCut[0]}
upperLimit=${centralityCut[1]}

collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_C"$lowerLimit"_"$upperLimit"_Nev"$Nevents

#=====================================
# Run Pythia (if desired)
(

	cd $PYTHIA_DIRECTORY

	if $runPythia
	then

		# make sure results directory exists
		if [ ! -d "$PYTHIA_RESULTS_DIRECTORY" ]
		then
			    mkdir $PYTHIA_RESULTS_DIRECTORY
		fi

		# time and run
		./run_mainHIC.sh $projectile $target $beamEnergy \
							$lowerLimit $upperLimit $Nevents \
							"$PYTHIA_RESULTS_DIRECTORY/"
		#cp $PYTHIA_RESULTS_DIRECTORY/* $RESULTS_DIRECTORY/
	fi

	# Get the filenames which need to be processed
	recordOfOutputFilenames_Sxp=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_S_x_p_filenames.dat"
	recordOfOutputFilename_mult=$PYTHIA_RESULTS_DIRECTORY/`echo $collisionSystemStem`"_total_N_filename.dat"
	rm $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	for line in `cat $recordOfOutputFilenames_Sxp`
	do
		readlink -f $PYTHIA_RESULTS_DIRECTORY/$line >> $HBT_EVENT_GEN_DIRECTORY/catalogue.dat
	done
	# Set particle catalogue
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_EVENT_GEN_DIRECTORY/particle_catalogue.dat
	readlink -f $PYTHIA_RESULTS_DIRECTORY/HBT_particle.dat > $HBT_FITCF_DIRECTORY/particle_catalogue.dat
	# Set ensemble catalogue
	echo $projectile $target $beamEnergy $lowerLimit $upperLimit $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	readlink -f $recordOfOutputFilename_mult >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
	#don't need to do this directory
	#echo $projectile$target"_"`echo $beamEnergy`"GeV_C"$lowerLimit"_"$upperLimit"_Nev"$Nevents > $HBT_FITCF_DIRECTORY/ensemble_catalogue.dat
	#readlink -f $recordOfOutputFilename_mult >> $HBT_FITCF_DIRECTORY/ensemble_catalogue.dat

	# do sorting inside HBT code, not here
	#recordOfOutputFilename_multSorted=`cat $recordOfOutputFilename_mult | sed 's/total_N/total_N_Sorted/g'`
	#sort -k2,2 $PYTHIA_RESULTS_DIRECTORY/$recordOfOutputFilename_mult > $recordOfOutputFilename_multSorted

	#nev=`wc -l $recordOfOutputFilename_multSorted`
	#calcString=`echo 'print int(0.01*'$lowerLimit'*'$nev')'`
	#firstSortedEvent=`python -c $calcString`
	#calcString=`echo 'print int(0.01*'$upperLimit'*'$nev')'`
	#lastSortedEvent=`python -c $calcString`

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

	cp ../parameters.dat .

	# time and run
	nohup time ./run_HBT_event_generator.e \
			centrality_minimum=$lowerLimit \
			centrality_maximum=$upperLimit \
			1> HBT_event_generator.out \
			2> HBT_event_generator.err
	cp HBT_event_generator.[oe]* \
		./results/* $RESULTS_DIRECTORY/

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

	cp ../parameters.dat .

	# time and run
	nohup time ./run_fit_correlation_function.e \
			1> fit_correlation_function.out \
			2> fit_correlation_function.err
	cp fit_correlation_function.[oe]* ./results/* $RESULTS_DIRECTORY/

)

zip -r $HOME_DIRECTORY/`echo $collisionSystemStem`"_results.zip" $RESULTS_DIRECTORY

# Clean-up HBT directories (but not Pythia results directory!!!)
rm -rf $HBT_EVENT_GEN_DIRECTORY/*HBT_event_generator.[oe]* $HBT_EVENT_GEN_DIRECTORY/results\
       $HBT_FITCF_DIRECTORY/*fit_correlation_function.[oe]* $HBT_FITCF_DIRECTORY/results

echo 'Finished everything!'


# End of file
