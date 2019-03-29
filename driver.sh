#! /usr/bin/env bash

#=====================================
# Header info
HOME_DIRECTORY=~/HBT_event_generator
# Pythia
PYTHIA_DIRECTORY=~/pythia8235/examples
ALT_PYTHIA_DIRECTORY=/scratch/blixen/plumberg
PYTHIA_RESULTS_DIRECTORY=$ALT_PYTHIA_DIRECTORY/check_old_BEeffects
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HOME_DIRECTORY/fit_correlation_function

#============================
# Some function definitions
#============================
check_success () {
	if [ "$2" -eq "0" ]
	then
		echo '===================================================='
		echo '== '$1': run completed successfully!'
		echo '===================================================='
	else
		echo '===================================================='
		echo '== '$1': problems encountered!'
		echo '===================================================='
		exit $2
	fi
} 

#===================
# Main calculation
#===================
runPythia=true

projectile="p"
target="p"
beamEnergy="13000.0"	#GeV
Nevents="100000"
QRefValue="0.2"			#GeV
#PythiaExecutableToUse=""

echo 'Processing Nevents =' $Nevents $projectile'+'$target 'collisions at' $beamEnergy 'GeV'

nCC=0
for centralityCutString in "0-100%"
do
	success=0
	echo '  -- analyzing centrality class' $centralityCutString

	centralityCut=(`echo $centralityCutString | sed 's/-/ /g' | sed 's/%//g'`)
	lowerLimit=${centralityCut[0]}
	upperLimit=${centralityCut[1]}

	collisionSystemStem=$projectile$target"_"`echo $beamEnergy`"GeV_Nev"$Nevents
	collisionSystemCentralityStem=$projectile$target"_"`echo $beamEnergy`"GeV_C"$lowerLimit"_"$upperLimit"_Nev"$Nevents

	# Clean out previous results
	RESULTS_DIRECTORY=$HOME_DIRECTORY/results
	if [ -d "$RESULTS_DIRECTORY" ]; then
		rm -rf $RESULTS_DIRECTORY
	fi

	mkdir $RESULTS_DIRECTORY

	#=====================================
	# Run Pythia (if desired)
	#PythiaSuccess=$(
	(

		cd $PYTHIA_DIRECTORY
		echo '     Now in '`pwd`

		if $runPythia && [ "$nCC" -eq "0" ]
		then

			# make sure results directory exists
			if [ ! -d "$PYTHIA_RESULTS_DIRECTORY" ]
			then
					mkdir $PYTHIA_RESULTS_DIRECTORY
			fi

			# time and run
			#./run_mainHIC.sh $projectile $target $beamEnergy \
			#					$Nevents $PYTHIA_RESULTS_DIRECTORY
			#./run_testBEeffects.sh $projectile $target $beamEnergy \
			#					$Nevents $PYTHIA_RESULTS_DIRECTORY \
			#					$QRefValue
			./run_BEeffects.sh $projectile $target $beamEnergy \
								$Nevents $PYTHIA_RESULTS_DIRECTORY \
								$QRefValue $lowerLimit $upperLimit

			# check and report whether run was successful
			runSuccess=`echo $?`
			check_success 'Pythia' $runSuccess

			# copy results
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
		#echo $projectile $target $beamEnergy 0 100 $Nevents > $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat
		readlink -f $PYTHIA_RESULTS_DIRECTORY/`cat $recordOfOutputFilename_mult` >> $HBT_EVENT_GEN_DIRECTORY/ensemble_catalogue.dat

		#exit $runSuccess
	)

	#=====================================
	# Run HBT_event_generator
	#HBTegSuccess=$(
	(

		cd $HBT_EVENT_GEN_DIRECTORY
		echo '     Now in '`pwd`

		# using OpenMP (leave a couple cores free)
		export OMP_NUM_THREADS=12

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
		# N.B. - centralities now determined in Pythia

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'HBT_event_generator' $runSuccess

		# copy results
		cp HBT_event_generator.[oe]* \
			./results/* $RESULTS_DIRECTORY/

		readlink -f ./results/HBT_pipiCF.dat > $HBT_FITCF_DIRECTORY/catalogue.dat

		#exit $runSuccess
	)

	#=====================================
	# Run fit_correlation_function
	#fitCFSuccess=$(
	(

		cd $HBT_FITCF_DIRECTORY
		echo '     Now in '`pwd`

		if [ ! -d "./results" ]
		then
			mkdir results
		fi

		cp ../parameters.dat .

		# time and run
		nohup time ./run_fit_correlation_function.e \
				1> fit_correlation_function.out \
				2> fit_correlation_function.err

		# check and report whether run was successful
		runSuccess=`echo $?`
		check_success 'fit_correlation_function' $runSuccess

		# copy results
		cp fit_correlation_function.[oe]* ./results/* $RESULTS_DIRECTORY/

		#exit $runSuccess
	)

	# only store results and clean up if run was successful
	#success=$[PythiaSuccess + HBTegSuccess + fitCFSuccess]

	# check if all codes were successful
	#if [ "$success" -eq "0" ]
	#then
		#add a few more files
		cp ./parameters.dat $RESULTS_DIRECTORY

		zipFilename=$HOME_DIRECTORY/`echo $collisionSystemCentralityStem`"_results.zip"

		zip -r `get_filename $zipFilename` $RESULTS_DIRECTORY

		# Clean-up HBT directories (but not Pythia results directory!!!)
		#rm -rf $HBT_EVENT_GEN_DIRECTORY/*HBT_event_generator.[oe]* $HBT_EVENT_GEN_DIRECTORY/results\
		#	   $HBT_FITCF_DIRECTORY/*fit_correlation_function.[oe]* $HBT_FITCF_DIRECTORY/results
	#fi

	# move on to next centrality class
	nCC=$[nCC+1]

done	# all centralities finished

echo 'Finished everything!'


# End of file
