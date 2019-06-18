#! /usr/bin/env bash

# using OpenMP (leave a couple cores free)
#export OMP_NUM_THREADS=`nproc --all`
export OMP_NUM_THREADS=12

rm -rf results auto
mkdir auto

for mult in 4 5 6 10 100 1000
do
	for nLoops in 100 10000 1000000
	do
		for bw in 0.5 0.1 0.05 0.025 0.01 0.005
		do

			mkdir results

			# time and run
			nohup time ./run_HBT_event_generator.e \
					file_mode=0 RNG_mult=$mult \
					RNG_nLoops=$nLoops bin_epsilon=$bw \
					1> HBT_event_generator.out \
					2> HBT_event_generator.err

			# make sure results directory exists
			DIRECTORY=auto/results_mult`echo $mult`_nLoops`echo $nLoops`_bw`echo $bw`
			mv results $DIRECTORY

			cp ../parameters.dat .
			cp ../parameters.dat $DIRECTORY

		done
	done
done


# End of file
