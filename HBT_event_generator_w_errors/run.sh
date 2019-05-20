#! /usr/bin/env bash

# using OpenMP (leave a couple cores free)
#export OMP_NUM_THREADS=`nproc --all`
export OMP_NUM_THREADS=12

# make sure results directory exists
DIRECTORY=results
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

#lowerLimit=$1
#upperLimit=$2

cp ../parameters.dat .
cp ../parameters.dat ./results/


(
for mult in 4 5 6 8 10 100 1000
do
	for nLoops in 10 100 1000 10000 100000 1000000
	do
		for bw in 0.5 0.1 0.05 0.25 0.1 0.05 0.025 0.01 0.005
		do

			# time and run
			nohup time ./run_HBT_event_generator.e
					file_mode=0 RNG_mult=$mult \
					RNG_nLoops=$nLoops
					1> HBT_event_generator.out \
					2> HBT_event_generator.err

		done
	done
done
) &


# End of file
