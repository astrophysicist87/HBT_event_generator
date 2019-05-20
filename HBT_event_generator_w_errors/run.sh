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

# time and run
nohup time ./run_HBT_event_generator.e file_mode=0 \
		1> HBT_event_generator.out \
		2> HBT_event_generator.err &
