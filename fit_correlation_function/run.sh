#! /usr/bin/env bash

# using OpenMP (leave a couple cores free)
export OMP_NUM_THREADS=1

# make sure results directory exists
DIRECTORY=results
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

# time and run
nohup time ./run.e 1> run.out 2> run.err &
