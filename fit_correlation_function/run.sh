#! /usr/bin/env bash

# using OpenMP (leave a couple cores free)
export OMP_NUM_THREADS=1

# make sure results directory exists
DIRECTORY=results
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

# time and run
nohup time ./run_fit_correlation_function.e \
				1> fit_correlation_function.out \
				2> fit_correlation_function.err
