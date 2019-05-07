#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY=results
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

cp ../parameters.dat .

# time and run
nohup time ./SV.e \
		1> SV_record.out \
		2> SV_record.err &
