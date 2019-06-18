#! /usr/bin/env bash

rm -rf results
mkdir results

gmake distclean && gmake all && cp ../parameters.dat . && OMP_NUM_THREADS=10 nohup ./run_HBT_event_generator.e file_mode=0 &> run.out &
