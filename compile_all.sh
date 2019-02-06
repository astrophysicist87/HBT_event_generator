#! /usr/bin/env bash

#=====================================
# Header info
HOME_DIRECTORY=~/HBT_event_generator
# Pythia
PYTHIA_DIRECTORY=~/pythia8235/examples
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HOME_DIRECTORY/fit_correlation_function

success=0

#=====================================
# Compile Pythia
echo '#====================================='
echo '# Compiling Pythia'
echo '#====================================='
cd $PYTHIA_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
rm mainHIC main_testBEeffects
make mainHIC main_testBEeffects
success=$[success+`echo $?`]

#=====================================
# Compile HBT_event_generator
echo '#====================================='
echo '# Compiling HBT_event_generator'
echo '#====================================='
cd $HBT_EVENT_GEN_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
export OMP_NUM_THREADS=1
gmake distclean
gmake all
success=$[success+`echo $?`]

#=====================================
# Compile fit_correlation_function
echo '#====================================='
echo '# Compiling fit_correlation_function'
echo '#====================================='
cd $HBT_FITCF_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
gmake distclean
gmake all
success=$[success+`echo $?`]
echo '#====================================='

#=====================================
# Check success
cd $HOME_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
if [ "$success" -eq "0" ]
then
	echo 'END RESULT: Everything compiled successfully!'
else
	echo 'END RESULT: There were problems compiling.'
	exit $success
fi

# End of file
