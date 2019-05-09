#! /usr/bin/env bash

#=====================================
# Header info
HOME_DIRECTORY=~/HBT_event_generator
# Pythia
PYTHIA_DIRECTORY=~/pythia8235_MYCOPY/examples
# HBT event generator
HBT_EVENT_GEN_DIRECTORY=$HOME_DIRECTORY/HBT_event_generator_w_errors
# Fit correlation function
HBT_FITCF_DIRECTORY=$HOME_DIRECTORY/fit_correlation_function
# Source variances/HBT radii
HBT_SV_DIRECTORY=$HOME_DIRECTORY/source_variances

success=0

#=====================================
# Compile Pythia
echo '#====================================='
echo '# Compiling Pythia'
echo '#====================================='
cd $PYTHIA_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
rm main_BEeffects
make main_BEeffects
success=$[success+`echo $?`]

#=====================================
# Compile HBT_event_generator
echo '#====================================='
echo '# Compiling HBT_event_generator'
echo '#====================================='
cd $HBT_EVENT_GEN_DIRECTORY
echo 'In directory='`pwd`':'
echo '#====================================='
export OMP_NUM_THREADS=10
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
# Compile SV.e
echo '#====================================='
echo '# Compiling SV.e'
echo '#====================================='
cd $HBT_SV_DIRECTORY
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
