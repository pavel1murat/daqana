#!/usr/bin/bash
# setup: 
# ------------------------------------------------------------------------------
    rn=$1
nfiles=1000 ; if [ ".$2" != "." ] ; then nfiles=$2 ; fi

echo rn:$rn nfiles:$nfiles

input_file_list=/tmp/input_00$rn.txt.$$
        logfile=trk_fragment_ana.$rn.log ;
#------------------------------------------------------------------------------
# 1. make list of input files
#------------------------------------------------------------------------------
export DATA_DIR=/data/tracker/vst/mu2etrk_daquser_001/data
ls -al $DATA_DIR/* | grep $rn | awk '{print $9}' | sort | head -n $nfiles > $input_file_list
#------------------------------------------------------------------------------
# 2. process input files
#------------------------------------------------------------------------------
echo "input file_list: " >| $logfile
cat  $input_file_list    >> $logfile

# input_fcl=$SPACK_ENV/daqana/fcl/trk_fragment_ana_`printf "%06i" $rn`.fcl
input_fcl=$SPACK_ENV/daqana/fcl/trk_fragment_ana.fcl

mu2e -c $input_fcl -S $input_file_list >> $logfile 2>&1 &
