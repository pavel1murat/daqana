#!/usr/bin/bash
# ------------------------------------------------------------------------------
# example: make_digi_ntuples.sh vst n002 107236 1
#------------------------------------------------------------------------------
    idsid=$1
ntuple_id=$2   # ie 'n002'
       rn=$3
   nfiles=1000 ; if [ ".$4" != "." ] ; then nfiles=$4 ; fi

echo idsid:$idsid rn:$rn nfiles:$nfiles  ntuple_id:$ntuple_id

input_file_list=/tmp/input_00$rn.txt.$$
        logfile=${idsid}.make_${ntuple_id}.$rn.log ;
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

input_fcl=$SPACK_ENV/daqana/fcl/make_${ntuple_id}.fcl

cmd="mu2e -c $input_fcl -S $input_file_list >> $logfile 2>&1 &"
echo "cmd:$cmd";

mu2e -c $input_fcl -S $input_file_list >> $logfile 2>&1 &
