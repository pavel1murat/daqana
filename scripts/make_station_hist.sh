#!/usr/bin/bash
# ------------------------------------------------------------------------------
# example: make_station_hist.sh vst 107236
#------------------------------------------------------------------------------
   idsid=$1
      rn=$2
min_edep=0.006 ; if [ ".$3" != "." ] ; then min_edep=$3; fi
      
echo idsid:$idsid rn:$rn min_edep:$min_edep

input_file_list=/tmp/input_00$rn.txt.$$
        logfile=${idsid}.$rn.make_station_hist.log ;
#------------------------------------------------------------------------------
# 1. make list of input files
#------------------------------------------------------------------------------
export DATA_DIR=/data/tracker/vst/mu2etrk_daquser_001/data
ls -al $DATA_DIR/* | grep $rn | awk '{print $9}' | sort >| $input_file_list
#------------------------------------------------------------------------------
# 2. process input files
#------------------------------------------------------------------------------
# input_fcl=$SPACK_ENV/daqana/fcl/make_station_hist_$thr.fcl
  tmp_fcl=`mktemp /tmp/make_station_hist-XXXX.fcl`
input_fcl=$SPACK_ENV/daqana/fcl/make_station_hist.fcl

cat  $input_fcl                                         >| $tmp_fcl
echo "physics.analyzers.StationAna.minEDep : $min_edep" >> $tmp_fcl

echo "input file_list: " >| $logfile
cat  $input_file_list    >> $logfile
echo ">>> input FCL:"    >> $logfile
cat  $tmp_fcl            >> $logfile

cmd="mu2e -c $tmp_fcl -S $input_file_list >> $logfile 2>&1 &"
echo "cmd:$cmd";

mu2e -c $tmp_fcl -S $input_file_list >> $logfile 2>&1 &
