#!/usr/bin/bash
# setup: tstation, roctower, ...
# ------------------------------------------------------------------------------
         rn=$1
input_files=/tmp/input_${setup}_00$rn.txt.$$
    logfile=trkvst.trk_fragment_ana.00$rn.log ;
#------------------------------------------------------------------------------
# 1. make list of input files
#------------------------------------------------------------------------------
ls -al /exp/mu2e/data/projects/vst/datasets/raw.mu2e.trkvst.$setup.art/* | grep $rn | awk '{print $9}' > $input_files
#------------------------------------------------------------------------------
# 2. process input files
#------------------------------------------------------------------------------
echo "input files: " >| $logfile
cat  $input_files    >> $logfile

input_fcl=daqana/fcl/trk_fragment_ana_`printf "%06i" $rn`.fcl

mu2e -c $input_fcl -S $input_files >> $logfile 2>&1 &
