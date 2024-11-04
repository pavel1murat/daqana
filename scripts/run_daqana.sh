#!/usr/bin/bash
# setup: tstation, roctower, ...
# ------------------------------------------------------------------------------
setup=$1
rn=$2
input_files=/tmp/input_${setup}_00$rn.txt.$$

ls -al /exp/mu2e/data/projects/tracker/vst/datasets/raw.mu2e.trkvst.$setup.art/* | grep $rn | awk '{print $9}' > $input_files
mu2e -c ./test_trk_fragment_ana.fcl -S $input_files >| trkvst.$setup.trk_fragment_ana.00$rn.log 2>&1 &
return $?

