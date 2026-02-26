#!/usr/bin/bash
#------------------------------------------------------------------------------
# call example: ./rdm_copy.sh 120340 120341
#------------------------------------------------------------------------------
RBASE=/pnfs/mu2e/persistent/users/mu2epro/RDM/remote/trk
# cp $FN $RBASE/output/$FN ; mv $RBASE/output/$FN $RBASE/upload/$FN

rn1=`echo $1 | awk -F : '{print $1}'`
rn2=`echo $1 | awk -F : '{print $2}'`

if [ ".$rn2" == "." ] ; then rn2=$rn1 ; fi

doit=$2

echo rn1=$rn1 rn2=$rn2 doit=$doit
# return

# name_stub=${USER}_`echo $MU2E_DAQ_DIR | awk -F / '{print $NF}'`

for rn in `seq $rn1 $rn2` ; do
    irn=`printf "%06i" $rn`
    for f in `ls /exp/mu2e/data/projects/vst/datasets/raw.mu2e.trk.vst.art/raw.mu2e.trk.vst.*.art | grep $rn` ; do
        bn=`basename $f`
        dsconf=`echo $bn | awk -F . '{print $4}'`
        dsid=raw.mu2e.trk.vst.art
        # cmd="scp $f murat@mu2egpvm06:/exp/mu2e/data/projects/vst/datasets/$dsid/."
        cmd="cp $f $RBASE/output/$bn" 
        echo "$cmd"
        if [ ".$doit" != "." ] ; then 
            # echo doit=$doit
            $cmd ; echo rc=$? ;
            mv $RBASE/output/$bn $RBASE/upload/$bn
        fi
    done
done
