#!/usr/bin/env python
# example :
#        v001/daqana/scripts/make_digi_ntuple.py --rn=107976 --idsid=vst --ntid=n002 --nsbl=10 --calib=01
#------------------------------------------------------------------------------
import subprocess, shutil, datetime
import sys, string, getopt, glob, os, time, re, array
import json
import inspect

class SubmitJob:
    
    def __init__(self):
        self.calib      = None;
        self.idsid      = None;
        self.nfiles     = 1;
        self.ntid       = 'n001';
        self.nsbl       = None;
        self.run_number = None
        self.diag_level = 0;
        
# ---------------------------------------------------------------------
    def Print(self,Name,level,Message):
        if (level > self.diag_level): return 0;
        now = time.strftime('%Y/%m/%d %H:%M:%S',time.localtime(time.time()))
        message = now+' [ SubmitJob::'+Name+' ] '+Message
        print(message)

#------------------------------------------------------------------------------
    def parse_parameters(self):
        name = 'parse_parameters'
        
        self.Print(name,2,'Starting')
        self.Print(name,2, '%s' % sys.argv)

        try:
            optlist, args = getopt.getopt(sys.argv[1:], '',
                     ['calib=', 'diag_level=', 'idsid=', 'ntid=', 'rn=', 'nfiles=', 'nsbl=' ] )
 
        except getopt.GetoptError:
            self.Print(name,0,'%s' % sys.argv)
            self.Print(name,0,'Errors arguments did not parse')
            return 110

        for key, val in optlist:

            # print('key,val = ',key,val)

            if   (key == '--calib'):
                self.calib = val
            elif (key == '--diag_level'):
                self.diag_level = int(val)
            elif (key == '--rn'):
                self.run_number = int(val)
            elif (key == '--idsid'):
                self.idsid = val
            elif (key == '--nfiles'):
                self.nfiles = int(val)
            elif (key == '--ntid'):
                self.ntid = val
            elif (key == '--nsbl'):
                self.nsbl = int(val)

        self.Print(name,1,f'self.diag_level: {self.diag_level}')
        self.Print(name,0,f'self.idsid     : {self.idsid}'     )
        self.Print(name,0,f'self.calib     : {self.calib}'     )
        self.Print(name,0,f'self.rn        : {self.run_number}')
        self.Print(name,0,f'self.ntid      : {self.ntid}'      )
        self.Print(name,0,f'self.nsbl      : {self.nsbl}'      )

        print(f'idsid:{self.idsid}')
#        if (self.fProject == None) :
#            self.Print(name,0,'Error: Project not defined - exiting !')
#            sys.exit(1)

        self.Print(name,1,'------------------------------------- Done')
        return 0

#------------------------------------------------------------------------------
# print statistics reported by a given artdaq process
#------------------------------------------------------------------------------
    def run(self):
        name      = 'run';
#------------------------------------------------------------------------------
# form the input fcl
#------------------------------------------------------------------------------
        input_fcl    = f'/tmp/make_digi_ntuple_{os.getpid()}.fcl';
        template_fcl = os.environ.get('SPACK_ENV')+'/daqana/fcl/make_'+self.ntid+'.fcl';
        print(f'000:template_fcl:{template_fcl}');

        cmd = f'cat {template_fcl} >| {input_fcl}';
#------------------------------------------------------------------------------
# overrides, calib: 'v1'
#------------------------------------------------------------------------------
        if (self.calib):
            cmd += f'; cat {template_fcl} | sed s/calibration_set_v0/calibration_set_v1/ > {input_fcl}';

        if (self.nsbl):
            cmd += f'; echo "physics.analyzers.MakeDigiNtuple.nSamplesBL : {self.nsbl}" >> {input_fcl}';

        print('001:cmd:',cmd);
        
        p   = subprocess.Popen(cmd,
                               executable="/bin/bash",
                               shell=True,
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               encoding="utf-8")
        p.communicate();
#------------------------------------------------------------------------------
# form the input file list
#------------------------------------------------------------------------------
        input_file_list=f'/tmp/make_digi_ntuples_input.{self.run_number}.txt.{os.getpid()}'
        cmd  = "ls -al $RAW_DATA_DIR/* | awk '{print $9}'"
        cmd += f' | grep {self.run_number} | sort';
        if (self.nfiles):
            cmd += f' | head -n {self.nfiles}'

        cmd += f' >| {input_file_list}'

        print('001:cmd:',cmd);

        p = subprocess.Popen(cmd,
                             executable="/bin/bash",
                             shell=True,
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             encoding="utf-8")
        (out, err) = p.communicate();

        self.Print(name,0,f'input_file_list:{input_file_list}')
#------------------------------------------------------------------------------
# submit the job
#-------v----------------------------------------------------------------------
        logfile=f'{self.idsid}.make_{self.ntid}.{self.run_number}.log' ;
        os.system(f'cat {input_fcl} >| {logfile}');
        os.system(f'echo ---------------------------  >> {logfile}');

        cmd = f'mu2e -c {input_fcl} -S {input_file_list} >> {logfile} 2>&1 &'
        p   = subprocess.Popen(cmd,
                             executable="/bin/bash",
                             shell=True,
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             encoding="utf-8")
        (out, err) = p.communicate();
        self.Print(name,1,f'out_2:{out}')

        return;
    
#------------------------------------------------------------------------------
if __name__ == "__main__":

    x = SubmitJob();
    x.parse_parameters();

    x.run()

