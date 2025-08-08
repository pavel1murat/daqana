#!/usr/bin/env python
# v001/daqana/scripts/make_station_hist.py --rn=107976 --idsid=vst --min_edep=0.0005
#------------------------------------------------------------------------------
import subprocess, shutil, datetime
import sys, string, getopt, glob, os, time, re, array
import json
import inspect

class SubmitJob:
    
    def __init__(self):
        self.idsid      = 'vst';
        self.nfiles     = None;
        self.min_edep   = None;
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
                     ['idsid=', 'min_edep=', 'rn=', 'nsbl=' ] )
 
        except getopt.GetoptError:
            self.Print(name,0,'%s' % sys.argv)
            self.Print(name,0,'Errors arguments did not parse')
            return 110

        for key, val in optlist:

            # print('key,val = ',key,val)

            if   (key == '--diag_level'):
                self.diag_level = int(val)
            elif (key == '--rn'):
                self.run_number = int(val)
            elif (key == '--idsid'):
                self.idsid = val
            elif (key == '--nfiles'):
                self.nfiles = int(val)
            elif (key == '--min_edep'):
                self.min_edep = float(val)

        self.Print(name,1,'verbose   = %s' % self.diag_level)
        self.Print(name,0,'rn        = %s' % self.run_number)
        self.Print(name,0,'min_edep  = %s' % self.min_edep  )

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
        input_fcl    = f'/tmp/make_station_hist_{os.getpid()}.fcl';
        template_fcl = os.environ.get('SPACK_ENV')+'/daqana/fcl/make_station_hist.fcl';

        cmd = f'cat {template_fcl} >| {input_fcl}';
        if (self.min_edep):
            cmd += f'; echo "physics.analyzers.StationAna.minEDep : {self.min_edep}" >> {input_fcl}';

        print('000:cmd:',cmd);
        
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
        input_file_list=f'/tmp/make_station_hist_input.{self.run_number}.txt.{os.getpid()}'
        cmd  = "ls -al /data/tracker/vst/mu2etrk_daquser_001/data/* | awk '{print $9}'"
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
        logfile=f'{self.idsid}.make_station_hist.{self.run_number}.log' ;
        
        os.system(f'cat {input_fcl}                   >| {logfile}');
        os.system(f'echo ---------------------------  >> {logfile}');
        os.system(f'cat {input_file_list}             >| {logfile}');
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

