#!/usr/bin/env python
# example :
#        v001/daqana/scripts/submit_mu2e_job.py --c=a.fcl --rn=105935 --idsid=vst --calib=v9 --diag_level=10
#------------------------------------------------------------------------------
import subprocess, shutil, datetime, socket
import sys, string, getopt, glob, os, time, re, array
import json
import inspect

class SubmitJob:
    
    def __init__(self):
        self.calib       = None;
        self.idsid       = None;
        self.nev         = None;
        self.first_event = None;
        self.nfiles      = 1;
        self.fcl         = None;
        self.run_number  = None
        self.diag_level  = 0;
        self.source      = None;
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
                     ['calib=', 'diag_level=', 'c=', 'e=', 'idsid=', 'nev=', 'nfiles=', 'rn=', 's=' ])
                     
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
            elif (key == '--c'):
                self.fcl = val
            elif (key == '--e'):
                self.first_event = val   # '107995:1:10445'
            elif (key == '--idsid'):
                self.idsid = val
            elif (key == '--rn'):
                self.run_number = int(val)
            elif (key == '--nev'):
                self.nev = int(val)
            elif (key == '--nfiles'):
                self.nfiles = int(val)
            elif (key == '--s'):
                self.source = val

        self.Print(name,0,'self.diag_level = %s' % self.diag_level)
        self.Print(name,0,'self.calib      = %s' % self.calib     )
        self.Print(name,0,'self.rn         = %s' % self.run_number)
        self.Print(name,0,'self.fcl        = %s' % self.fcl       )
        self.Print(name,0,'self.nfiles     = %s' % self.nfiles    )

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
# make output directory and cd to there
#------------------------------------------------------------------------------
        fcl_job_stub  = os.path.splitext(os.path.basename(self.fcl))[0];
        now           = datetime.datetime.now()
        formatted_now = now.strftime("%Y-%m-%d-%H-%M")
        host          = socket.gethostname()
        output_dir    = f'results/{formatted_now}.{fcl_job_stub}.{host}.{os.getpid()}'

        cmd  = f'if [ ! -d {output_dir} ] ; then mkdir -p {output_dir} ; fi ;'
        
        print('0001:cmd:',cmd);
        
        p   = subprocess.Popen(cmd,
                               executable="/bin/bash",
                               shell=True,
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               encoding="utf-8")
        p.communicate();
#------------------------------------------------------------------------------
# form the input fcl
#------------------------------------------------------------------------------
        template_fcl = os.getcwd()+'/'+self.fcl;
        pid          = os.getpid();
        job_fcl      = f'{fcl_job_stub}.{pid}.fcl';
        print(f'000:template_fcl:{template_fcl} job_fcl:{job_fcl}');
#------------------------------------------------------------------------------
# overrides, calib: 'v1'
#------------------------------------------------------------------------------
        overrides_cmd = ''
        if (self.calib):
            overrides_cmd  = f' | sed s/calibration_set_v0/calibration_set_v{self.calib}/'
            overrides_cmd += ' | sed s/s\{...\}r\{..\}\{.\}/s\{1\}r\{2\}'+f'{self.calib}/'

#------------------------------------------------------------------------------
# redefinitions --> appends
#------------------------------------------------------------------------------
        os.system(f'cat {template_fcl} {overrides_cmd}                                             >  {output_dir}/{job_fcl}')
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')
        os.system(f'echo "#  overrides by submit_mu2e_job.py"                                      >> {output_dir}/{job_fcl}')  
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')

        x = f'outputs.defaultOutput.fileName: \\"rec.mu2e.trk.vst00s000r01{self.calib}n000.%06r_%06s.art\\"'
        print(f'0011:x:{x}')
        os.system(f'echo {x}                                                                       >> {output_dir}/{job_fcl}')
#------------------------------------------------------------------------------
# form the input file list
#------------------------------------------------------------------------------
        input_file_list=f'/tmp/submit_mu2e_job_input.{self.run_number}.txt.{os.getpid()}'
        cmd  = "ls -al $RAW_DATA_DIR/* | awk '{print $9}'"   ## list all raw files
        cmd += f' | grep {self.run_number} | sort';          ## grep the run number
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
# form the command to execute
#-------v----------------------------------------------------------------------
        input_dsid = self.idsid
        if (input_dsid == None):
            if (self.source):
                # input file defined, assume Mu2e naming conventions
                input_dsid = os.path.basename(self.source).split('.')[3]
                if (input_dsid == 'vst'): input_dsid='vst00s000r000n000'

        run_number = self.run_number
        if (run_number == None):
            if (self.source):
                # input file defined, assume Mu2e naming conventions
                run_number = os.path.basename(self.source).split('.')[4]

        fn = os.getenv("WORK_DIR")+'/.source_me';
        if (os.path.exists(fn)):
            logfile=f'log.mu2e.{input_dsid}.{fcl_job_stub}.{run_number}.log' ;

            cmd  = f'cd $WORK_DIR; source $WORK_DIR/.source_me ;'
            cmd += f' cd {output_dir}; mu2e -c {job_fcl}' # fcl file is in the output_dif
        
            if (self.source     ): cmd += f' -s {self.source}'
            if (self.first_event): cmd += f' -e {self.first_event}'
            if (self.nev        ): cmd += f' -n {self.nev}'

            cmd += f' >> {logfile} 2>&1 &'
#------------------------------------------------------------------------------
# print the command and log it, together with the FCL file, in the log file
#-----------v------------------------------------------------------------------
            self.Print(name,0,cmd)

            os.system(f'echo "cmd:{cmd}"                               >| {output_dir}/{logfile}')
            os.system(f'echo "---------------------------------------" >> {output_dir}/{logfile}')
            os.system(f'cat {output_dir}/{job_fcl}                     >> {output_dir}/{logfile}')
            os.system(f'echo "---------------------------------------" >> {output_dir}/{logfile}')
               
            p   = subprocess.Popen(cmd,
                                   executable="/bin/bash",
                                   shell=True,
                                   stderr=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   encoding="utf-8")
            (out, err) = p.communicate();
            self.Print(name,1,f'out_2:{out}')
            self.Print(name,0,f'err_2:{err}')
        else:
            self.Print(name,0,f'ERROR: {os.getenv("WORK_DIR")}/.source_me doesn\'t exist, BAIL OUT')
            
        return;
#------------------------------------------------------------------------------
if __name__ == "__main__":

    x = SubmitJob();
    x.parse_parameters();

    x.run()

