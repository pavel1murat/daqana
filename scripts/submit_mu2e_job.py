#!/usr/bin/env python
# example :
# ---------
#  v001/daqana/scripts/submit_mu2e_job.py --c=a.fcl --rn=105935 --idsid=vst --calib=v0 --diag_level=10
#
# --rn  : use files of a given run number
# --fcl : 
#------------------------------------------------------------------------------
import subprocess, shutil, datetime, socket
import sys, string, argparse, glob, os, time, re, array
import json
import inspect

class SubmitJob:
    
    def __init__(self):
        self.args = None;

# ---------------------------------------------------------------------
    def Print(self,Name,level,Message):
        if (level > self.args.diag_level): return 0;
        now = time.strftime('%Y/%m/%d %H:%M:%S',time.localtime(time.time()))
        message = now+' [ SubmitJob::'+Name+' ] '+Message
        print(message)

#------------------------------------------------------------------------------
    def parse_parameters(self):
        name = 'parse_parameters'
        
#        self.Print(name,2,'Starting')
#        self.Print(name,2, '%s' % sys.argv)
        
        parser = argparse.ArgumentParser()

        parser.add_argument("--calib_version"   , default=None,           help="Path to the configuration file")
        parser.add_argument("--diag_level"      , type=int, default=0,    help="Path to the configuration file")
        parser.add_argument('-c',"--fcl"        , default=None,           help="Path to the configuration file")
        parser.add_argument('-e',"--first_event", type=int, default=None, help="Path to the configuration file")
        parser.add_argument('--idsid'           , default=None,           help="input dataset ID")
        parser.add_argument('-n','--nevents'    , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('--nfiles'          , type=int, default=None, help="N(files) to process")
        parser.add_argument('--nskip'           , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('-r','--run_number' , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('-s','--source'     , default=None,           help="input file, as in art")
        parser.add_argument('-S','--Source'     , default=None,           help="input file list, as in art")

        self.args = parser.parse_args()

        self.Print(name,0,'self.diag_level = %s' % self.args.diag_level)
        self.Print(name,0,'self.calib      = %s' % self.args.calib_version)
        self.Print(name,0,'self.rn         = %s' % self.args.run_number)
        self.Print(name,0,'self.fcl        = %s' % self.args.fcl       )
        self.Print(name,0,'self.nfiles     = %s' % self.args.nfiles    )

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
        fcl_job_stub  = os.path.splitext(os.path.basename(self.args.fcl))[0];
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
# form input fcl
#------------------------------------------------------------------------------
        template_fcl = os.getcwd()+'/'+self.args.fcl;
        pid          = os.getpid();
        job_fcl      = f'{fcl_job_stub}.{pid}.fcl';
        print(f'000:template_fcl:{template_fcl} job_fcl:{job_fcl}');
#------------------------------------------------------------------------------
# overrides, calib: 'v1'
#------------------------------------------------------------------------------
        overrides_cmd = ''
        if (self.args.calib_version):
            overrides_cmd  = f' | sed s/calibration_set_v0/calibration_set_v{self.args.calib}/'
            overrides_cmd += ' | sed s/s\{...\}r\{..\}\{.\}/s\{1\}r\{2\}'+f'{self.args.calib}/'

#------------------------------------------------------------------------------
# redefinitions --> appends 
#------------------------------------------------------------------------------
        os.system(f'cat {template_fcl} {overrides_cmd}                                             >  {output_dir}/{job_fcl}')
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')
        os.system(f'echo "#  overrides by submit_mu2e_job.py"                                      >> {output_dir}/{job_fcl}')  
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')

        x = f'outputs.defaultOutput.fileName: \\"rec.mu2e.trk.vst00s000r01{self.args.calib_version}n000.%06r_%06s.art\\"'
        print(f'0011:x:{x}')
        os.system(f'echo {x}                                                                       >> {output_dir}/{job_fcl}')
#------------------------------------------------------------------------------
# form the input file list
#-------v----------------------------------------------------------------------
        input_file_list=None
        if (not self.args.source):
            input_file_list=f'/tmp/submit_mu2e_job_input.{self.args.run_number}.txt.{os.getpid()}'
            cmd  = "ls -al $RAW_DATA_DIR/* | awk '{print $9}'"   ## list all raw files
            cmd += f' | grep {self.args.run_number} | sort';          ## grep the run number
            if (self.args.nfiles):
                cmd += f' | head -n {self.args.nfiles}'

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
        input_dsid = self.args.idsid
        if (input_dsid == None):
            if (self.args.source):
                # input file defined, assume Mu2e naming conventions
                input_dsid = os.path.basename(self.args.source).split('.')[3]
                if (input_dsid == 'vst'): input_dsid='vst00s000r000n000'

        run_number = self.args.run_number
        if (run_number == None):
            if (self.args.source):
                # input file defined, assume Mu2e naming conventions
                run_number = os.path.basename(self.args.source).split('.')[4]

        fn = os.getenv("WORK_DIR")+'/.source_me';
        if (os.path.exists(fn)):
            logfile=f'log.mu2e.{input_dsid}.{fcl_job_stub}.{run_number}.log' ;

            cmd  = f'cd $WORK_DIR; source $WORK_DIR/.source_me ;'
            cmd += f' cd {output_dir}; mu2e -c {job_fcl}' # fcl file is in the output_dif
        
            if (self.args.source     ): cmd += f' -s {self.args.source}'
            else:                       cmd += f' -S {input_file_list}'
            
            if (self.args.first_event): cmd += f' -e {self.args.first_event}'
            if (self.args.nevents    ): cmd += f' -n {self.args.nevents}'

            cmd += f' >> {logfile} 2>&1 &'
#------------------------------------------------------------------------------
# print the command and log it, together with the FCL file, in the log file
#-----------v------------------------------------------------------------------
            self.Print(name,0,cmd)

            os.system(f'echo "cmd:{cmd}"                               >| {output_dir}/{logfile}')
            os.system(f'echo "---------------------------------------" >> {output_dir}/{logfile}')
            os.system(f'cat {output_dir}/{job_fcl}                     >> {output_dir}/{logfile}')
            os.system(f'echo "---------------------------------------" >> {output_dir}/{logfile}')
            if (input_file_list):
                os.system(f'echo "input file list:"                        >> {output_dir}/{logfile}')
                os.system(f'cat  {input_file_list}                         >> {output_dir}/{logfile}')
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
            self.Print(name,0,f'WARNING: {os.getenv("WORK_DIR")}/.source_me doesn\'t exist, BAIL OUT')
            
        return;
#------------------------------------------------------------------------------
if __name__ == "__main__":

    x = SubmitJob();
    x.parse_parameters();

    x.run()

