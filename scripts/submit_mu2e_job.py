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
import inspect, logging, signal

#------------------------------------------------------------------------------
# configure logger
#------------------------------------------------------------------------------
logger = logging.getLogger('submit_mu2e_job')

logger.setLevel(logging.DEBUG)

# File handler
fh = logging.FileHandler('submit_mu2e_job.log', mode='a')
fh.setLevel(logging.INFO)

# Console handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# Formatter with timestamps
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s [%(filename)s:%(lineno)d:%(funcName)s]: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
# fh.setFormatter(formatter)
ch.setFormatter(formatter)

# logger.addHandler(fh)
logger.addHandler(ch)

class SubmitJob:
    
    def __init__(self):
        self.rc = None;

# ---------------------------------------------------------------------
#    def Print(self,Name,level,Message):
#        if (level > self.args.diag_level): return 0;
#        now = time.strftime('%Y/%m/%d %H:%M:%S',time.localtime(time.time()))
#        message = now+' [ SubmitJob::'+Name+' ] '+Message
#        print(message)

#------------------------------------------------------------------------------
    def parse_parameters(self):
        name = 'parse_parameters'
        
        logger.info('Starting')
        logger.info(f'sys.argv:{sys.argv}')
        
        parser = argparse.ArgumentParser()

        parser.add_argument("--calib-ver"       , default=None,           help="calibration version, defaults to 0")
        parser.add_argument("--calib-run"       , default=None,           help="use calibrations keyed on a given run")
        parser.add_argument("--diag_level"      , type=int, default=0,    help="Path to the configuration file")
        parser.add_argument("--dry-run"         , action='store_true',    help="dry run, if specified")
        parser.add_argument('-c',"--fcl"        , default=None,           help="Path to the configuration file")
        parser.add_argument('-e',"--first_event", type=int, default=None, help="Path to the configuration file")
        parser.add_argument('--idsid'           , default=None,           help="input dataset ID")
        parser.add_argument('-n','--nevents'    , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('--nfiles'          , type=int, default=None, help="N(files) to process")
        parser.add_argument('--nskip'           , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('-r','--run_number' , type=int, default=None, help="Path to the configuration file")
        parser.add_argument('-s','--source'     , default=None,           help="input file, as in art")
        parser.add_argument('-S','--Source'     , default=None,           help="input file list, as in art")

        args = parser.parse_args()

        logger.info(f'self.diag_level = {args.diag_level}'   )
        logger.info(f'self.calib_ver  = {args.calib_ver}'    )
        logger.info(f'self.calib_run  = {args.calib_run}'    )
        logger.info(f'self.rn         = {args.run_number}'   )
        logger.info(f'self.fcl        = {args.fcl}'          )
        logger.info(f'self.nfiles     = {args.nfiles}'       )

#        if (self.fProject == None) :
#            self.Print(name,0,'Error: Project not defined - exiting !')
#            sys.exit(1)

        logger.info(f'------------------------------------- Done')
        return args
   
#------------------------------------------------------------------------------
# print statistics reported by a given artdaq process
#------------------------------------------------------------------------------
    def run(self,args):
        name      = 'run';
#------------------------------------------------------------------------------
# make output directory and cd to there
#------------------------------------------------------------------------------
        fcl_job_stub  = os.path.splitext(os.path.basename(args.fcl))[0];
        now           = datetime.datetime.now()
        formatted_now = now.strftime("%Y-%m-%d-%H-%M")
        host          = socket.gethostname()
        output_dir    = f'results/{formatted_now}.{fcl_job_stub}.{host}.{os.getpid()}'

        cmd  = f'if [ ! -d {output_dir} ] ; then mkdir -p {output_dir} ; fi ;'
        
        logger.debug(f'0001:cmd:{cmd}');
        
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
        template_fcl = os.getcwd()+'/'+args.fcl;
        pid          = os.getpid();
        job_fcl      = f'{fcl_job_stub}.{pid}.fcl';
        logger.debug(f'000:template_fcl:{template_fcl} job_fcl:{job_fcl}');
#------------------------------------------------------------------------------
# overrides, calib: 'v1'
#------------------------------------------------------------------------------
        overrides_cmd = ''
        if (args.calib_ver):
            overrides_cmd  = f' | sed s/calibration_set_v0/calibration_set_v{args.calib_ver}/'
            # overrides_cmd += ' | sed s/s\{...\}r\{..\}\{.\}/s\{1\}r\{2\}'+f'{args.calib_set}/'

        if (args.calib_run):
            subdir = args.calib_run[0:3]+'000'
            if (overrides_cmd != ''):
                overrides_cmd  += f' | sed s!fcl/calibration_set!rundb/{subdir}/{args.calib_run}/calibration_set_{args.calib_run}!'
            else:
                overrides_cmd  = f' | sed s!fcl/calibration_set!rundb/{subdir}/{args.calib_run}/calibration_set_{args.calib_run}!'

        print(f'cmd:{overrides_cmd}')
#------------------------------------------------------------------------------
# redefinitions --> appends 
#------------------------------------------------------------------------------
        os.system(f'cat {template_fcl} {overrides_cmd}                                             >  {output_dir}/{job_fcl}')
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')
        os.system(f'echo "#  overrides by submit_mu2e_job.py"                                      >> {output_dir}/{job_fcl}')  
        os.system(f'echo "#----------------------------------------------------------------------" >> {output_dir}/{job_fcl}')

        if (args.dry_run):
            return
#        x = f'outputs.defaultOutput.fileName: \\"rec.mu2e.trk.vst00s000r01{args.calib_set}n000.%06r_%06s.art\\"'
#        print(f'0011:x:{x}')
#        os.system(f'echo {x}                                                                       >> {output_dir}/{job_fcl}')
#------------------------------------------------------------------------------
# form the input file list
#-------v----------------------------------------------------------------------
        input_file_list=None
        if (not args.source):
            input_file_list=f'/tmp/submit_mu2e_job_input.{args.run_number}.txt.{os.getpid()}'
            cmd  = "ls -al $RAW_DATA_DIR/* | awk '{print $9}'"   ## list all raw files
            cmd += f' | grep {args.run_number} | sort';          ## grep the run number
            if (args.nfiles):
                cmd += f' | head -n {args.nfiles}'

            cmd += f' >| {input_file_list}'

            logger.debug(f'001:cmd:{cmd}');

            p = subprocess.Popen(cmd,
                                 executable="/bin/bash",
                                 shell=True,
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 encoding="utf-8")
            (out, err) = p.communicate();

            logger.info(f'input_file_list:{input_file_list}')
#------------------------------------------------------------------------------
# form the command to execute
#-------v----------------------------------------------------------------------
        input_dsid = args.idsid
        if (input_dsid == None):
            if (args.source):
                # input file defined, assume Mu2e naming conventions
                input_dsid = os.path.basename(args.source).split('.')[3]
                if (input_dsid == 'vst'): input_dsid='vst00s000r000n000'

        run_number = args.run_number
        if (run_number == None):
            if (args.source):
                # input file defined, assume Mu2e naming conventions
                run_number = os.path.basename(args.source).split('.')[4]

        fn = os.getenv('WORK_DIR')+'/.source_me';
        if (os.path.exists(fn)):
            logfile=f'log.mu2e.{input_dsid}.{fcl_job_stub}.{run_number}.log' ;

            cmd  = f'cd $WORK_DIR; source $WORK_DIR/.source_me ;'
            cmd += f' cd {output_dir}; mu2e -c {job_fcl}' # fcl file is in the output_dif
        
            if (args.Source     ): cmd += f' -S {args.Source}'
            else:
                if (args.source     ): cmd += f' -s {args.source}'
                else                 : cmd += f' -S {input_file_list}'
            
            if (args.first_event): cmd += f' -e {args.first_event}'
            if (args.nevents    ): cmd += f' -n {args.nevents}'

            cmd += f' >> {logfile} 2>&1 &'
#------------------------------------------------------------------------------
# print the command and log it, together with the FCL file, in the log file
#-----------v------------------------------------------------------------------
            logger.info(f'cmd:{cmd}')

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
            logger.info(f'out_2:{out}')
            logger.info(f'err_2:{err}')
        else:
            logger.warning(f'{os.getenv("WORK_DIR")}/.source_me doesn\'t exist, BAIL OUT')
            
        return;
#------------------------------------------------------------------------------
if __name__ == "__main__":

    x    = SubmitJob();
    args = x.parse_parameters();

    x.run(args)

