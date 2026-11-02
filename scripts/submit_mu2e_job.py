#!/usr/bin/env python
# example :
# ---------
#  v001/daqana/scripts/submit_mu2e_job.py --fcl=a.fcl --run_number=105935 --idsid=vst --calib-ver=v0 --diag_level=10
#
# --run_number : use files of a given run number
# --fcl        : template fcl file
#------------------------------------------------------------------------------
import subprocess, datetime, socket
import sys, argparse, glob, os, shlex
import logging

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
fh.setFormatter(formatter)
ch.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(ch)

class SubmitJob:

#------------------------------------------------------------------------------
    def parse_parameters(self):
        name = 'parse_parameters'

        logger.info('Starting')
        logger.info(f'sys.argv:{sys.argv}')

        parser = argparse.ArgumentParser()

        parser.add_argument("-b","--book"       , default=""  ,           help="book name")
        parser.add_argument("--calib-ver"       , default=None,           help="calibration version, defaults to v0")
        parser.add_argument("--calib-run"       , default=None,           help="use calibrations keyed on a given run")
        parser.add_argument("-D", "--data-dir"  , default=None,           help="input data directory")
        parser.add_argument("--diag_level"      , type=int, default=0,    help="diagnostic level")
        parser.add_argument("--dry-run"         , action='store_true',    help="print the job setup, create nothing")
        parser.add_argument('-c',"--fcl"        , required=True,          help="template fcl file")
        parser.add_argument('-e',"--first_event", default=None,           help="event ID of the first event to process")
        parser.add_argument('--idsid'           , default=None,           help="input dataset ID")
        parser.add_argument('-n','--nevents'    , type=int, default=None, help="N(events) to process")
        parser.add_argument('--nfiles'          , type=int, default=None, help="N(files) to process")
        parser.add_argument('--nskip'           , type=int, default=None, help="unused")
        parser.add_argument('--odsid'           , default=None,           help="output dataset ID")
        parser.add_argument('-r','--run_number' , type=int, default=None, help="run number")
        parser.add_argument('-s','--source'     , default=None,           help="input file, as in art")
        parser.add_argument('-S','--Source'     , default=None,           help="input file list, as in art")

        args = parser.parse_args()

        logger.info(f'args.idsid      = {args.idsid}'        )
        logger.info(f'args.diag_level = {args.diag_level}'   )
        logger.info(f'args.calib_ver  = {args.calib_ver}'    )
        logger.info(f'args.calib_run  = {args.calib_run}'    )
        logger.info(f'args.rn         = {args.run_number}'   )
        logger.info(f'args.fcl        = {args.fcl}'          )
        logger.info(f'args.nfiles     = {args.nfiles}'       )
        logger.info(f'args.data_dir   = {args.data_dir}'     )

        logger.info(f'------------------------------------- Done')
        return args

#------------------------------------------------------------------------------
# set up the job directory, the job fcl and the log file, then submit the job
#------------------------------------------------------------------------------
    def run(self,args):
        q = shlex.quote
#------------------------------------------------------------------------------
# check the environment before anything is created
#------------------------------------------------------------------------------
        work_dir  = os.getenv('WORK_DIR')
        source_me = os.path.join(work_dir or '', '.source_me')
        if (not work_dir or not os.path.exists(source_me)):
            # a dry run only prints the setup, so it can go ahead without the environment
            if (not args.dry_run):
                logger.error(f'WORK_DIR={work_dir}: .source_me not found, BAIL OUT')
                sys.exit(1)
            logger.warning(f'WORK_DIR={work_dir}: .source_me not found, continuing the dry run')

        if (args.calib_ver == 'latest' and (not args.calib_run or not os.getenv('RUNDB_DIR'))):
            logger.error('--calib-ver=latest needs --calib-run and $RUNDB_DIR')
            sys.exit(1)

        if (args.idsid):
            muse_work_dir = os.getenv('MUSE_WORK_DIR')
            if (not muse_work_dir):
                logger.error('--idsid needs $MUSE_WORK_DIR, samdfc/ is a subdirectory of it')
                sys.exit(1)
            args.Source = os.path.join(muse_work_dir, 'samdfc', f'{args.idsid}.files')
#------------------------------------------------------------------------------
# form the input file list, input_dsid and run number come from it
#------------------------------------------------------------------------------
        input_dsid, run_number, file_list = self.resolve_inputs(args)
#------------------------------------------------------------------------------
# output directory and job names, absolute paths: the job starts inside WORK_DIR
#------------------------------------------------------------------------------
        fcl_job_stub  = os.path.splitext(os.path.basename(args.fcl))[0];
        pid           = os.getpid()
        now           = datetime.datetime.now()
        formatted_now = now.strftime("%Y-%m-%d-%H-%M")
        host          = socket.gethostname()
        output_dir    = os.path.abspath(f'results/{formatted_now}.{fcl_job_stub}.{host}.{pid}')
        template_fcl  = os.path.abspath(args.fcl)
        job_fcl       = os.path.join(output_dir, f'{fcl_job_stub}.{pid}.fcl')
        logfile       = f'log.mu2e.{input_dsid}.{fcl_job_stub}.{run_number}.log'
#------------------------------------------------------------------------------
# calibrations and overrides, see overrides_cmd
#------------------------------------------------------------------------------
        overrides = self.overrides_cmd(args, template_fcl, input_dsid)
        logger.info(f'overrides_cmd:{overrides}')

        if (args.dry_run):
            if (file_list is not None):
                logger.info(f'input file list:\n{file_list}')
            else:
                logger.info(f'input: -s {args.source}')
            return

        os.makedirs(output_dir, exist_ok=True)
#------------------------------------------------------------------------------
# process redefinitions and appends
#------------------------------------------------------------------------------
        self.run_shell(f'{overrides} > {q(job_fcl)}')
        with open(job_fcl, 'a', encoding="utf-8") as f:
            f.write('#----------------------------------------------------------------------\n')
            f.write('#  overrides by submit_mu2e_job.py\n')
            f.write('#----------------------------------------------------------------------\n')
#------------------------------------------------------------------------------
# write the input file list into the output directory
#------------------------------------------------------------------------------
        input_list = None
        if (file_list is not None):
            input_list = os.path.join(output_dir, f'submit_mu2e_job_input.{run_number}.txt')
            with open(input_list, 'w', encoding="utf-8") as f:
                f.write(file_list)
#------------------------------------------------------------------------------
# form the command to execute
#-------v----------------------------------------------------------------------
        cmd  = f'cd {q(work_dir)}; source {q(source_me)};'
        for name in ('.rootrc', 'rootlogon.C'):
            src  = os.path.join(work_dir, name)
            cmd += f' if [ -f {q(src)} ] ; then cp {q(src)} {q(output_dir)}/. ; fi ;'
        cmd += f' cd {q(output_dir)}; mu2e -c {q(job_fcl)}'

        if (input_list):
            cmd += f' -S {q(input_list)}'
        elif (args.source):
            cmd += f' -s {q(args.source)}'

        if (args.first_event): cmd += f' -e {q(args.first_event)}'
        if (args.nevents    ): cmd += f' -n {args.nevents}'

        cmd += f' >> {q(logfile)} 2>&1 &'
#------------------------------------------------------------------------------
# print the command and log it, together with the FCL file, in the log file
#-------v----------------------------------------------------------------------
        logger.info(f'cmd:{cmd}')

        with open(os.path.join(output_dir, logfile), 'w', encoding="utf-8") as log:
            log.write(f'cmd:{cmd}\n')
            log.write('---------------------------------------\n')
            with open(job_fcl, encoding="utf-8") as f:
                log.write(f.read())
            log.write('---------------------------------------\n')
            if (file_list is not None):
                log.write('input file list:\n' + file_list)
                log.write('---------------------------------------\n')
#------------------------------------------------------------------------------
# execute the command, the shell returns once mu2e is backgrounded,
# so this checks the setup steps only
#------------------------------------------------------------------------------
        self.run_shell(cmd)
        logger.info(f'submitted, log in {os.path.join(output_dir, logfile)}')

#------------------------------------------------------------------------------
# returns (input_dsid, run_number, file_list), file_list is None when --source is used
#------------------------------------------------------------------------------
    def resolve_inputs(self,args):
        input_dsid = args.idsid
        run_number = args.run_number

        if (input_dsid is None and args.source):
            # input file defined, assume Mu2e naming conventions
            input_dsid = os.path.basename(args.source).split('.')[2]
        if (run_number is None and args.source):
            run_number = os.path.basename(args.source).split('.')[4]

        if (args.Source):
            with open(args.Source, encoding="utf-8") as f:
                lines = f.readlines()
            if (args.run_number is not None):
                lines = [line for line in lines if self.has_run(line, args.run_number)]
            if (not lines):
                logger.error(f'no matching files in input file list {args.Source}')
                sys.exit(1)
            if (input_dsid is None):
                input_dsid = lines[0].split('.')[2]
            if (args.nfiles):
                lines = lines[:args.nfiles]
            return input_dsid, run_number, ''.join(lines)

        if (args.source):
            return input_dsid, run_number, None

        return input_dsid, run_number, self.raw_files(args)

#------------------------------------------------------------------------------
# field 4 (0-based) of the dot-separated name is formatted as %6r_%6s:
# the first six characters are the run number, the rest is the subrun
#------------------------------------------------------------------------------
    def has_run(self,path,run_number):
        fields = os.path.basename(path.strip()).split('.')
        return len(fields) > 4 and fields[4][:6] == f'{run_number:06d}'

#------------------------------------------------------------------------------
# list the raw data files of the run
#------------------------------------------------------------------------------
    def raw_files(self,args):
        data_dir = args.data_dir or os.getenv('RAW_DATA_DIR')
        if (not data_dir or args.run_number is None):
            logger.error('listing raw files needs --run_number and --data-dir (or $RAW_DATA_DIR)')
            sys.exit(1)

        files = sorted(f for f in glob.glob(os.path.join(data_dir, '*'))
                       if self.has_run(f, args.run_number))
        if (args.nfiles):
            files = files[:args.nfiles]
        if (not files):
            logger.error(f'no raw files for run {args.run_number} in {data_dir}')
            sys.exit(1)

        return ''.join(f'{f}\n' for f in files)

#------------------------------------------------------------------------------
# shell pipeline that prints the template fcl with the overrides applied
#------------------------------------------------------------------------------
    def overrides_cmd(self,args,template_fcl,input_dsid):
        q = shlex.quote

        cmd    = f'cat {q(template_fcl)}'
        subdir = f'{args.calib_run[0:3]}000' if args.calib_run else None

        # output ntuple name nts.mu2e.<input dataset ID>.<...>.stn: put the input dataset ID in the name
        if (input_dsid):
            cmd += ' | sed -E ' + q(f's|(histFileName *: *"nts\\.mu2e\\.)[^.]+(\\.)|\\1{input_dsid}\\2|')
        else:
            logger.warning('no input dataset ID, histFileName is left unchanged')

        if (args.calib_ver):
            ver = args.calib_ver
            if (ver == 'latest'):
                # read the latest calibration version of the run, it propagates to the file name as 'r100'
                fn = os.path.join(os.getenv('RUNDB_DIR'), 'rundb', subdir, args.calib_run, 'aaa_latest')
                with open(fn, encoding="utf-8") as f:
                    ver = f.readline().strip()
            cmd += ' | sed -E ' + q(f's|(calibration_set_).*(\\.fcl)|\\1{ver}\\2|')

        if (args.calib_run):
            cmd += ' | sed -E ' + q(f's|fcl/calibration_set|{subdir}/{args.calib_run}/calibration_set_{args.calib_run}|')

        if (args.odsid):
            cmd += ' | sed -E ' + q(f's|_odsid_|{args.odsid}|')

        return cmd

#------------------------------------------------------------------------------
# run a shell command, exit if it fails
#------------------------------------------------------------------------------
    def run_shell(self,cmd):
        p = subprocess.run(cmd,
                           executable="/bin/bash",
                           shell=True,
                           capture_output=True,
                           encoding="utf-8")
        if (p.returncode != 0):
            logger.error(f'rc={p.returncode} cmd:{cmd}\n{p.stderr}')
            sys.exit(1)
        return p.stdout

#------------------------------------------------------------------------------
if __name__ == "__main__":

    x    = SubmitJob();
    args = x.parse_parameters();

    x.run(args)
