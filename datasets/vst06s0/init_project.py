#!/usr/bin/python

from local_classes import *
# from mixing_inputs import *

class Project(ProjectBase):
    def init_datasets(self):
#------------------------------------------------------------------------------
# datasets of this family
#-------v----------------------------------------------------------------------
        self.add_dataset(Dataset('raw.mu2e.trk.vst.art'               ,'vst06s0s00r0000','local'))
        self.add_dataset(Dataset('rec.mu2e.vst06s0s10r0000.daqana.art','vst06s0s10r0000','local'))
#------------------------------------------------------------------------------
# a job always has an input dataset, but...
#------------------------------------------------------------------------------
        self.fInputDsID = None;
        if (self.fIDsID) : self.fInputDataset = self.fDataset[self.fIDsID];


    def __init__(self,idsid=None):
        
        ProjectBase.__init__(self,project='daqana',family_id='vst06s0',idsid=idsid);
        self.init_datasets();
#------------------------------------------------------------------------------
# stage 1
# -------
# s1:cosmic_reco : InputDsID is 'vst06s0s00r0000'
#                  reconstruction job has only one output stream
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s1');
        input_dsid                   = idsid;
        if (input_dsid == None): input_dsid = 'vst06s0s00r0000'
        
        job                          = s.new_job('cosmic_reco',input_dsid);

        job.fNInputFiles             = -1                     # number of segments defined by the input dataset
             
        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '3h'   
        job.fIfdh                    = 'ifdh' ## 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        # job.fInputDataset.print()
            
        output_stream                = job.fInputDataset.output_stream()

        odsid                        = self.fFamilyID+s.name()+output_stream+'r0000'; # i.e. +'s10r0000'

        job.fOutputStream            = ['defaultOutput'                ]
        job.fOutputDsID              = [odsid                          ]
        job.fOutputFnPattern         = ['rec.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = ['art'                          ]
#------------------------------------------------------------------------------
# s1:make_dgn : ntupling job has only one output stream
#------------------------------------------------------------------------------        
        if (input_dsid == None): input_dsid = 'vst06s0s10r0000'
        job                          = s.new_job('make_dgn',input_dsid);

        job.fNInputFiles             = -1                     # number of segments defined by the input dataset
             
        job.fMaxInputFilesPerSegment =  10
        # job.fNEventsPerSegment       =  100000
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '3h'   
        job.fIfdh                    = 'ifdh' ## 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        output_stream                = job.fInputDataset.output_stream()

        odsid                        = self.fFamilyID+s.name()+output_stream+'r0000';

        job.fOutputStream            = ['InitStntuple'     ]
        job.fOutputDsID              = [odsid              ]
        job.fOutputFnPattern         = ['nts.murat.'+odsid ]
        job.fOutputFormat            = ['root'             ]
#------------------------------------------------------------------------------
# s1:make_dgn02 : ntupling job has only one output stream
#------------------------------------------------------------------------------        
        if (input_dsid == None): input_dsid = 'vst06s0s10r0000'
        job                          = s.new_job('make_dgn02',input_dsid);

        job.fNInputFiles             = -1                     # number of segments defined by the input dataset
             
        job.fMaxInputFilesPerSegment =  10
        # job.fNEventsPerSegment       =  100000
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '3h'   
        job.fIfdh                    = 'ifdh' ## 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        output_stream                = job.fInputDataset.output_stream()

        odsid                        = self.fFamilyID+s.name()+output_stream+'r0002';

        job.fOutputStream            = ['InitStntuple'     ]
        job.fOutputDsID              = [odsid              ]
        job.fOutputFnPattern         = ['nts.murat.'+odsid ]
        job.fOutputFormat            = ['root'             ]
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------
