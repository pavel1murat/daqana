//

const char* Mu2eDatasetDir = "/data/tracker/vst";
//-----------------------------------------------------------------------------
TTree* get_digi_ntuple(const char* Dsid,int RunNumber) {
  TFile* f = TFile::Open(Form("%s/%s/nts.mu2e.trk.%s.%06i_000001.root",
                              Mu2eDatasetDir,Dsid,Dsid,RunNumber));
  TTree* t = (TTree*) f->Get("/MakeDigiNtuple/digis");
  return t;
}
