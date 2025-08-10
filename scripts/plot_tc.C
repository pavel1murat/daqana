//
#include <format>

class PlotTC {
public:
  TFile* fFile;
  TTree* fTree;
  int    fRunNumber;

  PlotTC(int RunNumber, int RecoVersion=0, const char* Fn = "");
  PlotTC(const char* Fn);
  ~PlotTC();

  int plot_panel_edep(int Panel);
  int plot_panel_dt  (int Panel1, int Panel2);
  int plot_plane_dt  ();
};



//-----------------------------------------------------------------------------
// by default, sue datasets
//-----------------------------------------------------------------------------
PlotTC::PlotTC(int RunNumber, int RecoVersion, const char* Fn) {
  
  fRunNumber = RunNumber;

  std::string fn(Fn);

  if (fn == "") {
    std::string top_dir = "/projects/mu2e/data/projects/vst/datasets";
    std::string dsid    = std::format("nts.mu2e.trk.vst00s000r{:03d}n002.root",RecoVersion);
    fn                  = std::format("{}/{}/nts.mu2e.trk.vst00s000r{:03d}n002.{:06d}_000001.root",
                                      top_dir,dsid,RecoVersion,RunNumber);
  }
  
  std::string daqana_lib = Form("%s/lib/libdaqana_obj.so",gSystem->Getenv("SPACK_VIEW"));
  if (not gInterpreter->IsLoaded(daqana_lib.data())) {
    gInterpreter->Load(daqana_lib.data());
  }


  fFile = (TFile*) gROOT->GetListOfFiles()->FindObject(fn.data());
  if (fFile == nullptr) {
    fFile = TFile::Open(fn.data());
  }
  if (fFile == nullptr) {
    std::cout << "ERROR: coudn;t open the file:" << fn << std::endl;
  }
  else {
    fTree = (TTree*) fFile->Get("/MakeDigiNtuple/digis");
  }
}

//-----------------------------------------------------------------------------
PlotTC::PlotTC(const char* Fn) {
  fFile = (TFile*) gROOT->GetListOfFiles()->FindObject(Fn);
  if (fFile == nullptr) {
    fFile = TFile::Open(Fn);
  }
  if (fFile == nullptr) {
    std::cout << "ERROR: coudn;t open the file:" << Fn << std::endl;
  }
  else {
    fTree = (TTree*) fFile->Get("/MakeDigiNtuple/digis");
  }

  fRunNumber = -1;
}

//-----------------------------------------------------------------------------
PlotTC::~PlotTC() {
}

//-----------------------------------------------------------------------------
int PlotTC::plot_panel_edep(int Panel) {
  std::string hname(Form("h_edep_%02d_%06i",Panel,fRunNumber));

  std::string sel  = Form("run == %d && tc.nh_panel(%d)>1",fRunNumber,Panel);
  std::string hist = Form("%s(120,0,0.006)",hname.data());
  std::string cmd  = Form("tc.edep_panel(%d)",Panel);
    
  fTree->Draw(Form("%s>>%s",cmd.data(),hist.data()),sel.data());

  TH1F* h = (TH1F*) gROOT->FindObject(hname.data());
  return 0;
}


//-----------------------------------------------------------------------------
int PlotTC::plot_panel_dt(int Panel1, int Panel2) {
  std::string hname(Form("h_%02d%02d_%06i",Panel1,Panel2,fRunNumber));

  std::string sel  = Form("run == %d && tc.nh_panel(%d)>1 && tc.nh_panel(%d)>1",fRunNumber,Panel1,Panel2);
  std::string hist = Form("%s(250,-250,250)",hname.data());
  std::string cmd  = Form("tc.time_panel(%d)-tc.time_panel(%d)",Panel1,Panel2);
    
  fTree->Draw(Form("%s>>%s",cmd.data(),hist.data()),sel.data());

  TH1F* h = (TH1F*) gROOT->FindObject(hname.data());
  h->Fit("gaus","");
  return 0;
}

//-----------------------------------------------------------------------------
int PlotTC::plot_plane_dt() {
  std::string hname(Form("h_plane_dt_%06i",fRunNumber));
  
  std::string var  = Form("tc.timep(0)-tc.timep(1)");
  std::string sel  = Form("run==%d && tc.nh_panel(0)>1 && tc.nh_panel(1)>1",fRunNumber);

  fTree->Draw(Form("%s>>%s(250,-250,250)",var.data(),hname.data()),sel.data());

  TH1F* h = (TH1F*) gROOT->FindObject(hname.data());
  h->Fit("gaus","","",-50.,50.);
  return 0;
}
