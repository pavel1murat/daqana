///////////////////////////////////////////////////////////////////////////////
// for now, process one station at a time
// station_time_calib :
// 1. fills edep histograms, masks off 'bad' panels: histograms with small number of entries
//     0x1: N(entries) = 0  ; 0x2: N(entries) < 200
// 2. for 'good' panels fills histograms in delta_t(i,j) = T(i)-T(j), i<j,
//    where T(i) and T(j) are the average per panel hit times within the time cluster
//    - optionally fit.
// 3. look at the histograms and determine the needed per-panel offsets
// 4. at this stage, don't need accuracy better than 1 ns
///////////////////////////////////////////////////////////////////////////////
/* Example:
root [0] .L v001/daqana/scripts/plot_tc.C
root [1] auto x = new PlotTC(120950,1,"/data/mu2e/mu2etrk/datasets/vst00s001r000n002/nts.mu2e.trk.vst00s001r000n002.120950_000001.root")
root [2] x->station_time_calib(11)
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
station:11 panel: 0 mnid:MN169 nent:  22909 panel_mask:0x0000
station:11 panel: 1 mnid:MN039 nent:     31 panel_mask:0x0002
station:11 panel: 2 mnid:MN043 nent:      0 panel_mask:0x0001
station:11 panel: 3 mnid:MN053 nent:  20669 panel_mask:0x0000
station:11 panel: 4 mnid:MN052 nent:  21283 panel_mask:0x0000
station:11 panel: 5 mnid:MN035 nent:  22479 panel_mask:0x0000
station:11 panel: 6 mnid:MN207 nent:  43942 panel_mask:0x0000
station:11 panel: 7 mnid:MN139 nent:  23577 panel_mask:0x0000
station:11 panel: 8 mnid:MN155 nent:  20922 panel_mask:0x0000
station:11 panel: 9 mnid:MN145 nent:  14917 panel_mask:0x0000
station:11 panel:10 mnid:MN132 nent:  26675 panel_mask:0x0000
station:11 panel:11 mnid:MN082 nent:  18300 panel_mask:0x0000

*** iteration 1: add 5 ns to T05; add 6 ns to T00; subtract 7 from T02       
|----------+----------+----------+----------+----------+----------+----------+----------+----------+--------+----------+----------+----------+-------|
| i1\i2    | 00:mn207 | 01:mn139 | 02:mn155 | 03:mn145 | 04:mn132 | 05:mn082 | 06:mn064 | 07:mn063 | 08:067 | 09:mn055 | 10:mn069 | 11:mn062 |  mean |
|----------+----------+----------+----------+----------+----------+----------+----------+----------+--------+----------+----------+----------+-------|
| 00:mn207 |          |      1.1 |      0.3 |          |          |      1.1 |          |      0.5 |    1.1 |      0.1 |      1.6 |      0.7 |  0.81 |
| 01:mn139 |          |          |      0.5 |          |          |     -1.5 |          |     -0.9 |   -0.5 |      0.2 |      0.1 |      0.6 | -0.21 |
| 02:mn155 |          |          |          |          |          |      0.2 |          |     -1.5 |   -1.1 |      0.4 |      0.2 |     -0.1 | -0.32 |
| 03:mn145 |          |          |          |          |          |          |          |          |        |          |          |          |  0.00 |
| 04:mn132 |          |          |          |          |          |          |          |          |        |          |          |          |  0.00 |
| 05:mn082 |          |          |          |          |          |          |          |     -0.7 |    0.1 |     -1.3 |      0.7 |     -1.6 | -0.56 |
| 06:mn064 |          |          |          |          |          |          |          |          |        |          |          |          |  0.00 |
| 07:mn063 |          |          |          |          |          |          |          |          |    0.1 |     -0.2 |      1.4 |     -0.2 |  0.28 |
| 08:mn067 |          |          |          |          |          |          |          |          |        |     -1.0 |      0.8 |      0.0 | -0.07 |
| 09:mn055 |          |          |          |          |          |          |          |          |        |          |      1.8 |      1.3 |  1.55 |
| 10:mn069 |          |          |          |          |          |          |          |          |        |          |          |      0.1 |  0.10 |
|----------+----------+----------+----------+----------+----------+----------+----------+----------+--------+----------+----------+----------+-------|
| #ERROR   |     0.00 |     1.10 |     0.40 |     0.00 |     0.00 |    -0.07 |     0.00 |    -0.65 |  -0.06 |    -0.30 |     0.94 |     0.10 |  0.14 |
|----------+----------+----------+----------+----------+----------+----------+----------+----------+--------+----------+----------+----------+-------|
#+TBLFM: @13=vmean(@2..@12);%.2f::$14=vmean($2..$13);%.2f::

*/
#include <format>
#include "daqana/obj/TrkPanelMap_t.hh"

class PlotTC {
public:

  struct Hist_t {
    TH1F* h_edep[12];
    TH1F* h_dt  [12][12]; // dt_ij = t[i]-t[j], i<j, only upper triangle is populated 
  } fHist;
  
  TFile* fFile;
  TTree* fTree;
  int    fRunNumber;

  int    fPanelMask[12];

  TrkPanelMap_t*  fTpm;

  PlotTC(int RunNumber, int RecoVersion=0, const char* Fn = "");
  PlotTC(const char* Fn);
  ~PlotTC();

  TH1F* plot_panel_edep(int Station,int Panel);
                                        // by defalt, fit, however sometimes do not want that
  TH1F* plot_panel_dt  (int Station1, int Panel1, int Station2, int Panel2, int PerformFit = 1);
  int   plot_plane_dt  (int Station);
  
  int   station_time_calib(int Station);
};



//-----------------------------------------------------------------------------
// by default, sue datasets
//-----------------------------------------------------------------------------
PlotTC::PlotTC(int RunNumber, int RecoVersion, const char* Fn) {
  
  fRunNumber = RunNumber;

  std::string fn(Fn);

  if (fn == "") {
    //    std::string top_dir = "/projects/mu2e/data/projects/vst/datasets";
    // std::string dsid    = std::format("nts.mu2e.trk.vst00s001r{:03d}n002.root",RecoVersion);
    std::string top_dir = "/data/mu2e/mu2etrk/datasets";
    std::string dsid    = std::format("vst00s001r{:03d}n002",RecoVersion);
    fn                  = std::format("{}/{}/nts.mu2e.trk.vst00s001r{:03d}n002.{:06d}_000001.root",
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

  fTpm = TrkPanelMap_t::Instance(RunNumber);
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
  fTpm = nullptr;
}

//-----------------------------------------------------------------------------
PlotTC::~PlotTC() {
}

//-----------------------------------------------------------------------------
TH1F* PlotTC::plot_panel_edep(int Station, int Panel) {
  std::string hname(Form("h_edep_%02d_%06i",Panel,fRunNumber));

  std::string sel  = Form("run == %d && tc.nh_panel(%d,%d)>1",fRunNumber,Station,Panel);
  std::string hist = Form("%s(120,0,0.006)",hname.data());
  std::string cmd  = Form("tc.edep_panel(%d,%d)",Station,Panel);
    
  fTree->Draw(Form("%s>>%s",cmd.data(),hist.data()),sel.data());

  return (TH1F*) gROOT->FindObject(hname.data());
}


//-----------------------------------------------------------------------------
// within one station, panels are numbered from 0 to 11
//-----------------------------------------------------------------------------
TH1F* PlotTC::plot_panel_dt(int Station1, int Panel1, int Station2, int Panel2, int PerformFit) {
  std::string hname(Form("h_%02d%02d_%02d%02d_%06i",Station1,Panel1,Station2,Panel2,fRunNumber));

  std::string sel  = Form("run == %d && tc.nh_panel(%d,%d)>1 && tc.nh_panel(%d,%d)>1",fRunNumber,Station1,Panel1,Station2,Panel2);
  std::string hist = Form("%s(250,-250,250)",hname.data());
  std::string cmd  = Form("tc.time_panel(%d,%d)-tc.time_panel(%d,%d)",Station1,Panel1,Station2,Panel2);
    
  fTree->Draw(Form("%s>>%s",cmd.data(),hist.data()),sel.data());

  TH1F* h = (TH1F*) gROOT->FindObject(hname.data());
                                        // sometimes fitting needs to be more intelligent...
  if (PerformFit) h->Fit("gaus","",-50,50);
  
  return h;
}

//-----------------------------------------------------------------------------
// for a given 'Station' plots delta_t = T(plane=1)-T(plane=0)
//-----------------------------------------------------------------------------
int PlotTC::plot_plane_dt(int Station, int PerformFit) {
  std::string hname(Form("h_plane_dt_%06i",fRunNumber));
  
  std::string var  = Form("tc.timep(1)-tc.timep(0)");
  std::string sel  = Form("run==%d && tc.nh_panel(%d,1)>1 && tc.nh_panel(%d,0)>1",
                          fRunNumber,Station,Station);

  fTree->Draw(Form("%s>>%s(250,-250,250)",var.data(),hname.data()),sel.data());

  TH1F* h = (TH1F*) gROOT->FindObject(hname.data());
  if (PerformFit) h->Fit("gaus","","",-50.,50.);
  return 0;
}


//-----------------------------------------------------------------------------
// all numbering : offline
// 'Station' : offline slot number
// 'plane'   : 2*Station+i/6
// 'i'       : offline panel number in a plane
//-----------------------------------------------------------------------------
int PlotTC::station_time_calib(int Station) {
  int rc(0);

  int   nent[12];

  for (int i=0; i<12; i++) {
    fPanelMask[i] = 0;
    int offline_plane = 2*Station + i/6;        // 20 and 21 for Station=10
    int offline_panel = i % 6;
    
    TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(offline_plane,offline_panel);
    
    fHist.h_edep[i] = plot_panel_edep(Station,i);
    nent        [i] = fHist.h_edep[i]->Integral();
//-----------------------------------------------------------------------------
// just in case, try to distinguish between panels which have been excluded from the readout
// and panels which have been read out but had low HV
//-----------------------------------------------------------------------------
    if      (nent[i] ==  0) fPanelMask[i] = 0x1;
    else if (nent[i] < 200) fPanelMask[i] = 0x2;
    
    std::cout << std::format("station:{:2d} panel:{:2d} mnid:MN{:03d} nent:{:7d} panel_mask:0x{:04x}\n",
                             Station,i,tpmd->mnid,nent[i],fPanelMask[i]);
  }

  std::cout << std::format("... fill dt histograms ...\n");
  
  for (int i1=0; i1<12; i1++) {
                                        // reset all dt histograms 
    for (int i2=0; i2<12; i2++) {
      fHist.h_dt[i1][i2] = nullptr;
    }
    
    if (fPanelMask[i1] != 0) continue;

    for (int i2=i1+1; i2<12; i2++) {
      if (fPanelMask[i2] != 0) continue;
                                        // do not fit, just fill the histograms
      int perform_fit(0);
      fHist.h_dt[i1][i2] = plot_panel_dt(Station, i1, Station, i2, perform_fit);
    }
  }
 
  return rc;
}
