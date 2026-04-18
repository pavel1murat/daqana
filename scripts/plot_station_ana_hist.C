//-----------------------------------------------------------------------------
// plot histograms made by StationAna module
// plot by plane - 6 panels at a time
// c->Print("trk_fragment_ana_run_002386.figure_0001.pdf")
//
// gSystem->Load("v001/.spack-env/view/lib/libdaqana_obj.so");
// 
//-----------------------------------------------------------------------------
#include <format>

#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"

#include "daqana/obj/TrkPanelMap_t.hh"

// #include "trace.h"
// #define TRACE_NAME "plot_station_ana_hist"

int MaxNHits   (-1);
int MaxChannel (-1);
int MaxEvent   (-1);
int MaxY       (-1);

namespace {
  //  const char* HistDir   = "/exp/mu2e/data/projects/tracker/vst/hist" ;
  //  const char* HistDir   = "/data/tracker/vst/hist" ;
  const char*    HistDir   = "./" ;        // local
  int            _figure   = 0;
  TCanvas*       _canvas(nullptr);

  TrkPanelMap_t* _trkPanelMap (nullptr);
};

// assume histogramming starter from the first fine (first subrun=1)
TFile* open_file(const char* Fn, int RunNumber) {
  TFile* f(nullptr);
  
  if (Fn == nullptr) f = TFile::Open(Form("%s/hst.mu2e.vst00s000r000n000.make_station_hist.%06i_000001.root",HistDir,RunNumber));
  else               f = TFile::Open(Fn);

  return f;
}

//-----------------------------------------------------------------------------
int plot_edep(int RunNumber, int Slot, const char* Fn = nullptr) {

  TFile* f = open_file(Fn,RunNumber);

  TrkPanelMap_t* tpm = TrkPanelMap_t::Instance(RunNumber);

  std::string cname = std::format("c_slot_{:02d}",Slot);
  TCanvas* c = new TCanvas(cname.data(),cname.data(),1800,1050);
  c->Divide(4,3);

  for (int pln=0; pln<2; pln++) {
    int plane = 2*Slot+pln;
    for (int panel=0; panel<6; panel++) {
      TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
      int mnid = tpmd->mnid;
      // TLOG(TLVL_DEBUG) << std::format("Slot:{:02d} Plane:{} panel:{:2d} mnin:MN{:03d}",Slot,Plane,plane,mnid);
      std::cout << std::format("Slot:{:02d} plane:{} panel:{:2d} plane:{} mnid:MN{:03d}\n",Slot,pln,plane,panel,mnid);
      
      c->cd(6*pln+panel+1);
      gPad->SetLogy(1);
      std::string hname = std::format("//StationAna/slot_{:02d}/MN{:03d}/edep",Slot,mnid);
      TH1* h1 = (TH1*) f->Get(hname.data());
      if (MaxY > 0) h1->GetYaxis()->SetRangeUser(0.1,MaxY);
      h1->Draw();
    }
  }

  _canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
int plot_occupg(int RunNumber, int Slot, const char* Fn = nullptr) {

  TFile* f = open_file(Fn,RunNumber);

  TrkPanelMap_t* tpm = TrkPanelMap_t::Instance(RunNumber);

  std::string cname = std::format("c_slot_{:02d}",Slot);
  TCanvas* c = new TCanvas(cname.data(),cname.data(),1800,1050);
  c->Divide(4,3);

  for (int pln=0; pln<2; pln++) {
    int plane = 2*Slot+pln;
    for (int panel=0; panel<6; panel++) {
      TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
      int mnid = tpmd->mnid;
      // TLOG(TLVL_DEBUG) << std::format("Slot:{:02d} Plane:{} panel:{:2d} mnin:MN{:03d}",Slot,Plane,plane,mnid);
      std::cout << std::format("Slot:{:02d} plane:{} panel:{:2d} plane:{} mnid:MN{:03d}\n",Slot,pln,plane,panel,mnid);
      
      c->cd(6*pln+panel+1);
      gPad->SetLogy(1);
      std::string hname = std::format("//StationAna/slot_{:02d}/MN{:03d}/occupg",Slot,mnid);
      TH1* h1 = (TH1*) f->Get(hname.data());
      if (MaxY > 0) h1->GetYaxis()->SetRangeUser(0.1,MaxY);
      h1->Draw();
    }
  }

  _canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
// if Fn = nullptr, use defaults
// Plane = 0 or 1
//-----------------------------------------------------------------------------
int plot(int RunNumber, int Figure, int Slot, const char* Fn = nullptr, int Print = 0) {
  int rc(0);
  
  _figure = Figure;
  
  if      (Figure == 1) rc = plot_edep  (RunNumber, Slot, Fn);
  else if (Figure == 2) rc = plot_occupg(RunNumber, Slot, Fn);
  else {
    std::cout << std::format("ERROR: undefined figure {}\n",Figure);
    return -1;
  }

  if ((rc == 0) and (Print > 0)) {
    // everything looks ok, print
    char fn[1000];
    sprintf(fn,"station_ana_run_%06i.figure_%05i.pdf",RunNumber,Figure);
    _canvas->Print(fn);
  }

  return rc;
}
