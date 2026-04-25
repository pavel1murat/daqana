#ifndef plot_n001_hist_h
#define plot_n001_hist_h

///////////////////////////////////////////////////////////////////////////////
/*
  gSystem->Load("v001/.spack-env/view/lib/libdaqana_obj.so");
  .L v001/daqana/scripts/plot_n001_hist_001.C
  auto x07 = new plot_n001_hist(120807);
  x07->Loop(151,44,229,44,20000);
  
*/
#include <format>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.
#include "daqana/obj/DaqEvent.hh"
#include "TObject.h"
#include "daqana/obj/DaqStrawDigi.hh"
#include "daqana/obj/DaqStrawHit.hh"
#include "daqana/obj/DaqComboHit.hh"
#include "daqana/obj/DaqTimeCluster.hh"
#include "daqana/obj/DaqTrack.hh"
#include "daqana/obj/DaqTrkStrawHit.hh"
#include "daqana/obj/DaqSegment.hh"

class plot_n001_hist {
public :
  
  struct Hist_t {
    TH2F* fDtVsEvt;
    TH1F* fDt10k;
    TH1F* fNhCh1;
    TH1F* fNhCh2;
    TH1F* fGenDt1;
    TH1F* fGenDt2;
  };

  Hist_t         fHist;
  int            fRunNumber;
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // will the same directory work ?
#include "daqana_nt_format.hh"
  
  plot_n001_hist(int RunNumber);
  plot_n001_hist(const char* Fn = "nts.mu2e.trk.vst00s000r000n001.120718_000001.root");
  virtual ~plot_n001_hist();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents, double DtMin = 0, double DtMax = -1);
  void             LoopScan      (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents = -1);
  int              BookHistograms(Hist_t* Hist, int Mnid1, int Ch1, int Mnid2, int Ch2, float DtMin, float DtMax);
};

#endif

//-----------------------------------------------------------------------------
plot_n001_hist::plot_n001_hist(int RunNumber) : fChain(0) {

  std::string dir("/data/mu2e/mu2etrk/datasets/vst00s000r000n001"); 
  std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000n001.{:06d}_000001.root",dir,RunNumber);
  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
  }

  // plane 20: MN 122, 080, 133, 150, 199, 278
  // plane 21: MN 112, 072, 202, 193, 200, 270
  // plane 22: MN 053, 052, 169, 035, 039, 043
  // plane 23: MN 151, 254, 225, 229, 275, 274
  // plane 24: MN 132, 145, 082, 207, 139, 155
  // plane 25: MN 055, 069, 064, 062, 063, 067

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  Init(tree);

  fRunNumber         = RunNumber;
  TrkPanelMap_t* fTpm = TrkPanelMap_t::Instance(RunNumber);

  fHist.fDtVsEvt = nullptr;
}
//-----------------------------------------------------------------------------
plot_n001_hist::plot_n001_hist(const char* Fn) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  std::string fn = std::format("/data/mu2e/mu2etrk/datasets/vst00s000r000n001/{}",Fn);
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
  }

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  Init(tree);

  fHist.fDtVsEvt = nullptr;
}


//-----------------------------------------------------------------------------
plot_n001_hist::~plot_n001_hist() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n001_hist::BookHistograms(Hist_t* Hist, int Mnid1, int Ch1, int Mnid2, int Ch2, float DtMin, float DtMax) {
  float dt_min(-100000.), dt_max(100000.);

  if (DtMax > DtMin) {
    dt_min = DtMin;
    dt_max = DtMax;
  }

  std::string prefix = std::format("h_r_{:06d}_mn{:03d}_{:02d}_mn{:03d}_{:02d}",fRunNumber,Mnid1,Ch1,Mnid2,Ch2);

  std::string name  = std::format("{}_dt_vs_evt",prefix);
  std::string title = std::format("run:{} {}",fRunNumber,name);
  Hist->fDtVsEvt = new TH2F(name.data(),title.data(),600,0,2400000,2000,dt_min,dt_max);
  Hist->fDtVsEvt->GetYaxis()->SetTitle("time, ns");

  name  = std::format("{}_dt10k",prefix);
  title = std::format("run:{} {}",fRunNumber,name);
  Hist->fDt10k = new TH1F(name.data(),title.data(),20000,0,20000);

  name  = std::format("{}_nh_ch1",prefix);
  title = std::format("run:{} {}",fRunNumber,name);
  Hist->fNhCh1 = new TH1F(name.data(),title.data(),50,0,50);

  name  = std::format("{}_nh_ch2",prefix);
  title = std::format("run:{} {}",fRunNumber,name);
  Hist->fNhCh2 = new TH1F(name.data(),title.data(),50,0,50);
  
  name  = std::format("{}_gen_dt1",prefix);
  title = std::format("run:{} {}",fRunNumber,name);
  Hist->fGenDt1 = new TH1F(name.data(),title.data(),2000,16300,16500);

  name  = std::format("{}_gen_dt2",prefix);
  title = std::format("run:{} {}",fRunNumber,name);
  Hist->fGenDt2 = new TH1F(name.data(),title.data(),2000,16300,16500);
  
  return 0;
}

//-----------------------------------------------------------------------------
void plot_n001_hist::Loop(int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents, double DtMin, double DtMax) {

   BookHistograms(&fHist,Mnid1,Ch1,Mnid2,Ch2,DtMin,DtMax);

   fHist.fDtVsEvt->Reset();

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   int nev = NEvents;
   if (NEvents <= 0) nev = nentries;
   
   for (int jentry=0; jentry<nev; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      std::vector<int> ind1;
      std::vector<int> ind2;

      for (int i=0; i<nsdtot; i++) {
        int ich = sd_sid[i] & 0x7f;
        if ((sd_mnid[i] == Mnid1) and (ich == Ch1)) {
          ind1.push_back(i);
        }
        if ((sd_mnid[i] == Mnid2) and (ich == Ch2)) {
          ind2.push_back(i);
        }
      }
      
      int nh_ch1 = ind1.size();
      int nh_ch2 = ind2.size();
      
      if ((nh_ch1 > 0) and (nh_ch2 > 0)) {
        int i1 = ind1[0];
        int i2 = ind2[0];
        
        float t1 = sd_tdc0[i1]*5./256;
        float t2 = sd_tdc0[i2]*5./256;
        float dt = t1-t2;                     // units of ns
        fHist.fDtVsEvt->Fill(evn,dt);
        fHist.fDt10k->Fill(evn,dt);
      }

      if (nh_ch1 > 1) {
        int i1 = ind1[0];
        int i2 = ind1[1];
        float t1 = sd_tdc0[i1]*5./256;
        float t2 = sd_tdc0[i2]*5./256;
        float dt = t2-t1;                     // units of ns
        fHist.fGenDt1->Fill(dt);
      }

      if (nh_ch2 > 1) {
        int i1 = ind2[0];
        int i2 = ind2[1];
        float t1 = sd_tdc0[i1]*5./256;
        float t2 = sd_tdc0[i2]*5./256;
        float dt = t2-t1;                     // units of ns
        fHist.fGenDt2->Fill(dt);
      }

      fHist.fNhCh1->Fill(nh_ch1);
      fHist.fNhCh2->Fill(nh_ch2);
   }
}

//-----------------------------------------------------------------------------
void plot_n001_hist::LoopScan(int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents) {
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   int nev = NEvents;
   if (NEvents <= 0) nev = nentries;
   
   for (int jentry=0; jentry<nev; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      std::vector<int> ind1;
      std::vector<int> ind2;

      for (int i=0; i<nsdtot; i++) {
        int ich = sd_sid[i] & 0x7f;
        if ((sd_mnid[i] == Mnid1) and (ich == Ch1)) {
          ind1.push_back(i);
        }
        if ((sd_mnid[i] == Mnid2) and (ich == Ch2)) {
          ind2.push_back(i);
        }
      }
      
      int nh_ch1 = ind1.size();
      int nh_ch2 = ind2.size();
      
      if ((nh_ch1 > 0) and (nh_ch2 > 0)) {
        int i1 = ind1[0];
        int i2 = ind2[0];
        
        float t1 = sd_tdc0[i1]*5./256;
        float t2 = sd_tdc0[i2]*5./256;
        float dt = t1-t2;                     // units of ns

        std::cout << std::format(" {:6d}:{:06d}:{:8d}",run,srn,evn)
                  << std::format("{:2d} {:2d} {:10.3f} {:10.3f} {:10.3f}",nh_ch1,nh_ch2,t1,t2,dt)
                  << std::endl;
      }
   }
}

//-----------------------------------------------------------------------------
Int_t plot_n001_hist::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
Long64_t plot_n001_hist::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

#include "daqana/scripts/daqana_nt_init.C"
