#ifndef plot_n002_hist_h
#define plot_n002_hist_h

///////////////////////////////////////////////////////////////////////////////
/*
  gSystem->Load("v001/.spack-env/view/lib/libdaqana_obj.so");
  .L v001/daqana/scripts/plot_n002_hist_001.C
  auto x07 = new plot_n002_hist(120807);
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
#include "daqana/obj/TrkPanelMap_t.hh"

#include "booking.C"

class plot_n002_hist { 
public :
  
  struct StrDigiHist_t {
    TH1F* h_tcal;
    TH1F* h_bl;
    TH1F* h_ph;
    TH1F* h_dt;
  };
  
  struct PanelHist_t {
    StrDigiHist_t  ch[96];
    TH1F* h_occup;
  };
  
  struct StationHist_t {
    PanelHist_t  panel[12];
  };
  
  struct Hist_t {
    StationHist_t  station[18];
  } fHist;

  TFolder*       fTopFolder;
  TFolder*       fRunFolder;
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;
  
  int            fRunNumber;
  int            fStation1;
  int            fStation2;
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

                                        // will the same directory work ?
#include "daqana_nt_format.hh"
  
  plot_n002_hist(int RunNumber, int Station1, int Station2);
  
  plot_n002_hist(const char* Fn = "nts.mu2e.trk.vst00s000r000n001.120718_000001.root");

  virtual ~plot_n002_hist();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  
  int              BookStrDigiHistograms(StrDigiHist_t* Hist, int Station, int Panel, int Channel, TFolder* Folder);
  int              BookPanelHistograms  (PanelHist_t*   Hist, int Station, int Panel, TFolder* Folder);
  int              BookStationHistograms(StationHist_t* Hist, int Station, TFolder* Folder);
  int              BookHistograms       (Hist_t* Hist, TFolder* Folder);

  int              FillStrDigiHistograms(StrDigiHist_t* Hist, int Ind);
  int              FillPanelHistograms  (PanelHist_t*   Hist, int Ind);

  int              ResetHistograms();
  int              SaveHistograms(const char* Filename);
};

//-----------------------------------------------------------------------------
plot_n002_hist::plot_n002_hist(int RunNumber, int Station1, int Station2) : fChain(0) {

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

  fRunNumber = RunNumber;
  fStation1  = Station1;
  fStation2  = Station2;

  fTopFolder = (TFolder*) gROOT->GetRootFolder()->FindObject("n002");
  if (fTopFolder == nullptr) {
    fTopFolder = gROOT->GetRootFolder()->AddFolder("n002","n002");
  }

  std::string rns = std::to_string(fRunNumber);
  
  fRunFolder = fTopFolder->AddFolder(rns.data(),rns.data());
  
  std::string folder_name = std::format("n002_{:6d}",fRunNumber);
//-----------------------------------------------------------------------------
// allow histograms in different folders to have the same name
//-----------------------------------------------------------------------------
  TH1::AddDirectory(0);
  
  fTpm       = TrkPanelMap_t::Instance(RunNumber);
  fBook      = new Booking(fRunFolder);

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  Init(tree);
  BookHistograms(&fHist,fRunFolder);

}


//-----------------------------------------------------------------------------
plot_n002_hist::plot_n002_hist(const char* Fn) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  std::string fn = std::format("/data/mu2e/mu2etrk/datasets/vst00s000r000n001/{}",Fn);
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
  }

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  Init(tree);

  BookHistograms(&fHist,fRunFolder);
}


//-----------------------------------------------------------------------------
plot_n002_hist::~plot_n002_hist() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n002_hist::BookStrDigiHistograms(StrDigiHist_t* Hist, int Station, int Panel, int Channel, TFolder* Folder) {
  
  std::string prefix = std::format("s{:02d}_p{:02d}_{:02d}",Station,Panel,Channel);
  std::string name, title;
  
  name  = std::format("h_{}_tcal",prefix);
  title = std::format("run:{} {} : tcal",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_tcal,name.data(),title.data(),1000,0,100000,Folder);
  
  name  = std::format("h_{}_bl",prefix);
  title = std::format("run:{} {} : bl",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,1000,Folder);
  
  name  = std::format("h_{}_ph",prefix);
  title = std::format("run:{} {} : ph",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),500,0,1000,Folder);
  
  name  = std::format("h_{}_dt",prefix);
  title = std::format("run:{} {} : T(cal)-T(hv)",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_dt,name.data(),title.data(),500,-25,25,Folder);
  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookPanelHistograms(PanelHist_t* Hist, int Station, int Panel, TFolder* Folder) {

  std::string prefix = std::format("s{:02d}_p{:02d}",Station,Panel);
  std::string name, title;
  
  name  = std::format("h_{}_occup",prefix);
  title = std::format("run:{} {} : occupancy",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_occup,name.data(),title.data(),96,0,96,Folder);
  
  for (int i=0; i<96; i++) {
    std::string folder_name = std::format("ch_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookStrDigiHistograms(&Hist->ch[i],Station,Panel,i,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookStationHistograms(StationHist_t* Hist, int Station, TFolder* Folder) {

  for (int pnl=0; pnl<12; pnl++) {
    std::string folder_name = std::format("pnl_{:02d}",pnl);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookPanelHistograms(&Hist->panel[pnl],Station,pnl,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  for (int i=fStation1; i<=fStation2; i++) {
    std::string folder_name = std::format("stn_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookStationHistograms(&Hist->station[i],i,fol);
  }
 
  return 0;
}

//-----------------------------------------------------------------------------
// 'Ind' : straw digi index in an array
//-----------------------------------------------------------------------------
int plot_n002_hist::FillStrDigiHistograms(StrDigiHist_t* Hist, int Ind) {
  Hist->h_ph->Fill(sd_ph[Ind]);
  Hist->h_bl->Fill(sd_bl[Ind]);
  Hist->h_tcal->Fill(sd_tdc0[Ind]);

  float dt = (sd_tdc0[Ind]-sd_tdc1[Ind])*5./256.;

  Hist->h_dt->Fill(dt);

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::FillPanelHistograms(PanelHist_t* Hist, int Ind) {

  int straw = sd_sid[Ind] & 0x7f;
  Hist->h_occup->Fill(straw);

  return 0;
}


//-----------------------------------------------------------------------------
Int_t plot_n002_hist::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n002_hist::Init(TTree *tree) {
#include "daqana/scripts/daqana_nt_init.C"
}

//-----------------------------------------------------------------------------
Long64_t plot_n002_hist::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

//-----------------------------------------------------------------------------
void plot_n002_hist::Loop(int NEvents) {

  // fHist.fDtVsEvt->Reset();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;
   
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // std::cout << std::format("----------- event:{:8} nsdtot:{:6}\n",evn,nsdtot);

    for (int i=0; i<nsdtot; i++) {
      TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_mnid(sd_mnid[i]);
      if (tpmd == nullptr) {
        std::cout << std::format("event:{:08d} ERROR: for digi #{} : wrong MNID:{}, SKIP\n",evn,i,sd_mnid[i]);
        continue;
      }
      int ipln = tpmd->plane;
      int ist  = ipln / 2;
      int ipnl = tpmd->panel + 6*(ipln % 2);

      if ((ist < fStation1) or (ist > fStation2)) continue;
      
      int ich  = sd_sid[i] & 0x7f;
        
      PanelHist_t*   panel_hist = &fHist.station[ist].panel[ipnl];
      StrDigiHist_t* sd_hist    = &fHist.station[ist].panel[ipnl].ch[ich];
        
      // std::cout << std::format("----------- ist:{} ipnl:{} ich:{}\n",ist,ipnl,ich);
      if (sd_ph[i] > 50.) {
        FillPanelHistograms  (panel_hist,i);
        FillStrDigiHistograms(sd_hist   ,i);
      }
    }
  }
}

//-----------------------------------------------------------------------------
int plot_n002_hist::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n002_hist::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

#endif


//-----------------------------------------------------------------------------
int make_calib_hist(int Run1, int Run2, int Station1, int Station2) {
  int k = 0;
  plot_n002_hist* p[100]; // not more than 8 calibration runs in a series
  for (int run=Run1; run<Run2+1; run++) {
    p[i] = new plot_n002_hist(run,Station1,Station2);
    p[i]->Loop();
  }

  std::string fn = std::format("pulse_injection_{:6d}_{:6d}.hist",Run1,Run2);
  p[0]->SaveHistograms(fn.data());
  
  return 0;
}
