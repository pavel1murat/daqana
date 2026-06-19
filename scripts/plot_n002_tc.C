///////////////////////////////////////////////////////////////////////////////
// run on digi_ntuple , process data of pulse injection runs
// for each straw, fill 4 histograms:
// - tcal : T(CAL),
// - bl   : baseline from the waveform fit
// - ph   : pulse height from the waveform fit
// - dt   : T(CAL)-T(HV)
// also for each panel plot "occupancy" histogram to see which channels
// have been pulsed in a given run
// for a given run, the histogtams are stored in the folder with the run number
// reference channel can be different for different panels and it also depends
// on the run number
// for each of 7-8 runs each panel has its own reference channel and then they will
// linked together
// 
// |------------+------------------+-------------|
// | run number | pulsed channels  | ref channel |
// |------------+------------------+-------------|
// |     121170 | 0, 4,  8, 12     |          12 |
// |     121171 | 1, 4,  9, 12     |          12 |
// |     121172 | 1, 5,  9, 13 ... |           9 |
// |     121173 | 2, 5, 10, 13     |          13 |
// |     121174 | missing          |             |
// |     121175 | 2, 6, 10, 14     |             |
// |     121176 | 3, 6, 11, 14     |             |
// |     121177 | 3, 7, 11, 15 ... |             |
// |------------+------------------+-------------|
// 
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include <format>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.
// #include "daqana/obj/DaqEvent.hh"
// #include "TObject.h"
// #include "daqana/obj/DaqStrawDigi.hh"
// #include "daqana/obj/DaqStrawHit.hh"
// #include "daqana/obj/DaqComboHit.hh"
// #include "daqana/obj/DaqTimeCluster.hh"
// #include "daqana/obj/DaqTrack.hh"
// #include "daqana/obj/DaqTrkStrawHit.hh"
// #include "daqana/obj/DaqSegment.hh"
#include "daqana/obj/TrkPanelMap_t.hh"

#include "daqana/ana/booking.hh"

class plot_n002_tc: public TNamed {
public :

  enum {
    kNStations         = 18,
    kNPanelsPerStation = 12,
  };
//-----------------------------------------------------------------------------
// fit results - data structures
//-----------------------------------------------------------------------------
  struct fit_result_t {
    double p[3];                        // gaussian fit
    double e[3];
    double chi2dof;
  };

//-----------------------------------------------------------------------------
// data structures
//-----------------------------------------------------------------------------
  struct PanelData_t {
    int    ref_channel;
    int    ref_tdc;                     // in units of 5./256 ns , approx 20 ps
  };

  struct TrackerData_t {
    PanelData_t   panel[kNStations][kNPanelsPerStation];
  } fTrackerData;
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct Hist_t {
    //    StationHist_t  station[kNStations];
    TH2F*          h_plane_dt;
    TH1F*          h_dt20;
  } fHist;
//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*       fTopFolder;
  TFolder*       fRunFolder;
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;

  fit_result_t   fFr[36]; 
  
  int            fRunNumber;
  int            fStation1;
  int            fStation2;
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

                                        // will the same directory work ?
#include "daqana_nt_format.hh"
                                        // for independent runs, the name should eb the same..
                                        // make it different to process the same run with different refence channels
  
  plot_n002_tc(const char* Name, int RunNumber, const char* Fn = nullptr);
  
  virtual ~plot_n002_tc();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  
  int              BookHistograms       (Hist_t* Hist, TFolder* Folder);

  int              CloneFitResults      ();

  int              ResetHistograms();
  int              SaveHistograms(const char* Filename);
  void             fit_plane_dt();
  
};

//-----------------------------------------------------------------------------
plot_n002_tc::plot_n002_tc(const char* Name, int RunNumber, const char* Fn = 0) : TNamed(Name,Name), fChain(0) {

  std::string dir("/data/mu2e/mu2etrk/datasets/vst00s000r000n002"); 
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000n002.{:06d}_000001.root",dir,RunNumber);
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
    if (!f || !f->IsOpen()) {
      f = new TFile(fn.data());
    }
  }

  // plane 20: MN 122, 080, 133, 150, 199, 278
  // plane 21: MN 112, 072, 202, 193, 200, 270
  // plane 22: MN 053, 052, 169, 035, 039, 043
  // plane 23: MN 151, 254, 225, 229, 275, 274
  // plane 24: MN 132, 145, 082, 207, 139, 155
  // plane 25: MN 055, 069, 064, 062, 063, 067

  fRunNumber = RunNumber;

  fTopFolder = (TFolder*) gROOT->GetRootFolder()->FindObject(Name);
  
  if (fTopFolder == nullptr) {
    fTopFolder = gROOT->GetRootFolder()->AddFolder(Name,Name);
  }

  std::string rns = std::to_string(fRunNumber);
  
  fRunFolder = fTopFolder->AddFolder(rns.data(),rns.data());
//-----------------------------------------------------------------------------
// allow histograms in different folders to have the same name
//-----------------------------------------------------------------------------
  TH1::AddDirectory(0);
  
  fTpm       = TrkPanelMap_t::Instance(RunNumber);
  fBook      = new Booking(fRunFolder);

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  BookHistograms(&fHist,fRunFolder);
}

//-----------------------------------------------------------------------------
plot_n002_tc::~plot_n002_tc() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n002_tc::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;
  
  name  = std::format("h_plane_dt");
  title = std::format("run:{} {} : plane dt",fRunNumber,prefix);
  fBook->HBook2F(Hist->h_plane_dt,name.data(),title.data(),500,-250,250,36,0,36,Folder);
  
  name  = std::format("h_dt20_{}",GetName());
  title = std::format("run:{} : plane dt20",fRunNumber);
  fBook->HBook1F(Hist->h_dt20,name.data(),title.data(),36,0,36,Folder);
  
  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n002_tc::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n002_tc::Init(TTree *tree) {
#include "daqana/scripts/daqana_nt_init.C"
}

//-----------------------------------------------------------------------------
Long64_t plot_n002_tc::LoadTree(Long64_t entry) {
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
void plot_n002_tc::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;
   
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // std::cout << std::format("----------- event:{:8} ntc:{:6}\n",evn,ntc);
//-----------------------------------------------------------------------------
// event initialization
//-----------------------------------------------------------------------------
    for (int itc=0; itc<ntc; itc++) {
      // use only 'good' time clusters
      if (tc_ngh[itc] == 0)                                 continue;
      
      for (int stn=0; stn<18; stn++) {
        for (int ip=0; ip<2; ip++) {
          int pln = 2*stn+ip;
                                        // ask for at least two hits in the plane
          if (tc__nhp[itc][stn][ip] > 1) {
                                        // use reference plane : slot 10 plane 0 (20)
            if (tc__nhp[itc][10][0] > 1) {
              float dt = tc__timep[itc][stn][ip]-tc__timep[itc][10][0];
              fHist.h_plane_dt->Fill(dt,pln);
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
int plot_n002_tc::ResetHistograms() {
  fHist.h_plane_dt->Reset();
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n002_tc::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
void plot_n002_tc::fit_plane_dt() {

  // auto x385 = new plot_n002_tc("x385",122385);
  // x385->Loop();

  for (int i=0; i<36; i++) {
    std::cout << std::format(".... fitting projection i:{}\n",i);
    
                                        // initialization
    for (int ip=0; ip<3; ip++) {
      fFr[i].p[ip] = 0.;
      fFr[i].e[ip] = 0.;
    }
    fFr[i].chi2dof = -1.;
    
    TH1D* hp =  fHist.h_plane_dt->ProjectionX(Form("_px%02d",i),i+1,i+1);
    if (hp->Integral() < 10) {
      continue;
    }
    
    TFitResultPtr tfr = hp->Fit("gaus","s","");

    if ((! tfr->IsValid()) or tfr->IsEmpty()) {
      printf("# FIT ERROR: channel: %2i \n",i);
      continue;
    }

    fFr[i].chi2dof = tfr->Chi2()/tfr->Ndf();
    double sf     = sqrt(fFr[i].chi2dof);

    for (int ip=0; ip<3; ip++) {
      fFr[i].p[ip] = tfr->Parameter(ip);
      if (i == 20) fFr[i].e[ip] = 0.;
      else         fFr[i].e[ip] = tfr->Error(ip)*sf;
    }
  }

  for (int i=0; i<36; i++) {
    fHist.h_dt20->SetBinContent(i+1,fFr[i].p[1]);
    fHist.h_dt20->SetBinError  (i+1,fFr[i].e[1]);
  }
 
  return;
}
