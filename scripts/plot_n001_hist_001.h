//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 31 14:42:22 2026 by ROOT version 6.32.06
// from TTree digis/digis
// found on file: nts.mu2e.trk.vst00s000r000n001.120718_000001.root
//////////////////////////////////////////////////////////

#ifndef plot_n001_hist_h
#define plot_n001_hist_h

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

#include "daqana/obj/daqana_nt_format.hh"
  
  plot_n001_hist(int RunNumber);
  plot_n001_hist(const char* Fn = "nts.mu2e.trk.vst00s000r000n001.120718_000001.root");
  virtual ~plot_n001_hist();
  
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual bool     Notify();
  virtual void     Show(Long64_t entry = -1);

  void             Loop          (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents, double DtMin = 0, double DtMax = -1);
  void             LoopScan      (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents = -1);
  int              BookHistograms(Hist_t* Hist, int Mnid1, int Ch1, int Mnid2, int Ch2, float DtMin, float DtMax);
};

#endif
