#ifndef __daqana_ana_plot_n001_noise_hh__
#define __daqana_ana_plot_n001_noise_hh__

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
#include "daqana/obj/DaqEvent.hh"

// #include "TObject.h"
// #include "daqana/obj/DaqStrawDigi.hh"
// #include "daqana/obj/DaqStrawHit.hh"
// #include "daqana/obj/DaqComboHit.hh"
// #include "daqana/obj/DaqTimeCluster.hh"
// #include "daqana/obj/DaqTrack.hh"
// #include "daqana/obj/DaqTrkStrawHit.hh"
// #include "daqana/obj/DaqSegment.hh"
#include "daqana/obj/TrkPanelMap_t.hh"

#include "ana/ana/booking.hh"

class plot_n001_noise: public TNamed {
public :

  enum {
    kNStations         = 18,
    kNPlanes           = 36,
    kNPanelsPerStation = 12,
  };

  //-----------------------------------------------------------------------------
  // fit results - data structures
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // data structures
  //-----------------------------------------------------------------------------

  struct noise_t {
    int evmin;
    int evmax;
    double sumx;
    double sumx2;
    int nevents;

    noise_t(int FirstEvent = 0) {
      evmin   = FirstEvent;
      evmax   = 0;
      nevents = 0;
      sumx    = 0;
      sumx2   = 0;
    }

    double sigm() {
      double s(0);
      if (nevents > 0) {
        s = sumx2/nevents - (sumx/nevents)*(sumx/nevents);
      }
      return s;
    }
    
  };

  std::vector<noise_t> dat; // for one panel;

  std::vector<noise_t>*    fDat[500];   // sparse

  std::vector<std::vector<noise_t>*> fListOfDat;    // compact
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct Hist_t {
    TH1D* h_noise[12];
  } fHist;
  
//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*       fTopFolder;
  TFolder*       fRunFolder;
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;

  //  fit_result_t   fFr[36]; 
  
  int            fRunNumber;

  int            fRefChannel; // 21
  int            fMaxEvent;   // for X-axis truncation
  
  DaqStrawDigi*  fSdr[216];
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

                                        // will the same directory work ?
  DaqEvent*       fEvent; // #include "daqana_nt_format.hh"
  TClonesArray*   fSd;
                                        // for independent runs, the name should eb the same..
                                        // make it different to process the same run with different refence channels
  
  plot_n001_noise(const char* Name, int RunNumber, const char* Fn = nullptr);
  
  virtual ~plot_n001_noise();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  int              BookHistograms (Hist_t* Hist, TFolder* Folder);
  int              ResetHistograms();
  int              SaveHistograms (const char* Filename);
  void             fill_histograms();
};
#endif
