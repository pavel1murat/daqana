#ifndef __daqana_ana_plot_crv_tc_dt_hh__
#define __daqana_ana_plot_crv_tc_dt_hh__

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

class plot_crv_tc_dt: public TNamed {
public :

  enum {
    kNStations         = 18,
    kNPlanes           = 36,
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
  struct Index_t {
    int sel;
    int slot;                           // 0-17
    int plane;                          // offline
    int panel;                          // offline
    int pnl12;                          // panel index within the station (0-11)
    int mnid;
    int ch;
  };
  
  struct CrvIndex_t {
    int sel  {-1};
    int sbid {-1};
    int roc  {-1};                           // 0-17
    int feb  {-1};                           // 
    int ch   {-1};                            // 
  };
  
  struct RunData_t {
    int run_number;
    int n_pulsed_channels;
    int pulsed_channel[96];             // only n_pulsed_clannels are used
  };
    
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct CrvcHist_t {
    TH1F* h_dt;
  };
  
  struct ChannelHist_t {
    TH1F* h_ch;
    TH1F* h_dt;
  };

  struct CrvpHist_t {
    TH1F* h_ch;
    TH1F* h_feb;
    TH1F* h_dt;
    TH2F* h_dt_vs_feb;
  };

  struct FebHist_t {
    ChannelHist_t ch[64];
    CrvpHist_t*   crvp;               // for the whole FEB
    TH1F*         h_sbid;
    TH1F*         h_dt;
  };
  
  struct RocHist_t {
    CrvpHist_t*   crvp;               // for the whole ROC
    FebHist_t*    feb [30];
    TH1F*         h_sbid;
    TH2F*          h_feb_vs_dt;
  };
  
  struct Hist_t {
    CrvpHist_t* crvp[10];               // for all
    CrvcHist_t* crvc[10];               // for all
    RocHist_t*  roc [10];               // by ROC
    TH1F*       h_sbid;
    TH2F*       h_feb_vs_ch;
    TH2F*       h_dt_vs_sbid;
    TH2F*       h_feb_vs_sbid[2];       // one per ROC
  };

//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*       fTopFolder;
  TFolder*       fRunFolder;

  Hist_t*        fHist;
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;

  int            fRunNumber;

  int            fMaxEvent;             // for X-axis truncation
  int            fNEvents;
  
  TTree          *fChain;               //! pointer to the analyzed TTree or TChain
  Int_t           fCurrent;             //!current Tree number in a TChain

                                        // will the same directory work ?
  DaqEvent*       fEvent;               // #include "daqana_nt_format.hh"
  TClonesArray*   fSh;
  TClonesArray*   fTc;

  // float           t05 [36];
  // int             n05 [36];
  // float           dt05[36][36];
//-----------------------------------------------------------------------------
                                        // for independent runs, the name should be the same..
                                        // make it different to process the same run with different refence channels
                                        // subrun number - chooses the file, assuming only one file
  
  plot_crv_tc_dt(int RunNumber, int SubrunNumber, const char* Label = "n002");

                                        // one file in an arbitrary place
  
  plot_crv_tc_dt(int RunNumber, const char* Fn);
  
  virtual ~plot_crv_tc_dt();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  int              BookCrvpHistograms (CrvpHist_t*  Hist, CrvIndex_t* Index, TFolder* Folder);
  int              BookCrvcHistograms (CrvcHist_t*  Hist, CrvIndex_t* Index, TFolder* Folder);
  int              BookFebHistograms  (FebHist_t*   Hist, CrvIndex_t* Index, TFolder* Folder);
  int              BookRocHistograms  (RocHist_t*   Hist, CrvIndex_t* Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*      Hist, TFolder* Folder);

  int              FillHistograms       ();

  int              ResetHistograms();
  int              SaveHistograms (const char* Filename);
  
};
#endif
