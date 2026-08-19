#ifndef __daqana_ana_plot_n006_time_res_hh__
#define __daqana_ana_plot_n006_time_res_hh__

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
#include "daqana/obj/RunDb.hh"

#include "ana/ana/booking.hh"

class plot_n006_time_res: public TNamed {
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
  
  RunDb::Data_t     fRunInfo;
    
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  enum {
    kMaxSelHistSets = 100,
  };
  
  struct SelHist_t {
    TH1F*        h_ecal;                        // 
    TH1F*        h_edisk[2];                     // 
  };

  struct Hist_t {
    TH1F*        h_dt_calh_tc;                        // 
    TH1F*        h_dt_crvc_tc[2];                     //
    TH1F*        h_n_close_calh;
    TH1F*        h_n_close_crvc;
    
    SelHist_t*   sel[kMaxSelHistSets];
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
  TClonesArray*   fSd;

  int             fNCloseCalh;
  int             fNCloseCrvc;
//-----------------------------------------------------------------------------
                                        // for independent runs, the name should be the same..
                                        // make it different to process the same run with different refence channels
  
  plot_n006_time_res(int RunNumber, const char* Fn = nullptr, const char* Label="n006");
  
  virtual ~plot_n006_time_res();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init     (TTree   *tree);

  void             Loop          (int NEvents = -1);
  
                                        // returns 1 if event loaded, 0 otherwise
  int              LoadEvent     (int Run, int Subrun, int Event);

  int              BookSelHistograms    (SelHist_t* Hist, Index_t* Index, TFolder* Folder);
  int              BookHistograms       (TFolder* Folder);

  int              FillSelHistograms    (SelHist_t* Hist);
  
  int              FillHistograms       ();

  int              ResetHistograms();
  int              SaveHistograms (const char* Filename);
  
  int              PrintHistograms   (int ISet);
};
#endif
