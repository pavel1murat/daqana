#ifndef __daqana_ana_plot_calo_tc_dt_hh__
#define __daqana_ana_plot_calo_tc_dt_hh__

#include <format>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.
#include "daqana/obj/DaqEvent.hh"

#include "daqana/obj/TrkPanelMap_t.hh"
#include "daqana/obj/CalChannelMap_t.hh"

#include "ana/ana/booking.hh"

class plot_calo_tc_dt: public TNamed {
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
  
  struct CaloIndex_t {
    int sel   {-1};
    int crate {-1};
    int disk  {-1};
    int cid   {-1};                           // channel within the FEB (0-63)
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

  struct DiskHist_t {
    BoardHist_t*  board[30];
    TH1F*         h_ch;
    TH2F*         h_board_vs_dt;
  };
  
  struct Hist_t {
    DiskHist_t* disk[2];                // 2 disks
    TH2F*       h_dt_vs_sipmid;
    TH1F*       h_cid;                  // occupancy offline channel
  };

//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*         fTopFolder;
  TFolder*         fRunFolder;

  Hist_t*          fHist;
  
  Booking*         fBook;

  TrkPanelMap_t*   fTpm;
  
  CalChannelMap_t* fCalCm;

  int            fRunNumber;

  int            fMaxEvent;             // for X-axis truncation
  int            fNEvents;
  
  TTree          *fChain;               //! pointer to the analyzed TTree or TChain
  Int_t           fCurrent;             //!current Tree number in a TChain

                                        // will the same directory work ?
  DaqEvent*       fEvent;               // #include "daqana_nt_format.hh"
  // TClonesArray*   fSh;
  // TClonesArray*   fTc;

  // fit_result_t     fFr[10][25];
  // fit_result_t*    fFrRef;

  // float           t05 [36];
  // int             n05 [36];
  // float           dt05[36][36];
//-----------------------------------------------------------------------------
                                        // for independent runs, the name should be the same..
                                        // make it different to process the same run with different refence channels
                                        // subrun number - chooses the file, assuming only one file
  
  plot_calo_tc_dt(int RunNumber, int SubrunNumber, const char* Label = "n002");

                                        // one file in an arbitrary place
  
  plot_calo_tc_dt(int RunNumber, const char* Fn);
  
  virtual ~plot_calo_tc_dt();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  int              BookCalhHistograms (CalhHist_t*   Hist, CalIndex_t* Index, TFolder* Folder);
  int              BookDiskHistograms (DiskHist_t*   Hist, CalIndex_t* Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*      Hist, TFolder* Folder);

  int              FillDiskHistograms (DiskHist_t* Hist, DaqCalHit* Calh);
  int              FillHistograms     ();

  int              FitHistogram       (TH1* Hist, fit_result_t* Fp, float XMin = 1, float XMax = -1, int NMin = 100);
  
                                        // if TMax > TMin, use them as limits, otherwise determine them automatically
  int              FitTimeOffsets     (float TMin = 1., float TMax = -1.);
  int              PrintTimeCorrections();

  int              ResetHistograms();
  int              SaveHistograms (const char* Filename);
  
};
#endif
