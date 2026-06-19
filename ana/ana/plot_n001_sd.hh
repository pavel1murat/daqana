#ifndef __daqana_ana_plot_n001_sd_hh__
#define __daqana_ana_plot_n001_sd_hh__

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

class plot_n001_sd: public TNamed {
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
  struct RunData_t {
    int run_number;
    int n_pulsed_channels;
    int pulsed_channel[96];             // only n_pulsed_clannels are used
  };
    
  struct PanelData_t {
    TrkPanelMap_t::Data_t* tpm;
                                        // SD index in the event
    int i21;                            // index of the first hit in channel 21 (assume pulsing channels 5+8*i...)
  };
  
  struct PlaneData_t {
    PanelData_t panel[6];
  };
    
  struct Data_t {
    PlaneData_t plane[36];
  };

//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct Hist_t {
    TH1F*     h_tdc0;
    TH1F*     h_ch;                       // straw
    TH1F*     h_ph;                       // pulse height
    TH1F*     h_plane;
    TH2F*     h_panel_dt;
    TH2F*     h_panel_dt_111;             // wrt first hit in panel 111
    TProfile* h_panel_dt_111_vs_evn[216]; // wrt first hit in panel 111
    TH1F*     h_dt20[6];
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
  
  plot_n001_sd(const char* Name, int RunNumber, const char* Fn = nullptr);
  
  virtual ~plot_n001_sd();
  
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
