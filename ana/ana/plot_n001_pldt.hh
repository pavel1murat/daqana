#ifndef __daqana_ana_plot_n001_pldt_hh__
#define __daqana_ana_plot_n001_pldt_hh__

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

class plot_n001_pldt: public TNamed {
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
  struct ChannelHist_t {
    TH1F*     h_tdc0;
    TH1F*     h_dt10;
    TH1F*     h_ph;                       // pulse height
    TH1F*     h_bl;                       // pulse height
    TH1F*     h_fs;                       // pulse height
    TH1F*     h_edep;                     // 
  };
  
  struct PanelHist_t {
    ChannelHist_t* ch[96];
    TH1F*          h_occup;                    // pulse height
    TH1F*          h_edep;                     // 
  } ;
  
  struct SlotHist_t {
    PanelHist_t* panel[12];
    TH1F*        h_edep;                        // 
  };
  
  struct Hist_t {
    TH1F*       h_occup;                       // occupancy per panel, offline indexing
    TH2F*       h_occup_2d;                    // occupancy
    TH1F*       h_dt05[36];
    SlotHist_t* slot[18];
  };

//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*       fTopFolder;
  TFolder*       fRunFolder;

  Hist_t*        fHist[10];
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;

  int            fRunNumber;

  int            fMaxEvent;             // for X-axis truncation
  int            fNEvents;
  
  DaqStrawDigi*  fSdr[216];
  
  TTree          *fChain;               //! pointer to the analyzed TTree or TChain
  Int_t           fCurrent;             //!current Tree number in a TChain

                                        // will the same directory work ?
  DaqEvent*       fEvent;               // #include "daqana_nt_format.hh"
  TClonesArray*   fSd;

  float           t05 [36];
  float           dtpl[36];
  int             n05 [36];
//-----------------------------------------------------------------------------
                                        // for independent runs, the name should eb the same..
                                        // make it different to process the same run with different refence channels
  plot_n001_pldt(int RunNumber, const char* Fn = nullptr);
  
  virtual ~plot_n001_pldt();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  int              BookChannelHistograms(ChannelHist_t* Hist, Index_t* Index, TFolder* Folder);
  int              BookPanelHistograms  (PanelHist_t*   Hist, Index_t* Index, TFolder* Folder);
  int              BookSlotHistograms   (SlotHist_t*    Hist, Index_t* Index, TFolder* Folder);
  int              BookSelHistograms    (Hist_t*        Hist, Index_t* Index, TFolder* Folder);
  int              BookHistograms       (TFolder* Folder);

  int              FillChannelHistograms(ChannelHist_t* Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh);
  int              FillPanelHistograms  (PanelHist_t*   Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh);
  int              FillSlotHistograms   (SlotHist_t*    Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh);
  int              FillHistograms       ();

  int              ResetHistograms();
  int              SaveHistograms (const char* Filename);
  
  int              PrintHistograms   (int ISet);
  int              PrintNoisyChannels(float Percentage = 0.01);
};
#endif
