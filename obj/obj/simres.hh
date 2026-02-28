///////////////////////////////////////////////////////////////////////////////
// DebugMode:
//   bit0: printouts in the start and the end of each function
///////////////////////////////////////////////////////////////////////////////
#ifndef __simres_hh__
#define __simres_hh__

#include <iostream>
#include <format>
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEllipse.h"
#include "TRandom3.h"
#include "TFolder.h"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "daqana/obj/TrkSegment.hh"
#include "daqana/obj/SegmentFit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"

class simres {
public:
  struct Straw_t {
    mu2e::StrawId  id;                          // 
    double         fX;
    double         fY;
    double         fR;                            // straw radius
    double         fDoca;                         // track to wire distance
    int            fMask;

    Straw_t(mu2e::StrawId Id) {
      id      = Id;
      fMask   = 0;
    }
  };

  enum {
    kNHitHistSets     = 1000,           // used up to 109
    kNSegmentHistSets =   10,
    kNEventHistSets   =   10,
  };

  struct Panel_t {
    std::vector<Straw_t*>        fListOfStraws;
    std::vector<Straw_t*>        fListOfHitStraws;
    std::vector<mu2e::ComboHit*> fListOfComboHits; // list of faked ComboHits

    double fStrawOffsetX[96];           // provision for simulating misalignments
    double fStrawOffsetY[96];

    void Clear();

    int NHits() { return fListOfHitStraws.size(); }
  };

  struct HitHist_t {
    TH1F* r;                           // drift radius
    TH1F* doca;                        // [signed] track-to-wire distance
    TH1F* dr;
    TH1F* drho;
    TH2F* dr_vs_r;
    TH2F* drho_vs_r;
    TH2F* dr_vs_is;
    TH2F* drho_vs_is;
  };

  struct EventHist_t {
    TH1F* nhits;
  };

struct SegmentHist_t {
    TH1F* dy;
    TH1F* ddy;
    TH1F* dt0;
  };

  struct Hist_t {
    EventHist_t*   fEvt[kNEventHistSets  ];
    SegmentHist_t* fSeg[kNSegmentHistSets];
    HitHist_t*     fHit[kNHitHistSets    ];
  };

  struct HitData_t {
    int    index;                       // hit index within the segment
    int    is;                          // straw number
    double r;                           // drift radius
    double dr;                          // residual
    double drho;                        // rdrift - track_doca
    double doca;                        // 
  } fHitData;

  struct EventData_t {
    int    number;
    int    nhits;
    double dy;
    double ddy;
    double dt0;
  } fEvent;

  Hist_t  fHist;
  Panel_t fPanel;

  int     fDebugMode;
  double  fStrawRadius;

  double  fSigmaT0;                     // timing resolution
  int     fBeamMode;
                                        // parameters of the line
  double fX0;
  double fY0;
  double fDyDx;
  double fNx;
  double fNy;
  double fNux;
  double fNuy;

  double fMaxGoodDy;

  double      fVDrift;                       // drift velocity, mm/ns
  double      fT0;                           // T0 of the simulated particle
  TrkSegment  fSegment;                 // input for the fitter
  SegmentFit* fFitter;

  TRandom3    fTrn3;

  TFolder*    fFolder;
//-----------------------------------------------------------------------------
// functions
//------------------------------------------------------------------------------
  simres(int DebugMode = 0);

                                        // convert drift distance into the drift time
  double D2T(double R);
                                        // convert [reconstructed] into hit radius used by the fit
  double T2D(double T);

  Straw_t*        Straw   (int I) { return fPanel.fListOfStraws   [I]; }
  Straw_t*        HitStraw(int I) { return fPanel.fListOfHitStraws[I]; }
  mu2e::ComboHit* ComboHit(int I) { return fPanel.fListOfComboHits[I]; }

  int init_geometry  ();

  int book_event_histograms  (EventHist_t*   Hist, int ISet, const char* Folder);
  int book_hit_histograms    (HitHist_t*     Hist, int ISet, const char* Folder);
  int book_segment_histograms(SegmentHist_t* Hist, int ISet, const char* Folder);
  int book_histograms        ();

  int fill_event_histograms  (EventHist_t*   Hist);
  int fill_hit_histograms    (HitHist_t*     Hist, HitData_t* Hd);
  int fill_segment_histograms(SegmentHist_t* Hist);
  int fill_histograms        (Hist_t*        Hist);

  int simulate_line_parameters();
  int simulate_hits();
  int process_event(int IEvent);
  int run(int NEvents = 1000, int DebugMode = 0);

  int display_event  ();

  void  print_(const std::string&  Message,
               const std::source_location& location = std::source_location::current());

//-----------------------------------------------------------------------------
// the following helper methods allow to save 1 line per request, which in 
// case of 100's histograms booked is a non-negligible number
//-----------------------------------------------------------------------------
  void  DeleteHistograms(TFolder* Folder = (TFolder*) -1);

  void  AddHistogram(TObject* hist, const char* FolderName = "Hist");

  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, const float* LowEdge, 
		const char* FolderName = "Hist");

  void  HBook1D(TH1D*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook1D(TH1D*& Hist, const char* Name, const char* Title,
		Int_t Nx, const double* LowEdge, 
		const char* FolderName = "Hist");

  void  HBook2F(TH2F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		Int_t Ny, Double_t YMin, Double_t YMax,
		const char* FolderName = "Hist");

  void  HProf (TProfile*& Hist, const char* Name, const char* Title,
	       Int_t Nx, Double_t XMin, Double_t XMax,
	       Double_t YMin, Double_t YMax,
	       const char* FolderName = "Hist");

  static int  SaveFolder(TFolder* Folder, TDirectory* Dir);
  void        SaveHist  (const char* Filename);
};

#endif
