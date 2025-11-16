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
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "daqana/obj/TrkSegment.hh"
#include "daqana/obj/SegmentFit.hh"

class simres {
public:
  struct Straw_t {
    int    fNumber;
    int    fLayer;
    double fX;
    double fY;
    double fR;                            // straw radius
    double fDoca;                         // track to wire distance
    int    fMask;

    Straw_t(int StrawNumber) {
      fMask   = 0;
      fNumber = StrawNumber / 2;             // index in the layer
      fLayer  = StrawNumber % 2;
    }
  };

  struct Panel_t {
    std::vector<Straw_t*>       fListOfStraws;
    int                         fNHits;
    std::vector<Straw_t*>       fListOfHitStraws;
    std::vector<mu2e::ComboHit> fListOfComboHits; // list of faked ComboHits
  };

  struct Hist_t {
    TH1F* fNHits;
  };
  
  Hist_t  fHist;
  Panel_t fPanel;

  int     fDebugMode;
  double  fStrawRadius;
                                        // parameters of the line
  double fX0;
  double fY0;
  double fDyDx;
  double fNx;
  double fNy;
  double fNux;
  double fNuy;

  double fVDrift;                       // drift velocity, mm/ns
  double fT0;                           // T0 of the simulated particle
  TrkSegment  fSegment;                 // input for the fitter
  SegmentFit* fFitter;
//-----------------------------------------------------------------------------
// functions
//------------------------------------------------------------------------------
  simres(int DebugMode = 0);

                                        // convert drift distance into the drift time
  double D2T(double R);
                                        // convert [reconstructed] into hit radius used by the fit
  double T2D(double T);

  Straw_t* Straw   (int I) { return fPanel.fListOfStraws   [I]; }
  Straw_t* HitStraw(int I) { return fPanel.fListOfHitStraws[I]; }

  int init_geometry  ();
  int book_histograms();
  int fill_histograms(Hist_t* Hist);
  int run            (int NEvents = 1000, int DebugMode = 0);
  int display_event  ();
};

#endif
