//
#ifndef __test3_hh__
#define __test3_hh__

#include <iostream>
#include <format>
#include "TH2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/obj/TrkSegment.hh"

#include "daqana/obj/SegmentFit.hh"
//-----------------------------------------------------------------------------
// global variables
//-----------------------------------------------------------------------------
class test3 {
public:
  TCanvas*    fCanvas;
  TH2F*       fH2;
  TrkSegment* fSegment;
  SegmentFit* fSfitter;

  test3() {
    fSegment = nullptr;
    fSfitter = nullptr;
    fCanvas  = nullptr;
    fH2      = nullptr;
  }

  int init_segment_geometry(TrkSegment* Seg);
  int test_fit_line        (const char* Fn, int Plane, int Panel, int NIter=0);
};

#endif
