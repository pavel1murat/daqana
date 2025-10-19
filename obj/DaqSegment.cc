//

#include "daqana/obj/DaqSegment.hh"

ClassImp(DaqSegment)

//-----------------------------------------------------------------------------
DaqSegment::DaqSegment() : TObject() {
  nh    = -1;
  for (int i=0; i<2; ++i) {
    nghl[i] = -1;
    nmhl[i] = -1;
  }

  chi2d = -1.;
  t0    = 1.e6;
  y0    = 1.e6;
  ymean = 1.e6;
  dzdy  = 1.e6;
  y0t   = 1.e6;
  dzdyt = 1.e6;
}

//-----------------------------------------------------------------------------
DaqSegment::~DaqSegment() {
}
