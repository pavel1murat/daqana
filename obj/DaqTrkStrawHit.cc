//

#include "daqana/obj/DaqTrkStrawHit.hh"

ClassImp(DaqTrkStrawHit)

//-----------------------------------------------------------------------------
DaqTrkStrawHit::DaqTrkStrawHit() : TObject() {

  sid  = -1;
  time = -1;
  dt   = -1;
  tot0 = -1.;
  tot1 = -1.;
  edep = -1;
                                        // 
  rdrift =  1.e6;                       //
  dr     =  1.e6;
  doca   =  1.e6;                       // vdot(dist.track-wire,dist.track-hit.unit)
  drho   =  1.e6;                       // residual sign: (Rtrk-Rhit)*vdrift (unit vectors)

  iseg   = -1;                          // 
  itrk   = -1;
  ihit   = -1;
}

//-----------------------------------------------------------------------------
DaqTrkStrawHit::~DaqTrkStrawHit() {
}
