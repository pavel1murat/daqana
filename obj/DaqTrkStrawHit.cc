//

#include "daqana/obj/DaqTrkStrawHit.hh"

ClassImp(DaqTrkStrawHit)

//-----------------------------------------------------------------------------
DaqTrkStrawHit::DaqTrkStrawHit() : DaqStrawHit() {
                                        // 
  rdrift =  1.e6;                       // 
  doca   =  1.e6;                       // vdot(dist.track-wire,dist.track-hit.unit)
  drho   =  1.e6;                       // residual sign: (Rtrk-Rhit)*vdrift (unit vectors)

  iseg   = -1;                          // 
  itrk   = -1;
}

//-----------------------------------------------------------------------------
DaqTrkStrawHit::~DaqTrkStrawHit() {
}
