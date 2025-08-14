//

#include "daqana/obj/DaqComboHit.hh"

ClassImp(DaqComboHit)

//-----------------------------------------------------------------------------
DaqComboHit::DaqComboHit() : TObject(),
  sid(-1),nsh(-1),time(-1.),dtime(1.e6),x(1.e6),y(1.e6),z(1.e6),
  ux(1.e6),uy(1.e6),ures(1.e6),vres(1.e6),edep(-1.) {
}

//-----------------------------------------------------------------------------
DaqComboHit::~DaqComboHit() {
}
