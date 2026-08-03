//

#include "daqana/obj/DaqCaloCluster.hh"

ClassImp(DaqCaloCluster)

//-----------------------------------------------------------------------------
DaqCaloCluster::DaqCaloCluster() : TObject() {
  disk     =  0;
  size     = -1;
  split    =  0;
  time     = -1.;
  sigt     = -1;
  edep     = -1.;
  sige     = -1;
  x        = -9999;
  y        = -9999;
  z        = -9999;
}

//-----------------------------------------------------------------------------
DaqCaloCluster::~DaqCaloCluster() {
}

//-----------------------------------------------------------------------------
// does nothing
//-----------------------------------------------------------------------------
int DaqCaloCluster::Init(int Nd) {
  return 0;
}

//-----------------------------------------------------------------------------
void DaqCaloCluster::Clear(const char* Opt) {
}

