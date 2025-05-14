//
#include <vector>

#include "daqana/obj/DaqEvent.hh"
// #include "DaqEvent.hh"

// ClassImp(DaqEvent)

//-----------------------------------------------------------------------------
DaqEvent::DaqEvent() { // : TObject () {
  nsd    = 0;
  sd     = new TClonesArray("DaqStrawDigi",100);
  nshtot = 0;
  sh     = new TClonesArray("DaqStrawHit" ,100);
  ntc    = 0;
  tc     = new TClonesArray("DaqTimeCluster",10);
}

//-----------------------------------------------------------------------------
DaqEvent::~DaqEvent() {
  sd->Delete();
  delete sd;
  sh->Delete();
  delete sh;
  tc->Delete();
  delete tc;
}
