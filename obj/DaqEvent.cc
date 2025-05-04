//
#include <vector>

#include "daqana/obj/DaqEvent.hh"
// #include "DaqEvent.hh"

// ClassImp(DaqEvent)

//-----------------------------------------------------------------------------
DaqEvent::DaqEvent() { // : TObject () {
  nsd = 0;
  sd = new TClonesArray("DaqStrawDigi",100);
  nsh = 0;
  sh = new TClonesArray("DaqStrawHit" ,100);
}

//-----------------------------------------------------------------------------
DaqEvent::~DaqEvent() {
  sd->Delete();
  delete sd;
  sh->Delete();
  delete sh;
}
