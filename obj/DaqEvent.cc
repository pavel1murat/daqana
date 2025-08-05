//
#include <vector>

#include "daqana/obj/DaqEvent.hh"

// ClassImp(DaqEvent)

//-----------------------------------------------------------------------------
DaqEvent::DaqEvent() { // : TObject() {
  nsdtot = 0; sd = new TClonesArray("DaqStrawDigi"  ,100);
  nshtot = 0; sh = new TClonesArray("DaqStrawHit"   ,100);
  ntc    = 0; tc = new TClonesArray("DaqTimeCluster", 10);

  Clear();
}

//-----------------------------------------------------------------------------
DaqEvent::~DaqEvent() {
  sd->Delete(); delete sd;
  sh->Delete(); delete sh;
  tc->Delete(); delete tc;
}


//-----------------------------------------------------------------------------
void DaqEvent::Clear(const char* Opt) {
  for (int i=0; i<2; i++) {
    for (int link=0; link<6; link++) {
      nsh[i][link] = 0;
    }
  }

  for (int i=0; i<2; i++) {
    for (int link=0; link<6; link++) {
      for (int ich=0; ich<96; ich++) {
        nsd[i][link][ich] = 0;
      }
    }
  }

  sd->Clear();
  sh->Clear();
  tc->Clear();
}
