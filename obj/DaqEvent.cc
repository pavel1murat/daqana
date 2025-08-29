//
#include <vector>

#include "daqana/obj/DaqEvent.hh"

// ClassImp(DaqEvent)

//-----------------------------------------------------------------------------
DaqEvent::DaqEvent() { // : TObject() {
  nsdtot  = 0; sd  = new TClonesArray("DaqStrawDigi"  ,100);
  nshtot  = 0; sh  = new TClonesArray("DaqStrawHit"   ,100);
  nch     = 0; ch  = new TClonesArray("DaqComboHit"   ,100);
  ntc     = 0; tc  = new TClonesArray("DaqTimeCluster", 10);
  ntracks = 0; trk = new TClonesArray("DaqTrack", 10);

  Clear();
}

//-----------------------------------------------------------------------------
DaqEvent::~DaqEvent() {
  sd->Delete(); delete sd;
  sh->Delete(); delete sh;
  ch->Delete(); delete ch;
  tc->Delete(); delete tc;
  trk->Delete(); delete trk;
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

  nsdtot  = 0;
  nshtot  = 0;
  nch     = 0;
  ntc     = 0;
  ntracks = 0;
  sd->Clear();
  sh->Clear();
  ch->Clear();
  tc->Clear();
  trk->Clear();
}
