//
#include <vector>

#include "daqana/obj/DaqEvent.hh"

// ClassImp(DaqEvent)

//-----------------------------------------------------------------------------
DaqEvent::DaqEvent() { // : TObject() {
  sd    = new TClonesArray("DaqStrawDigi"  ,100);
  sh    = new TClonesArray("DaqStrawHit"   ,100);
  ch    = new TClonesArray("DaqComboHit"   ,100);
  tc    = new TClonesArray("DaqTimeCluster", 10);
  trk   = new TClonesArray("DaqTrack", 10);
  seg   = new TClonesArray("DaqSegment"    ,  10);
  segsh = new TClonesArray("DaqTrkStrawHit", 100);
  trksh = new TClonesArray("DaqTrkStrawHit", 100);

  Clear();
}

//-----------------------------------------------------------------------------
DaqEvent::~DaqEvent() {
  sd->Delete(); delete sd;
  sh->Delete(); delete sh;
  ch->Delete(); delete ch;
  tc->Delete(); delete tc;
  trk->Delete(); delete trk;
  seg->Delete(); delete seg;
  segsh->Delete(); delete segsh;
  trksh->Delete(); delete trksh;
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

  nsdtot  = 0; sd->Clear();
  nshtot  = 0; sh->Clear();
  nch     = 0; ch->Clear();
  ntc     = 0; tc->Clear();
  ntrk    = 0; trk->Clear();
  ntrksh  = 0; trksh->Clear();
  nseg    = 0; seg->Clear();
  nsegsh  = 0; segsh->Clear();

}
