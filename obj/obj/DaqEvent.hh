#ifndef __daqana_DaqEvent_hh__
#define __daqana_DaqEvent_hh__

#include "TClonesArray.h"
#include "TObject.h"

#include "daqana/obj/DaqStrawDigi.hh"
#include "daqana/obj/DaqStrawHit.hh"
#include "daqana/obj/DaqComboHit.hh"
#include "daqana/obj/DaqTimeCluster.hh"
#include "daqana/obj/DaqTrack.hh"
#include "daqana/obj/DaqSegment.hh"
#include "daqana/obj/DaqTrkStrawHit.hh"

class DaqEvent { // : public TObject {
public:
  int            run;
  int            srn;                   // subrun number
  int            evn;                   // event number
  int            nsdtot;                // number of straw digis in event
  short          nsd[2][6][96];         // by [dtc][link][channel]
  TClonesArray*  sd;
  int            nshtot;                // total number of straw digis in event
  int            nsh[2][6];             // by [dtc][link]
  TClonesArray*  sh;
  float          maxEdep;               // max hit Edep in the event
  int            nch;
  TClonesArray*  ch;
  int            ntc;                   // N(time clusters)
  TClonesArray*  tc;
  int            ntrk;                  // N(tracks)
  TClonesArray*  trk;
  int            ntrksh;                // N(segment straw hits)
  TClonesArray*  trksh;                 // straw hits associated with tracks
  int            nseg;                  // N(standalone segments)
  TClonesArray*  seg;                   // 
  int            nsegsh;                // N(segment straw hits)
  TClonesArray*  segsh;                 // straw hits associated with segments

  // int            ncalodigis;         // number of calo digis
  // TClonesArray*  calodigis;
  // int            ncrvdigis;          // number of crv digis
  // TClonesArray*  crvdigis;
  // int            nstmdigis;          // number of stm digis
  // TClonesArray*  stmdigis;

  DaqStrawDigi*    Sd(int I) { return (DaqStrawDigi*  ) sd->At(I); }
  DaqStrawHit*     Sh(int I) { return (DaqStrawHit*   ) sh->At(I); }
  DaqTimeCluster*  Tc(int I) { return (DaqTimeCluster*) tc->At(I); }
  int              Nsh(int I, int J) { return nsh[I][J]; }
  int              Nsd(int Dtc, int Link, int Channel) { return nsd[Dtc][Link][Channel]; }

  /* virtual */ void     Clear(const char* Opt = "") ; // override;
  
  DaqEvent();
  ~DaqEvent();

  //   ClassDef(DaqEvent,1)
};

#endif
