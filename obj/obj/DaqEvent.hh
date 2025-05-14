#ifndef __daqana_DaqEvent_hh__
#define __daqana_DaqEvent_hh__

#include "TClonesArray.h"
#include "TObject.h"

#include "daqana/obj/DaqStrawDigi.hh"
#include "daqana/obj/DaqStrawHit.hh"
#include "daqana/obj/DaqTimeCluster.hh"

class DaqEvent { // : public TObject {
public:
  int            run;
  int            srn;                   // subrun number
  int            evn;                   // event number
  int            nsd;                   // number of straw digis in event
  TClonesArray*  sd;
  int            nshtot;                // total number of straw digis in event
  int            nsh[2][6];             // by [dtc][link]
  TClonesArray*  sh;
  int            ntc;                   // N(time clusters)
  TClonesArray*  tc;
  // int            ncalodigis;         // number of calo digis
  // TClonesArray*  calodigis;
  // int            ncrvdigis;          // number of crv digis
  // TClonesArray*  crvdigis;
  // int            nstmdigis;          // number of stm digis
  // TClonesArray*  stmdigis;

  //   static std::vector<DaqStrawDigi> fgSd;

  DaqStrawDigi*    Sd(int I) { return (DaqStrawDigi*  ) sd->At(I); }
  DaqStrawHit*     Sh(int I) { return (DaqStrawHit*   ) sh->At(I); }
  DaqTimeCluster*  Tc(int I) { return (DaqTimeCluster*) tc->At(I); }
  
  DaqEvent();
  ~DaqEvent();

  // ClassDef(DaqEvent,1)
};

#endif
