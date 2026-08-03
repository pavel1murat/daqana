#ifndef __daqana_DaqEvent_hh__
#define __daqana_DaqEvent_hh__

#include "TClonesArray.h"
#include "TObject.h"

#include "daqana/obj/DaqCaloDigi.hh"
#include "daqana/obj/DaqCaloHit.hh"
#include "daqana/obj/DaqCaloCluster.hh"

#include "daqana/obj/DaqCrvDigi.hh"
#include "daqana/obj/DaqCrvRecoPulse.hh"
#include "daqana/obj/DaqCrvCoincidenceCluster.hh"

#include "daqana/obj/DaqStrawDigi.hh"
#include "daqana/obj/DaqStrawHit.hh"
#include "daqana/obj/DaqComboHit.hh"
#include "daqana/obj/DaqFragment.hh"
#include "daqana/obj/DaqSegment.hh"
#include "daqana/obj/DaqTimeCluster.hh"
#include "daqana/obj/DaqTrack.hh"
#include "daqana/obj/DaqTrkStrawHit.hh"

class DaqEvent { // : public TObject {
public:
  int            run;
  int            srn;                   // subrun number
  int            evn;                   // event number within the subrun
  long           ewtag;                 // event window tag
  int            nsdtot;                // number of straw digis in event
  TClonesArray*  sd;
  int            nshtot;                // total number of straw digis in event
  int            nsh[36][6];            // per panel, by [dtc_id-1][link] , dtc_id, NOT pcie_arrd
  int            pmp[36];               // per plane, by [dtc_id-1]
  
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

  int            ncald;                 // number of calo digis
  TClonesArray*  cald;

  int            ncalh;                 // number of calo hits (normally, 2 digis per hit)
  TClonesArray*  calh;

  int            ncalc;                 // number of calo clusters
  TClonesArray*  calc;

  float          edisk[2];
  float          ecal;

  int            ncrvd;                 // number of crv digis
  TClonesArray*  crvd;

  int            ncrvp;                 // number of crv digis
  TClonesArray*  crvp;

  int            ncrvc;                 // number of crv digis
  TClonesArray*  crvc;

  int            nfrag;                 // number of ARTDAQ fragments (so far, only TRK)
  TClonesArray*  frag;

  // int            nstmdigis;          // number of stm digis
  // TClonesArray*  stmdigis;

  DaqStrawDigi*    Sd(int I) { return (DaqStrawDigi*  ) sd->At(I); }
  DaqStrawHit*     Sh(int I) { return (DaqStrawHit*   ) sh->At(I); }
  DaqTimeCluster*  Tc(int I) { return (DaqTimeCluster*) tc->At(I); }
  DaqSegment*      Seg(int I) { return (DaqSegment*) seg->At(I); }
  int              Nsh(int Plane, int Panel) { return nsh[Plane][Panel]; }
  int              Pmp(int DtcID) { return pmp[DtcID] ; }

  DaqStrawDigi* NewSD(int I) {
    DaqStrawDigi* digi = new ((*sd)[I]) DaqStrawDigi();
    return digi;
  }
  
  DaqFragment* NewFragment(int I) {
    DaqFragment* df = new ((*frag)[I]) DaqFragment();
    return df;
  }
  
  /* virtual */ void     Clear(const char* Opt = "") ; //  override ;
  
  DaqEvent();
  ~DaqEvent();

  //  ClassDefOverride(DaqEvent,1)
};

#endif
