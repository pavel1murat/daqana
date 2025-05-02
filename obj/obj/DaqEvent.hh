#ifndef __daqana_DaqEvent_hh__
#define __daqana_DaqEvent_hh__

#include "TClonesArray.h"
#include "TObject.h"

#include "otsdaq-mu2e-tracker/Nt/DaqStrawDigi.hh"

class DaqEvent { // : public TObject {
public:
  int            run;
  int            srn;                   // subrun number
  int            evn;                   // event number
  int            nsd;         // number of straw digis in event
  TClonesArray*  sd;
  // int            ncalodigis;          // number of calo digis
  // TClonesArray*  calodigis;
  // int            ncrvdigis;           // number of crv digis
  // TClonesArray*  crvdigis;
  // int            nstmdigis;           // number of stm digis
  // TClonesArray*  stmdigis;

  //   static std::vector<DaqStrawDigi> fgSd;

  DaqStrawDigi* sdigi(int I) { return (DaqStrawDigi*) sd->At(I); }
  
  DaqEvent();
  ~DaqEvent();

  // ClassDef(DaqEvent,1)
};

#endif
