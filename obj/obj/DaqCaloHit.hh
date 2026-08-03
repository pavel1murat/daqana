#ifndef __daqana_obj_DaqCaloHit_hh__
#define __daqana_obj_DaqCaloHit_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCaloHit : public TObject {
public:
  short  cid;                           // crystal ID
  short  nd;                            // N(digis)
  short  nsipms;
  float  time;
  float  sigt;
  float  edep;                          // peak position
  float  sige;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCaloHit();
  virtual ~DaqCaloHit();

  int     Init(int Ns);

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCaloHit,1);
};

#endif
