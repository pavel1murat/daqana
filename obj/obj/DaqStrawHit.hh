#ifndef __daqana_obj_DaqStrawHit_hh__
#define __daqana_obj_DaqStrawHit_hh__

#include "TClonesArray.h"
#include "TObject.h"

class DaqStrawHit : public TObject {
public:
  int     sid;
  float   time;                            // 0:CAL
  float   dt;
  float   tot0;
  float   tot1;
  float   edep;

  DaqStrawHit();

  void     Init(int Ns);
  virtual ~DaqStrawHit();

  ClassDef(DaqStrawHit,1);
};

#endif
