#ifndef __daqana_obj_DaqStrawHit_hh__
#define __daqana_obj_DaqStrawHit_hh__

#include "TClonesArray.h"
#include "TObject.h"

class DaqStrawHit : public TObject {
public:
  int     sid;
  int     zface;                           // z-ordered face
  int     mnid;                            // Minnesota panel ID 
  float   time;                            // 0:CAL
  float   dt;
  float   tot0;
  float   tot1;
  float   edep;

  DaqStrawHit();
                                        // low 7 bits
  int Straw  () { return (sid >>  0) & 0x007f ; }
  int Panel  () { return (sid >>  7) & 0x0007 ; }
  int Face   () { return (sid >>  7) & 0x0001 ; }
  int Plane  () { return (sid >> 10) & 0x003f ; }
  int Station() { return (sid >> 11) & 0x001f ; }
  
  virtual ~DaqStrawHit();

  ClassDef(DaqStrawHit,1);
};

#endif
