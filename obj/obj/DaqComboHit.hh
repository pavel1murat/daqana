#ifndef __daqana_obj_DaqComboHit_hh__
#define __daqana_obj_DaqComboHit_hh__

#include "TClonesArray.h"
#include "TObject.h"

class DaqComboHit : public TObject {
public:
  int     sid;
  int     nsh;
  int     zface;                        // z-ordered face
  int     mnid;                         // Minnesota panel ID 
  float   time;                         // correctedTime
  float   dtime;                        // drift time
  float   x;
  float   y;
  float   z;
  float   ux;                           // nz is always zero
  float   uy;
  float   ures;
  float   vres;
  float   edep;

  DaqComboHit();
                                        // low 7 bits
  int straw  () { return (sid >>  0) & 0x007f ; }
  int panel  () { return (sid >>  7) & 0x0007 ; }
  int face   () { return (sid >>  7) & 0x0001 ; }
  int plane  () { return (sid >> 10) & 0x003f ; }
  int station() { return (sid >> 11) & 0x001f ; }
  
  virtual ~DaqComboHit();

  ClassDef(DaqComboHit,1);
};

#endif
