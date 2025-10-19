///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqStrawHit_hh__
#define __daqana_obj_DaqStrawHit_hh__

#include "TObject.h"

class DaqStrawHit : public TObject {
public:
  int     sid;                          // hit id = sid | (segment #) << 16 | (track #) << 24
  int     zface;                        // z-ordered face ... can be deduced from sid...
  int     mnid;                         // Minnesota panel ID 
  float   time;                         // 0:CAL
  float   dt;
  float   tot0;
  float   tot1;
  float   edep;

  DaqStrawHit();

  int straw  () { return (sid >>  0) & 0x007f ; }
  int panel  () { return (sid >>  7) & 0x0007 ; }
  int face   () { return (sid >>  7) & 0x0001 ; } // not quite right
  int plane  () { return (sid >> 10) & 0x003f ; }
  int station() { return (sid >> 11) & 0x001f ; }

  virtual ~DaqStrawHit();

  ClassDef(DaqStrawHit,1);
};

#endif
