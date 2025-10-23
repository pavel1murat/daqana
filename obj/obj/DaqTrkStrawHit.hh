///////////////////////////////////////////////////////////////////////////////
// R - position
// r - drift radius (distance)
// iseg and itrk should default to -1
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqTrkStrawHit_hh__
#define __daqana_obj_DaqTrkStrawHit_hh__

#include "TObject.h"

class DaqTrkStrawHit : public TObject {
public:
  int     sid;                          // hit id = sid | (segment #) << 16 | (track #) << 24
  int     zface;                        // z-ordered face ... can be deduced from sid...
  int     mnid;                         // Minnesota panel ID 
  float   time;                         // 0:CAL
  float   dt;
  float   tot0;
  float   tot1;
  float   edep;
  int     iseg;
  int     itrk;
  float   rdrift;                       // drift distance
  float   doca;                         // vdot((track-to-hit dist),(track-to-wire distance).unit)
  float   drho;                         // unsigned: |track-to-wire dist| - rdrift

  DaqTrkStrawHit();
  
  virtual ~DaqTrkStrawHit();

  ClassDef(DaqTrkStrawHit,1);
};

#endif
