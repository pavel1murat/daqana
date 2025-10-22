///////////////////////////////////////////////////////////////////////////////
// R - position
// r - drift radius (distance)
// iseg and itrk should default to -1
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqTrkStrawHit_hh__
#define __daqana_obj_DaqTrkStrawHit_hh__

#include "daqana/obj/DaqStrawHit.hh"

class DaqTrkStrawHit : public DaqStrawHit {
public:
  float   rdrift;                       // drift distance
  float   doca;                         // vdot((track-to-hit dist),(track-to-wire distance).unit)
  float   drho;                         // unsigned: |track-to-wire dist| - rdrift
  int     iseg;
  int     itrk;

  DaqTrkStrawHit();
  
  virtual ~DaqTrkStrawHit();

  ClassDef(DaqTrkStrawHit,1);
};

#endif
