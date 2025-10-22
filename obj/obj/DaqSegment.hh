///////////////////////////////////////////////////////////////////////////////
// for now, assume the time cluster in one station, update later
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqSegment_hh__
#define __daqana_obj_DaqSegment_hh__

#include "TObject.h"

class DaqSegment : public TObject {
public:
  int     sid   ;                       // straw ID of the straw=0 of the panel
  int     nh    ;
  int     ngh   ;
  int     nghl[2];                      // N(Good Hits/Layer)
  int     nmhl[2];                      // N(Missing Hits/Layer)
  float   t0    ;
  float   chi2d ;                       // chi2/dof
  float   y0    ;
  float   z0    ;                       // Z0 - to calculate the track parameters
  float   ymean ;
  float   dzdy  ;
  float   y0t   ;                       // at z=Z(mid panel)
  float   dzdyt ;                       // dzdy of the track (local coord system)
  
  DaqSegment();

  virtual ~DaqSegment();

  int panel  () { return (sid >>  7) & 0x0007 ; }
  int plane  () { return (sid >> 10) & 0x003f ; }

  ClassDef(DaqSegment,1);
};

#endif
