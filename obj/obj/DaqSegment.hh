///////////////////////////////////////////////////////////////////////////////
// for now, assume the time cluster in one station, update later
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqSegment_hh__
#define __daqana_obj_DaqSegment_hh__

#include "TObject.h"

class DaqSegment : public TObject {
public:
  int     nhits ;
  float   t0    ;
  float   chi2  ;
  
  DaqSegment();

  virtual ~DaqSegment();

  ClassDef(DaqSegment,1);
};

#endif
