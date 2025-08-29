///////////////////////////////////////////////////////////////////////////////
// for now, assume the time cluster in one station, update later
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqTrack_hh__
#define __daqana_obj_DaqTrack_hh__

#include "TObject.h"

class DaqTrack : public TObject {
public:
  int     nhits ;
  float   t0    ;
  float   chi2  ;
  
  DaqTrack();

  virtual ~DaqTrack();

  ClassDef(DaqTrack,1);
};

#endif
