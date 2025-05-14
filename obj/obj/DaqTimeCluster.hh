#ifndef __daqana_obj_DaqTimeCluster_hh__
#define __daqana_obj_DaqTimeCluster_hh__

#include "TObject.h"

class DaqTimeCluster : public TObject {
public:
  int     nsh;
  int     nch;
  float   t0;
  float   tmin;
  float   tmax;

  DaqTimeCluster();

  virtual ~DaqTimeCluster();

  ClassDef(DaqTimeCluster,1);
};

#endif
