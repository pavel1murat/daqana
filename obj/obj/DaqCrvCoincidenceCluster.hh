#ifndef __daqana_obj_DaqCrvCoincidenceCluster_hh__
#define __daqana_obj_DaqCrvCoincidenceCluster_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCrvCoincidenceCluster : public TObject {
public:
  int      stype;                       // sector type
  float    tstart;
  float    tend;
  float    pes;
  float    time;
  float    x;
  float    y;
  float    z;
  int      nlayers;
  int      nsides;                      // 1 or 2
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCrvCoincidenceCluster();
  virtual ~DaqCrvCoincidenceCluster();

  int     Init();

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCrvCoincidenceCluster,1);
};

#endif
