#ifndef __daqana_obj_DaqCaloCluster_hh__
#define __daqana_obj_DaqCaloCluster_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"
#include "TVector3.h"

class DaqCaloCluster : public TObject {
public:
  uint8_t  disk;                         // disk
  uint8_t  split;                        // 
  short    size;
  float    time;                            // N(digis)
  float    sigt;
  float    edep;
  float    sige;
  float    x;                           // coordinates of the COG
  float    y;
  float    z;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCaloCluster();
  virtual ~DaqCaloCluster();

  int     Init(int Ns);

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCaloCluster,1);
};

#endif
