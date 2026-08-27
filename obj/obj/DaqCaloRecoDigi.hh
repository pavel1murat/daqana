// one-channel reco - calibration only
#ifndef __daqana_obj_DaqCaloRecoDigi_hh__
#define __daqana_obj_DaqCaloRecoDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCaloRecoDigi : public TObject {
public:
  int    sipmid;
  int    ndf;
  int    pileup;
  float  edep;
  float  edep_err;
  float  time;
  float  time_err;
  float  chi2;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCaloRecoDigi();
  virtual ~DaqCaloRecoDigi();

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCaloRecoDigi,1);
};

#endif
