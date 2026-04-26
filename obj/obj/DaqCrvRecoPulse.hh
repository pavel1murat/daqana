#ifndef __daqana_obj_DaqCrvRecoPulse_hh__
#define __daqana_obj_DaqCrvRecoPulse_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCrvRecoPulse : public TObject {
public:
  float    npes;
  float    pes_ph;
  float    time;

  float    ph;
  float    ped;                         // pedestal;
  float    beta;
  float    chi2;
  float    le_time;
  int      flags;
  
  float    npes_nofit;
  float    time_nofit;
  float    tstart;
  float    tend;

  int      sbid;                        // scintillator bar index
  int      sipm;
  int      roc;
  int      feb;
  int      ch;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCrvRecoPulse();
  virtual ~DaqCrvRecoPulse();

  int      Init();

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCrvRecoPulse,1);
};

#endif
