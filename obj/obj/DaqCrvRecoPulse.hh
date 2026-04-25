#ifndef __daqana_obj_DaqCrvRecoPulse_hh__
#define __daqana_obj_DaqCrvRecoPulse_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCrvRecoPulse : public TObject {
public:
  float    npe;
  float    pe_ph;
  float    pe_time;
  float    ph;
  float    beta;
  float    chi2;
  float    le_time;
  float    chan;
  int      flags;
  float    npe_nofit;
  float    time_nofit;
  float    tstart;
  float    tend;

  int      sbid;                        // scintillator bar index
  int      roc;
  int      feb;
  int      chan;

  float    ped;                         // pedestal;
  bool     ped_from_db;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCrvRecoPulse();i
  DaqCrvRecoPulse(int ns);
  virtual ~DaqCrvRecoPulse();

  void     Init(int Ns);

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCrvRecoPulse,1);
};

#endif
