#ifndef __daqana_obj_DaqCaloDigi_hh__
#define __daqana_obj_DaqCaloDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCaloDigi : public TObject {
public:
  int    _ns;
  int    sipmid;
  int    t0;
  int    ppos;                          // peak position
  std::vector<uint16_t> wf;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCaloDigi();
  DaqCaloDigi(int ns);
  virtual ~DaqCaloDigi();

  int   ns() { return _ns; } // return adc.size(); }
  
  int     Init(int Ns);

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCaloDigi,1);
};

#endif
