#ifndef __daqana_obj_DaqCrvDigi_hh__
#define __daqana_obj_DaqCrvDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqCrvDigi : public TObject {
public:
  int    sbid;
  int    tdc;
  int    nzs;
  int    odd_ts;
  int    sipm;
  int    roc;
  int    feb;
  int    chan;
  std::vector<uint16_t> adc;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqCrvDigi();
  DaqCrvDigi(int ns);
  virtual ~DaqCrvDigi();

  void     Init(int Ns);

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqCrvDigi,1);
};

#endif
