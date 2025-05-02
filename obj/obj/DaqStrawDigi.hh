#ifndef __daqana_obj_DaqStrawDigi_hh__
#define __daqana_obj_DaqStrawDigi_hh__

#include "TClonesArray.h"
#include "TObject.h"

class DaqStrawDigi : public TObject {
public:
  int    sid;
  int    tdc0;
  int    tdc1;
  int    tot0;
  int    tot1;
  int    pmp;
  int    flag;
  int    ns;        // N(ADC samples)
  //  short* adc;       // [ns] ns shorts

  DaqStrawDigi();
  DaqStrawDigi(int ns);

void     Init(int Ns);
  virtual ~DaqStrawDigi();

ClassDef(DaqStrawDigi,1);
};

#endif
