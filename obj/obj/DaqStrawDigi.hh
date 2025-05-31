#ifndef __daqana_obj_DaqStrawDigi_hh__
#define __daqana_obj_DaqStrawDigi_hh__

#include "TClonesArray.h"
#include "TObject.h"

class DaqStrawDigi : public TObject {
public:
  int    sid;
  int    mnid;                          // 
  int    tdc0;
  int    tdc1;
  int    tot0;
  int    tot1;
  int    pmp;
  int    flag;
  int    fs;
  float  bl;
  float  ph;
  int    ns;        // N(ADC samples)
  short* adc;       //[ns] ns shorts
  

  DaqStrawDigi();
  DaqStrawDigi(int ns);
  virtual ~DaqStrawDigi();

  void     Init(int Ns);

  int      Straw() {return (sid      ) & 0x7f; }
  int      Panel() {return (sid >>  7) & 0x07; } 
  int      Plane() {return (sid >> 10) & 0x3f; } 

  ClassDef(DaqStrawDigi,1);
};

#endif
