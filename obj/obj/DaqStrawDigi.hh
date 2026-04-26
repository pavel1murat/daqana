#ifndef __daqana_obj_DaqStrawDigi_hh__
#define __daqana_obj_DaqStrawDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class DaqStrawDigi : public TObject {
public:
  int    _ns;                   //     N(ADC samples)
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
  std::vector<uint16_t> adc;
  //  uint16_t adc[15];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  DaqStrawDigi();
  DaqStrawDigi(int ns);
  virtual ~DaqStrawDigi();

  int      InitSD(int Ns);
  // int      InitSD2(int Ns);

  int      straw() {return (sid      ) & 0x7f; }
  int      panel() {return (sid >>  7) & 0x07; } 
  int      plane() {return (sid >> 10) & 0x3f; }

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(DaqStrawDigi,1);
};

#endif
