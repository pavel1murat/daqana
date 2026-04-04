//

#include "daqana/obj/DaqStrawDigi.hh"
// #include "DaqStrawDigi.hh"

ClassImp(DaqStrawDigi)

//-----------------------------------------------------------------------------
DaqStrawDigi::DaqStrawDigi() : TObject() {
  //  adc = nullptr;
  ns = -1;
}

//-----------------------------------------------------------------------------
void DaqStrawDigi::Init(int Ns) {
  

  // if (cleared != 1) {
  //   ns  = Ns;
  //   adc = new uint16_t[Ns];
  // }
  if (ns < 0) {
    ns = Ns;
    adc = new short[Ns];
  }
}

//-----------------------------------------------------------------------------
DaqStrawDigi::~DaqStrawDigi() {
  //  delete adc; 
}

//-----------------------------------------------------------------------------
void DaqStrawDigi::Clear(const char* Opt) {
  //cleared = 1;
  // adc.clear();
}
