//

#include "daqana/obj/DaqStrawDigi.hh"
// #include "DaqStrawDigi.hh"

ClassImp(DaqStrawDigi)

//-----------------------------------------------------------------------------
DaqStrawDigi::DaqStrawDigi() : TObject() {
  ns = -1;
}

//-----------------------------------------------------------------------------
int DaqStrawDigi::Init(int Ns) {
  

  // if (cleared != 1) {
  //   ns  = Ns;
  //   adc = new uint16_t[Ns];
  // }
  if (ns < 0) {
    adc.resize(Ns);
    ns = Ns;
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
