//

#include "daqana/obj/DaqStrawDigi.hh"
// #include "DaqStrawDigi.hh"

ClassImp(DaqStrawDigi)

//-----------------------------------------------------------------------------
DaqStrawDigi::DaqStrawDigi() : TObject() {
  ns = -1;
  //  adc = nullptr;
}

//-----------------------------------------------------------------------------
DaqStrawDigi::DaqStrawDigi(int Ns) {
  Init(Ns);
}

//-----------------------------------------------------------------------------
void DaqStrawDigi::Init(int Ns) {
  ns = Ns;
  // if (adc == nullptr) {
  //   ns  = Ns;                             // shouldn't change within the job
  //   adc = new int16_t[ns];
  // }
}

//-----------------------------------------------------------------------------
DaqStrawDigi::~DaqStrawDigi() {
  //  if (adc) delete [] adc; 
}
