//

#include "daqana/obj/DaqStrawDigi.hh"

ClassImp(DaqStrawDigi)

//-----------------------------------------------------------------------------
DaqStrawDigi::DaqStrawDigi() : TObject() {
  _ns = -1;
  
  sid  = -1;
  mnid = -1;
}

//-----------------------------------------------------------------------------
DaqStrawDigi::~DaqStrawDigi() {
  //  delete adc; 
}

//-----------------------------------------------------------------------------
int DaqStrawDigi::InitSD(int Ns) {
  int rc(0);
  if (_ns < 0) {
    _ns = Ns;
    adc.resize(Ns);
  }
  return rc;
}

// //-----------------------------------------------------------------------------
// int DaqStrawDigi::InitSD2(int Ns) {
//   int rc(0);
//   // if (_ns < 0) {
//   //   _ns = Ns;
//   //   //    adc.resize(Ns);
//   // }
//   return rc;
// }

//-----------------------------------------------------------------------------
void DaqStrawDigi::Clear(const char* Opt) {
  //cleared = 1;
  // adc.clear();
}
