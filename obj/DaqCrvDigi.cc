//

#include "daqana/obj/DaqCrvDigi.hh"

ClassImp(DaqCrvDigi)

//-----------------------------------------------------------------------------
DaqCrvDigi::DaqCrvDigi() : TObject() {
  _ns = -1;
}

//-----------------------------------------------------------------------------
DaqCrvDigi::~DaqCrvDigi() {
}

//-----------------------------------------------------------------------------
int DaqCrvDigi::Init(int Ns) {
  if (_ns != Ns) {
    _ns = Ns;
    adc.resize(Ns);
  }
  return 0;
}

//-----------------------------------------------------------------------------
void DaqCrvDigi::Clear(const char* Opt) {
}
