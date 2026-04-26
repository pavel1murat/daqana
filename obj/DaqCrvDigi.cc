//

#include "daqana/obj/DaqCrvDigi.hh"

ClassImp(DaqCrvDigi)

//-----------------------------------------------------------------------------
DaqCrvDigi::DaqCrvDigi() : TObject() {
}

//-----------------------------------------------------------------------------
DaqCrvDigi::~DaqCrvDigi() {
}

//-----------------------------------------------------------------------------
int DaqCrvDigi::Init(int Ns) {
  adc.resize(Ns);
  return 0;
}

//-----------------------------------------------------------------------------
void DaqCrvDigi::Clear(const char* Opt) {
}

