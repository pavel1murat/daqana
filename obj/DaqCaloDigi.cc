//

#include "daqana/obj/DaqCaloDigi.hh"

ClassImp(DaqCaloDigi)

//-----------------------------------------------------------------------------
DaqCaloDigi::DaqCaloDigi() : TObject() {
  _ns = -1;
}

//-----------------------------------------------------------------------------
DaqCaloDigi::~DaqCaloDigi() {
}

//-----------------------------------------------------------------------------
int DaqCaloDigi::Init(int Ns) {
  int rc(0);
  if (_ns != Ns) {
    _ns = Ns;
    wf.resize(Ns);
  }
  return 0;
}

//-----------------------------------------------------------------------------
void DaqCaloDigi::Clear(const char* Opt) {
}

