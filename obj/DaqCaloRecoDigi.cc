//

#include "daqana/obj/DaqCaloRecoDigi.hh"

ClassImp(DaqCaloRecoDigi)

//-----------------------------------------------------------------------------
DaqCaloRecoDigi::DaqCaloRecoDigi() : TObject(),
  sipmid(-1),
  ndf(-1),
  pileup(-1),
  edep(-1.),
  edep_err(-1.),
  time(-1.),
  time_err(-1),
  chi2(-1.)
{
}

//-----------------------------------------------------------------------------
DaqCaloRecoDigi::~DaqCaloRecoDigi() {
}

//-----------------------------------------------------------------------------
void DaqCaloRecoDigi::Clear(const char* Opt) {
}

