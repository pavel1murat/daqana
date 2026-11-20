//

#include "daqana/art_obj/SegmentHit.hh"

//
double       SegmentHit::fgTOffset =  0.0;      // 1.6;       // ns
double       SegmentHit::fgVDrift  = 61.3e-3;   // mm/ns

SegmentHit::SegmentHit() {                        // default constructor needed by the I/O
  fIndex     = -1;
  fStraw     = -1;
  fComboHit  = nullptr;
  fDoca      = 0;
  fMask      = 0;
  x          = 0;
  y          = 0;
  t          = 0;
  fSigmaR    = 0.2;                   // 200 um, for simplicity
  drs        = -999;
  sign_fixed = 0;
}

//-----------------------------------------------------------------------------
SegmentHit::SegmentHit(int Index, mu2e::StrawId Sid) : SegmentHit() {
  fIndex = Index;
  fSid   = Sid.asUint16();
  fStraw = Sid.straw();
}


//-----------------------------------------------------------------------------
void SegmentHit::print(const mu2e::ComboHit* Ch, std::ostream& Stream) {
  Stream << std::format("Point2D: sid:0x{:04x} mask:0x{:04x}",Ch->strawId().asUint16(),fMask)
         << std::format(" x,y,time,tprop,tcorr, drs,fixed: {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:2d} {}",
                        x,y,Ch->time(),Ch->propTime(),Ch->correctedTime(),drs,sign_fixed)
         << std::endl;
}

