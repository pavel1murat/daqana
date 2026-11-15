//-----------------------------------------------------------------------------
// SegmentHit: a 2D hit of a track segment reconstructed in a single panel
//-----------------------------------------------------------------------------
#ifndef __daqana_art_obj_SegmentHit_hh_
#define __daqana_art_obj_SegmentHit_hh_
#include <format>
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
//-----------------------------------------------------------------------------
// representation of a 2D hit used in the segment fit
//-----------------------------------------------------------------------------
struct SegmentHit {
  int                     fIndex;                        // index of the ComboHit in the ComboHitCollection of the segment, -1: no ComboHit
  int                     fStraw;                        // straw number in the panel
  uint16_t                fSid;                          // integer part of Sid
  const mu2e::ComboHit*   fComboHit;                     // the ComboHit, set from fIndex by TrkSegment::setComboHits, transient, null if fIndex<0
  double                  fDoca;                         // distance between the segment and the wire, after the final fit
  int                     fMask;                         // 0: ok , 0x1: bad
  double                  x;
  double                  y;
  double                  t;                             // time-tprop (-fTMean)
  int                     drs;
  int                     sign_fixed;                    // if 1, the drift sign is not updated any more

  double                  fR;                            // updated
  double                  fSigmaR;                       // hit coordinate resolution
  double                  fChi2;                         // hit contribution to the segment chi^2

                                                         // drift velocity parameterization : R = A*(T-T0)
  static double           fgVDrift;                      // mm/ns, defined in TrkSegment.cc
  static double           fgTOffset;

  SegmentHit() {                        // default constructor needed by the I/O
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

  explicit SegmentHit(int Index, mu2e::StrawId Sid) : SegmentHit() {
    fIndex = Index;
    fSid   = Sid.asUint16();
    fStraw = Sid.straw();
  }

  // double        r() { return fR; } // ## TODO

  const mu2e::ComboHit* ComboHit() const { return fComboHit; }
  double                Doca    () const { return fDoca;     }
  int                   IsGood  () const { return (fMask == 0); }
  int                   Straw   () const { return fStraw; }
  int                   Layer   () const { return fStraw % 2 ; }
  mu2e::StrawId         Sid     () const { return mu2e::StrawId(fSid); }

  void print(const mu2e::ComboHit* Ch, std::ostream& Stream = std::cout) {
    Stream << std::format("Point2D: sid:0x{:04x} mask:0x{:04x}",Ch->strawId().asUint16(),fMask)
           << std::format(" x,y,time,tprop,tcorr, drs,fixed: {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:2d} {}",
                          x,y,Ch->time(),Ch->propTime(),Ch->correctedTime(),drs,sign_fixed)
           << std::endl;
  }

  static void SetTOffset(double T) { fgTOffset = T; }
  static void SetVDrift (double V) { fgVDrift  = V; }
//-----------------------------------------------------------------------------
// R in mm, t in ns ; implement simplest non-linearity
//-----------------------------------------------------------------------------
  static double T2D(double T) {
    double r = (T+fgTOffset)*fgVDrift;
    return r;
  }

};

#endif
