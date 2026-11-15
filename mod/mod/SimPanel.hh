///////////////////////////////////////////////////////////////////////////////
// simulation helpers shared by the SimHits and SimRes modules: the geometry of a
// single tracker panel, a simulated straight track, and the straws hit by the track
//
// DebugMode:
//   bit0: printouts of the line parameters and of the hit straws
//   bit1: hits
//   bit2: geometry
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_mod_SimPanel_hh__
#define __daqana_mod_SimPanel_hh__

#include <cmath>
#include <format>
#include <iostream>
#include <vector>

#include "cetlib_except/exception.h"
#include "Offline/DataProducts/inc/StrawId.hh"

namespace mu2e {
namespace simpanel {

struct Straw_t {
  StrawId id;
  double  fX;
  double  fY;
  double  fDoca;                        // track to wire distance, for the hit straws
};

struct Line_t {
  double fX0;
  double fY0;
  double fDyDx;
  double fT0;                           // T0 of the simulated track
  double fNx;
  double fNy;
  double fNux;
  double fNuy;
};

//-----------------------------------------------------------------------------
// straws of the panel, the layers are staggered; OffsetX and OffsetY are the
// per-straw misalignments, one value per straw, in mm
//-----------------------------------------------------------------------------
inline std::vector<Straw_t> straws(int DebugMode, std::vector<double> const& OffsetX,
                                   std::vector<double> const& OffsetY) {
  if ((OffsetX.size() != 96) or (OffsetY.size() != 96)) {
    throw cet::exception("CONFIG") << "simpanel::straws: strawOffsetX and strawOffsetY need 96 values, got "
                                   << OffsetX.size() << " and " << OffsetY.size() << "\n";
  }
  std::vector<Straw_t> result;
  double pitch = 6.25;                  // step between straws in the layer
  for (int is=0; is<96; ++is) {
    StrawId sid(0,0,is);
    int     lay = sid.layer();
    int     num = is/2;
    Straw_t s{sid, (num-24)*pitch + pitch/2*lay + OffsetX[is], 3.*(2*lay-1) + OffsetY[is], 0.};
    result.push_back(s);
    if (DebugMode & 0x4) {
      std::cout << std::format("is:{:2d} layer:{} num:{:2d} y:{:10.3f} x:{:10.3f}\n",is,lay,num,s.fX,s.fY);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// the simulated line: a fixed line for each event, BeamMode=1 gives the parameters of the odd-reco study
//-----------------------------------------------------------------------------
inline Line_t line(int BeamMode, int DebugMode) {
  Line_t l{};
  l.fT0   = 0.;
  l.fDyDx = 0.4;
  l.fX0   = 0.;
  double x = 0.;                        // [-10,10] mm would cover pretty much the full range
  l.fY0   = -x/l.fDyDx;                 // y(x=0)

  if (BeamMode == 1) {
                                        // these parameters lead to odd reco
    l.fY0   = -22.94303316;
    l.fDyDx = 0.4;
    l.fT0   = 0.;
  }

  l.fNx  = 1./sqrt(1+l.fDyDx*l.fDyDx);
  l.fNy  = l.fDyDx*l.fNx;
  l.fNux = -l.fNy;
  l.fNuy = l.fNx;

  if (DebugMode & 0x1) {
    std::cout << std::format("fY0:{:14.8f} fDyDx:{:14.8f} fT0:{:14.8f} nx:{:10.6f} ny:{:10.6f}\n",
                             l.fY0,l.fDyDx,l.fT0,l.fNx,l.fNy);
  }
  return l;
}

//-----------------------------------------------------------------------------
// the straws within StrawRadius of the line, in the straw order, with the track-to-wire distance
//-----------------------------------------------------------------------------
inline std::vector<Straw_t> hitStraws(std::vector<Straw_t> const& Straws, Line_t const& Line,
                                      double StrawRadius, int DebugMode) {
  std::vector<Straw_t> result;
  for (auto const& s : Straws) {
                                        // distance (signed) from the track to the wire
    double dr = (Line.fY0-s.fY)*Line.fNuy+(Line.fX0-s.fX)*Line.fNux;
    if (fabs(dr) < StrawRadius) {
      Straw_t h = s;
      h.fDoca = fabs(dr);
      result.push_back(h);
      if (DebugMode & 0x2) {
        std::cout << std::format("-- hit in straw:{:2d} x0:{:10.3f}  y0:{:10.3f} doca:{:10.3f}\n",
                                 s.id.straw(),s.fX,s.fY,h.fDoca);
      }
    }
  }
  return result;
}

}  // namespace simpanel
}  // namespace mu2e

#endif
