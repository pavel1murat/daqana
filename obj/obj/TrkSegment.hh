#ifndef __daqana_mod_TrkSegment_hh_
#define __daqana_mod_TrkSegment_hh_
#include <format>
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "TGeoMatrix.h"

//-----------------------------------------------------------------------------
// 'per-panel' track segments
//-----------------------------------------------------------------------------
struct Point2D {
  double x;
  double y;
  double r;
  int    drift_sign;
  int    sign_fixed;                  // if 1, the drift sign is not updated any more

  Point2D(double X, double Y, double R, int Sign = 0) {
    x          = X;
    y          = Y;
    r          = R;
    drift_sign = Sign;                // 0:undefined
    sign_fixed = 0;
  }

  void print(std::ostream& Stream = std::cout) {
    Stream << std::format("Point2D: x,y,r,s,fixed: : {:10.3f}: {:10.3f} {:10.3f} {:2d} {}",
                          x,y,r,drift_sign,sign_fixed) << std::endl;
  }
};

class TrkSegment {
public:
  const mu2e::Panel*   trkPanel;
  int                  plane;
  int                  panel;
  TGeoCombiTrans*      combiTrans;             // global to local transform
  std::vector<const mu2e::TrkStrawHitSeed*> hits;
  double               fT0;                     //
  double               ymean;
  double               xmean;
  std::vector<Point2D> points;                 // local points
  int                  ihit[2];                // indices of the two hits corresponding to the layer transition
  double               line[2];                // line tangent to the two key hits, [0]: intercept [1]: slope
  double               chi2;
  double               drho;

  static int           _debugMode;

  TrkSegment();

  void addPoint(double XLoc, double YLoc, double RDrift, int Sign) {
    Point2D pt(XLoc,YLoc,RDrift,Sign);
    points.push_back(pt);
  }

  void setSign(int I, int Sign) {
    points.at(I).drift_sign = Sign;
  }

  int      correctT0               ();
  int      determineDriftSigns     ();
  int      finalFit                ();
  int      findLine                ();
  int      findInitialApproximation();
  int      fitLine                 ();
  int      fit                     ();
  int      updateCoordinates       ();

  void     print   (std::ostream& Stream);

};


#endif
