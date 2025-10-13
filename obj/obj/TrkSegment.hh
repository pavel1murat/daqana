#ifndef __daqana_mod_TrkSegment_hh_
#define __daqana_mod_TrkSegment_hh_
#include <format>
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "TGeoMatrix.h"

#include "daqana/obj/ComboHitData_t.hh"

//-----------------------------------------------------------------------------
// 'per-panel' track segments
//-----------------------------------------------------------------------------
struct Point2D {
  int    sid;                           // integer part of straw ID
  int    fMask;                          // 0: ok , 0x1: bad
  double x;
  double y;
  double time;                          // measured time
  double tprop;                         // propagation time
  double t;                             // time-tprop (-fTMean)
  int    drs;
  int    sign_fixed;                    // if 1, the drift sign is not updated any more

  double fR;                            // updated 

  static double const fgVDrift;

  Point2D(int Sid, int Mask, double X, double Y, double Time, double TProp, int Sign = 0) {
    sid        = Sid;
    fMask      = Mask;
    x          = X;
    y          = Y;
    time       = Time;
    tprop      = TProp;
    drs        = Sign;                // 0:undefined
    sign_fixed = 0;
  }

  double        r() { return fR; } // ## TODO 

  int IsGood() { return (fMask == 0); }

  mu2e::StrawId strawId() const { return mu2e::StrawId(sid); }

  void print(std::ostream& Stream = std::cout) {
    // double rdrift = this->r();
    Stream << std::format("Point2D: sid:0x{:04x} mask:0x{:04x}",sid,fMask) 
           << std::format(" x,y,time,tprop,drs,fixed: {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:2d} {}",
                          x,y,time,tprop,drs,sign_fixed)
           << std::endl;
  }
};

class TrkSegment {
public:

  enum {
    kNTransitionsBit   = 0x0001,
    kNHitsBit          = 0x0002,
  };

  static double const  fgRStraw;

  const mu2e::Panel*   trkPanel;
  int                  fMask;           // 0: OK
  int                  plane;
  int                  panel;
  int                  fNTransitions;
  int                  fNGoodHits;
  TGeoCombiTrans*      fCombiTrans;                 // global to local transform
  std::vector<const mu2e::TrkStrawHitSeed*> hits;   // initialization in the ntuple making code
  double               fT0;                         // segment T0
  double               fXMean;
  double               fYMean;
  double               fTMean;                      // <t-tprop>
  double               fNx;
  double               fNy;
  // double               fNux;
  // double               fNuy;
  std::vector<Point2D> points;                 // local points
  int                  fIhit[2];                // indices of the two hits corresponding to the layer transition
  double               fY0;                  // Y(X=0) (local coordinate system of the segment)
  double               fDyDx;                // line tangent to the two key hits, [0]: intercept [1]: slope dydx
  double               fChi2Dof;
  //  double               fDrho;

  static int           _debugMode;

  TrkSegment(int Plane = -1, int Panel = -1);
  ~TrkSegment() { if (fCombiTrans) delete fCombiTrans ; }

  void addPoint(int Sid, int Flag, double XLoc, double YLoc, double Time, double TProp, int Sign) {
    Point2D pt(Sid, Flag, XLoc,YLoc, Time, TProp, Sign);
    points.push_back(pt);
  }

  int InitHits(std::vector<ComboHitData_t>& Hits, int UniquePlane = -1, int Panel = -1);

  void setSign(int I, int Sign) {
    points.at(I).drs = Sign;
  }

  int      nHits() { return points.size(); }

  int      determineDriftSigns     ();
  int      UpdateDriftRadii        ();
                                        // line goes through (x0(), y0())
  double   X0  () { return    0;  }
  double   Y0  () { return  fY0;  }
  double   DyDx() { return  fDyDx;}
  double   Nx  () { return  fNx;  }
  double   Ny  () { return  fNy;  }
  double   Nux () { return -fNy;  }
  double   Nuy () { return  fNx;  }
  double   Chi2Dof() { return  fChi2Dof;}

  double   R(int i) {
    Point2D* pt = &points[i];
    double r    = (pt->time-pt->tprop-fT0)*Point2D::fgVDrift;
    return r;
  }

  void     UpdateLineParameters(double DyDx, double Y0, double T0) {
    fY0   = Y0;
    fDyDx = DyDx;
    fNx   =  1.   /sqrt(1+fDyDx*fDyDx);
    fNy   =  fDyDx/sqrt(1+fDyDx*fDyDx);
    fT0   = T0;
  }

  int     CalculateChi2();

  void     print   (std::ostream& Stream = std::cout);

};


#endif
