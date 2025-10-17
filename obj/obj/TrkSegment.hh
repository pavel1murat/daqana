#ifndef __daqana_mod_TrkSegment_hh_
#define __daqana_mod_TrkSegment_hh_
#include <format>
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "TGeoMatrix.h"

#include "daqana/obj/ComboHitData_t.hh"

struct Par_t {
  double a;
  double b;
  double tau;
  double nx;
  double ny;
  double chi2dof;

  double DyDx() const { return  a ; }
  double Y0  () const { return  b ; }
  double X0  () const { return  0 ; }
  double Nx  () const { return  nx; }
  double Ny  () const { return  ny; }
  double Nux () const { return -ny; }
  double Nuy () const { return  nx; }
  double T0  () const { return tau; }
};

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
  double tcorr;
  double t;                             // time-tprop (-fTMean)
  int    drs;
  int    sign_fixed;                    // if 1, the drift sign is not updated any more

  double fR;                            // updated 

  static double const fgVDrift;

  Point2D(int Sid, int Mask, double X, double Y, double Time, double TProp, double TCorr, int Sign = 0) {
    sid        = Sid;
    fMask      = Mask;
    x          = X;
    y          = Y;
    time       = Time;
    tprop      = TProp;
    tcorr      = TCorr;
    drs        = Sign;                // 0:undefined
    sign_fixed = 0;
  }

  // double        r() { return fR; } // ## TODO 

  int IsGood() { return (fMask == 0); }

  mu2e::StrawId strawId() const { return mu2e::StrawId(sid); }

  void print(std::ostream& Stream = std::cout) {
    // double rdrift = this->r();
    Stream << std::format("Point2D: sid:0x{:04x} mask:0x{:04x}",sid,fMask) 
           << std::format(" x,y,time,tprop,tcorr, drs,fixed: {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:2d} {}",
                          x,y,time,tprop,tcorr,drs,sign_fixed)
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
  double               fXMean;
  double               fYMean;
  double               fTMean;                      // <t-tprop>
  std::vector<Point2D> points;               // local points
  int                  fIhit[2];             // indices of the two hits corresponding to the layer transition
                                             // // 
                                                // fits
  Par_t                fTangentLine;
  Par_t                fPar4;                   // best parameters after a 4-point fit
  Par_t                fPar;                    // best fit parameters

  static int           _debugMode;

  TrkSegment(int Plane = -1, int Panel = -1);
  ~TrkSegment() { if (fCombiTrans) delete fCombiTrans ; }

  void addPoint(int Sid, int Flag, double XLoc, double YLoc, double Time, double TProp, double TCorr, int Sign) {
    Point2D pt(Sid, Flag, XLoc,YLoc, Time, TProp, TCorr, Sign);
    points.push_back(pt);
  }

  int InitHits(std::vector<ComboHitData_t>& Hits, int UniquePlane = -1, int Panel = -1);

  int InitHits(); // does the same starting from TrkStrawHitSeeds

  void setSign(int I, int Sign) {
    points.at(I).drs = Sign;
  }

  int      nHits() { return points.size(); }
  int      Plane() { return plane; }
  int      Panel() { return panel; }

  //  int      DefineTangentLine();
//-----------------------------------------------------------------------------
// the segment parameters need always to be redefined consistently
//-----------------------------------------------------------------------------
                                        // line goes through (x0(), y0())
  double   X0     () const { return    0;  }
  double   Y0     () const { return  fPar.b;   }
  double   DyDx   () const { return  fPar.a;   }
  double   T0     () const { return  fPar.tau; }
  double   Nx     () const { return  fPar.nx;  }
  double   Ny     () const { return  fPar.ny;  }
  double   Nux    () const { return -fPar.ny;  }
  double   Nuy    () const { return  fPar.nx;  }
  double   Chi2Dof() const { return  fPar.chi2dof;}

  double  DriftTime (const Point2D* Pt, double T0)  const {
    double t    = Pt->t-T0;
    return t;
  }

  double  R(const Point2D* Pt, double T0)  const {
    double r    = DriftTime(Pt,T0)*Point2D::fgVDrift;
    if (r < 0) {
      std::cout << ">>> ERROR: negative radius: point:" << " r:" << r << std::endl;
      r = 0;
    }
    return r;
  }
                                        // parameters of the line have to be defined
  double Rho(const Point2D* Pt, const Par_t* Par) const {
    double rdrift = R(Pt,Par->T0());
    double rho  = (Par->X0()-Pt->x)*Par->Nux()+(Par->Y0()-Pt->y)*Par->Nuy() - rdrift*Pt->drs;
    return rho;
  }

  int      DefineTangentLine();
                                        // recalculate chi2 - don't have a formula for that 
  double      Chi2();
                                        // if Chi2Dof=-1, recalculate chi2
  void     UpdateParameters(double DyDx, double Y0, double T0, double Chi2Dof = -1);

  void     print   (const char* Message = nullptr, std::ostream& Stream = std::cout);

};


#endif
