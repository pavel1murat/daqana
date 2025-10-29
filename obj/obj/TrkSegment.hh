#ifndef __daqana_mod_TrkSegment_hh_
#define __daqana_mod_TrkSegment_hh_
#include <format>
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"

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
  int    fMask;                         // 0: ok , 0x1: bad
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
  static double const fgTOffset;

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
    kLargeDrhoBit      = 0x0004,
  };

  static double const  fgRStraw;

  const mu2e::Panel*   fTrkPanel;
  int                  fMask;           // 0: OK
  int                  fPlane;
  int                  fPanel;
  int                  fNTransitions;
  int                  fNGoodHits;                  // total number of good hits
  int                  fNghl[2];                    // # good hits in each layer
  int                  fNmhl[2];                    // # of straws w/o hits in each layer
  std::vector<const mu2e::ComboHit*> hits;          // initialization (from a time cluster) in the ntuple making code
  double               fXMean;
  double               fYMean;
  double               fTMean;                      // <t-tprop>
  std::vector<Point2D> points;               // local points
  int                  fIhit[2];             // indices of the two hits corresponding to the layer transition
                                                // fits
  Par_t                fTangentLine;
  Par_t                fPar4;                   // best parameters after a 4-point fit
  Par_t                fPar;                    // best fit parameters

  static int           fgDebugMode;

  TrkSegment(int Plane = -1, int Panel = -1);
  ~TrkSegment() { }

  void Clear();

  Point2D* Point(int I) { return &points[I] ; }

  void AddPoint(int Sid, int Mask, double XLoc, double YLoc, double Time, double TProp, double TCorr, int Sign) {
    Point2D pt(Sid, Mask, XLoc,YLoc, Time, TProp, TCorr, Sign);
    points.push_back(pt);
  }

  int InitHits(std::vector<const mu2e::ComboHit*>* Hits = nullptr, int UniquePlane = -1, int Panel = -1);

  // int InitHits(); // does the same starting from TrkStrawHitSeeds

  // void setSign(int I, int Sign) {
  //   points.at(I).drs = Sign;
  // }
                                        // any set bit marks the hit as bad, so far use 'dead'
  int      GoodHit(const mu2e::ComboHit* Ch) {
    int x  = std::stoi(Ch->_flag.hex(),nullptr,0);
    int good = (x == 0);
    return good;
  }

  int      nHits() { return (int) hits.size(); }
  int      Plane() { return fPlane; }
  int      Panel() { return fPanel; }

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

//-----------------------------------------------------------------------------
// by default, use the segment T0
//-----------------------------------------------------------------------------
  double  DriftTime (const Point2D* Pt, double HitT0=1.e12)  const {
    double t0 = HitT0;
    if (t0 > 1.e10) t0 = T0();

    double t    = Pt->t-t0 + Point2D::fgTOffset;
    return t;
  }

  double  R(const Point2D* Pt, double T0=1.e12)  const {
    static int nerr(0), last_printed(0), nprinted(0), step(1);
    double r    = DriftTime(Pt,T0)*Point2D::fgVDrift;
    if (r < 0) {
      nerr++;
      if (nerr-last_printed >= step) {
        std::cout << ">>> ERROR: nerr:" << std::setw(14) << nerr <<" negative radius: point:" << " r:" << r << std::endl;
        nprinted++;
      }
      if (nprinted > 10) {
        step      = step*2;
        nprinted  = 0;
      }
      r = 0;
    }
    return r;
  }
//-----------------------------------------------------------------------------
// signed distance between the drift circle and the segment (DOCA),
// definition of the sign - to be clarified
// use segment T0
// parameters of the line have to be defined
//-----------------------------------------------------------------------------
  double Rho(const Point2D* Pt, const Par_t* Par = nullptr) const {
    const Par_t* par(Par);
    if (par == nullptr) par = &fPar;

    double rdrift = R(Pt,par->T0());
    double rho  = (par->X0()-Pt->x)*par->Nux()+(par->Y0()-Pt->y)*par->Nuy() - rdrift*Pt->drs; // 
    return rho;
  }

//-----------------------------------------------------------------------------
// segment_to_wire_distance.vdot.normal_to_the_segment
//-----------------------------------------------------------------------------
  double SwDist(const Point2D* Pt, const Par_t* Par = nullptr) const {
    const Par_t* par(Par);
    if (par == nullptr) par = &fPar;

    double rho  = (par->X0()-Pt->x)*par->Nux()+(par->Y0()-Pt->y)*par->Nuy(); // 
    return rho;
  }

  int      DefineTangentLine();
                                        // recalculate chi2 - don't have a formula for that 
  double   Chi2();
                                        // if Chi2Dof=-1, recalculate chi2
  void     UpdateParameters(double DyDx, double Y0, double T0, double Chi2Dof = -1);

  void     print   (const char* Message = nullptr, std::ostream& Stream = std::cout);

};


#endif
