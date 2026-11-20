//-----------------------------------------------------------------------------
// TrkSegment: a track segment reconstructed in a single panel
//-----------------------------------------------------------------------------
#ifndef __daqana_art_obj_TrkSegment_hh_
#define __daqana_art_obj_TrkSegment_hh_
#include <format>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/art_obj/SegmentHit.hh"

//-----------------------------------------------------------------------------
class TrkSegment {
public:

  struct Par_t {
    double a;
    double b;
    double tau;
    double nx;                          // nx > 0 (convention)
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

  enum {
    kNTransitionsBit   = 0x0001,
    kNHitsBit          = 0x0002,
    kLargeDrhoBit      = 0x0004,
    kSubsegmentBit     = 0x0008,
    kNoComboHitBit     = 0x0010,        // a straw of the segment without a ComboHit
  };

  static int              fgDebugMode;
  static double const     fgRStraw;        // straw radius

  mu2e::Panel*            fTrkPanel;
  int                     fIndex;                  // segment number in the list
  int                     fMask;                   // 0: OK
  int                     fPlane;
  int                     fPanel;
  int                     fNTransitions;
  int                     fNGoodHits;               // total number of good hits
  int                     fNghl[2];                 // # good hits in each layer
  int                     fNmhl[2];                 // # of straws w/o hits in each layer
                                                    // SegmentHit replaces Point2D
  std::vector<SegmentHit> fListOfHits;              // initialization (from a time cluster) in the ntuple making code
  double                  fXMean;
  double                  fYMean;
  double                  fTMean;                   // <t-tprop>
  int                     fIhit[2];                 // indices of the two hits corresponding to the layer transition
                                                    // fits
  Par_t                   fTangentLine;
  Par_t                   fPar4;                    // best parameters after a 4-point fit
  Par_t                   fPar;                     // best fit parameters

  TrkSegment(int Plane = -1, int Panel = -1);
  ~TrkSegment() { }

  void Clear();

  SegmentHit*  Hit(int I) const { return (SegmentHit*) &fListOfHits[I]; }

  // this is a list of all straws in between the first and the last hits, some straws could be empty
  int          NStraws()    const { return fListOfHits.size(); }

                                        // set the ComboHit of each hit from its index in the collection;
                                        // a null collection clears the pointers
  void setComboHits(mu2e::ComboHitCollection const* Hits);
  int  NHitsWithComboHit() const;
                                        // one hit per straw between the first and the last hit (hits sorted by straw);
                                        // a straw without a ComboHit gets a hit with fIndex=-1
  void fillMissingStraws();
                                        // distance between the segment and the wire of each hit, after the final fit
  void UpdateHitParameters();

                                        // rotate by default (DoTransform=1)

  int InitHits(int DoTransform = 1);

                                        // any set bit marks the hit as bad, so far use 'dead'

  int      GoodHit(const mu2e::ComboHit* Ch) const {
    int x  = std::stoi(Ch->_flag.hex(),nullptr,0);
    int good = (x == 0);
    return good;
  }
                                        // returns next straw with combohit after I
  SegmentHit* NextHit(int I) const;

  int      Plane() const { return fPlane; }
  int      Panel() const { return fPanel; }
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

  double   t0_rec () const { return  fPar.tau+fTMean; }
  double   y0_rec () const { return  fYMean-fPar.a*fXMean+fPar.b; }
//-----------------------------------------------------------------------------
// by default, use the segment T0
//-----------------------------------------------------------------------------
  double  DriftTime (const SegmentHit* Pt, double HitT0=1.e12)  const {
    double t0 = HitT0;
    if (t0 > 1.e10) t0 = T0();

    double t    = Pt->t-t0 + SegmentHit::fgTOffset;
    return t;
  }

//-----------------------------------------------------------------------------
  double  R(const SegmentHit* Pt, double HitT0=1.e12)  const {
    static int nerr(0), last_printed(0), nprinted(0), step(1);

    double t0 = HitT0;
    if (t0 > 1.e10) t0 = T0();

    double tdrift = Pt->t-t0;
    double r      = SegmentHit::T2D(tdrift);

    if (r < 0) {
      nerr++;
      if (nerr-last_printed >= step) {
        // std::cout << ">>> ERROR: nerr:" << std::setw(14) << nerr <<" negative radius: point:" << " r:" << r << std::endl;
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
// (R_seg-R_wire).vdot.Nu_seg
// Nu_seg - normal to the segment, 90 deg counterclockwise from the segment direction)
//-----------------------------------------------------------------------------
  double Doca(const SegmentHit* Pt, const Par_t* Par = nullptr) const {
    const Par_t* par(Par);
    if (par == nullptr) par = &fPar;

    double doca = (par->X0()-Pt->x)*par->Nux()+(par->Y0()-Pt->y)*par->Nuy(); // 
    return doca;
  }

//-----------------------------------------------------------------------------
// signed distance between the drift circle and the segment (DOCA),
// definition of the sign - to be clarified, but looks projected to (Nux(),Nuy())
// use segment T0
// parameters of the line have to be defined
//-----------------------------------------------------------------------------
  double Dr(const SegmentHit* Pt, const Par_t* Par = nullptr) const {
    const Par_t* par(Par);
    if (par == nullptr) par = &fPar;

    double rdrift = R(Pt,par->T0());
    double dr     = (par->X0()-Pt->x)*par->Nux()+(par->Y0()-Pt->y)*par->Nuy() - rdrift*Pt->drs; // 
    return dr;
  }

  int      DefineTangentLine();
                                        // recalculate chi2 - don't have a formula for that 
  double   Chi2();

  void     SetIndex(int I) { fIndex = I; }
                                        // if Chi2Dof=-1, recalculate chi2
  void     UpdateParameters(double DyDx, double Y0, double T0, double Chi2Dof = -1);

  void     print   (const char* Message = nullptr, std::ostream& Stream = std::cout);

};

using TrkSegmentCollection = std::vector<TrkSegment>;

#endif
