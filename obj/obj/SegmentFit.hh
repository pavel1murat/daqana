#ifndef __SegmentFit_hh__
#define __SegmentFit_hh__

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "daqana/obj/TrkSegment.hh"

class SegmentFit {
public:

  static int fgDebugMode;
  static int fgDebugBits[100];

  TrkSegment* fSegment;                 // cached segment

  //  Par_t      fP0;

  LsqSums2   fSxy;
  LsqSums2   fSxs;
  LsqSums2   fSys;
  LsqSums2   fSts;

  SegmentFit(TrkSegment* Seg);
  ~SegmentFit();
//-----------------------------------------------------------------------------
// F = dChi2/dA (partial derivative)
// Df/Da = dF/dA, with the second derivatives over b and tau beging zeroes
//-----------------------------------------------------------------------------
  double F   (double A);
  double DfDa(double A);

  double B(double A) {
    double smn     = fSts.yMean();
    double vd      = SegmentHit::fgVDrift;

    double sig_yss = fSys.xMean()-fSys.yMean()*fSys.sigXY();
    double sig_xss = fSxs.xMean()-fSxs.yMean()*fSxs.sigXY();
    double x       = (sig_yss - A*sig_xss + vd*sqrt(1+A*A)*fSts.sigXY())/(1-smn*smn);

    return x;
  }

  double DbDa(double A) {
    double smn = fSts.yMean();
    double vd  = SegmentHit::fgVDrift;

    double sig_xss = fSxs.xMean()-fSxs.yMean()*fSxs.sigXY();

    // double x = (-1.)*(sig_xss+A*vd*fSts.sigXY()/sqrt(1+A*A))/(1-smn*smn);
    double x = (-sig_xss + A*vd*fSts.sigXY()/sqrt(1+A*A))/(1-smn*smn);
    return x;
  }

  double Tau(double A) {
    double smn = fSts.yMean();
    double vd  = SegmentHit::fgVDrift;
    //    double tau = fSts.xMean()+1./(1-smn*smn)*((A*fSxs.sigXY()-fSys.sigXY())/(vd*sqrt(1+A*A))-smn*fSts.sigXY());
    double tau = fSts.xMean() - 1./(1-smn*smn)*((A*fSxs.sigXY() - fSys.sigXY())/(vd*sqrt(1+A*A)) + smn*fSts.sigXY());
    return tau;
  }

  double DtauDa(double A) {
    double smn = fSts.yMean();
    double vd  = SegmentHit::fgVDrift;
    //    double x   = (fSxs.sigXY()+A*fSys.sigXY())/(vd*(1-smn*smn)*pow(1+A*A,3/2.));
    double x   = -(fSxs.sigXY() + A*fSys.sigXY())/(vd*(1-smn*smn)*pow(1+A*A,3/2.));
    return x;
  }
                                        // for debugging/validation

  double DChi2Da  (double A, double B, double Tau);
  double DChi2Db  (double A, double B, double Tau);
  double DChi2Dtau(double A, double B, double Tau);

  int    Fit(int NIterations, int DoCleanup, const TrkSegment::Par_t* Pin, TrkSegment::Par_t* Par);
  int    Init();
                                        // find the segment line parameters using two seed hits and two edge hits
  int         DefineDriftDirections(const TrkSegment::Par_t* Pin = nullptr);
  int         CalculateLsqSums();
  int         DisplaySegment();
                                        // if nullptr, use fTangentLine
  int         DefineTangentLine(TrkSegment::Par_t* Par = nullptr);

  static int  DebugMode()      { return fgDebugMode;    }
  static int  DebugBit (int I) { return fgDebugBits[I]; }

  static void SetDebugMode(int Mode) { fgDebugMode = Mode; }
  static void DebugBit(int I, int Val) { fgDebugBits[I] = Val; }
};

#endif
