#ifndef __SegmentFit_hh__
#define __SegmentFit_hh__

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "daqana/obj/TrkSegment.hh"

class SegmentFit {
public:

  struct Par_t {
    double a;
    double b;
    double tau;
    double chi2dof;
  };

  static int fgDebugMode;

  TrkSegment* fSegment;                 // cached segment

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
    double smn = fSts.yMean();
    double vd  = Point2D::fgVDrift;

    double sig_yss = fSys.xMean()-fSys.yMean()*fSys.sigXY();
    double sig_xss = fSxs.xMean()-fSxs.yMean()*fSxs.sigXY();

    double x = (sig_yss-A*sig_xss-vd*sqrt(1+A*A)*fSts.sigXY())/(1-smn*smn);

    return x;
  }

  double DbDa(double A) {
    double smn = fSts.yMean();
    double vd  = Point2D::fgVDrift;

    double sig_xss = fSxs.xMean()-fSxs.yMean()*fSxs.sigXY();

    double x = (-1.)*(sig_xss+A*vd*fSts.sigXY()/sqrt(1+A*A))/(1-smn*smn);
    return x;
  }

  double Tau(double A) {
    double smn = fSts.yMean();
    double vd  = Point2D::fgVDrift;
    double tau = fSts.xMean()+1./(1-smn*smn)*((A*fSxs.sigXY()-fSys.sigXY())/(vd*sqrt(1+A*A))-smn*fSts.sigXY());
    return tau;
  }

  double DtauDa(double A) {
    double smn = fSts.yMean();
    double vd  = Point2D::fgVDrift;
    double x   = (fSxs.sigXY()+A*fSys.sigXY())/(vd*(1-smn*smn)*pow(1+A*A,3/2.));
    return x;
  }

  int Fit(int NIterations,Par_t* Par);
  int Init();
                                        // find the segment line parameters using two seed hits and two edge hits
  int DefineDriftDirections();
  int FindTangentLine (TrkSegment* Seg);
  int CalculateLsqSums();
  int DisplaySegment();
};

#endif
