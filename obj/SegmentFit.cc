///////////////////////////////////////////////////////////////////////////////
// fgDebugMode = 1: enables debug printouts, some printouts require specific 
//                  fgDebugBits to be non-zero
// fgDebugBits[0] = 1 : pritnout from DisplaySegment
///////////////////////////////////////////////////////////////////////////////
#include "TLine.h"
#include "TEllipse.h"

#include "daqana/obj/SegmentFit.hh"
#include <memory>


int SegmentFit::fgDebugMode{0};
int SegmentFit::fgDebugBits[100] = {}; 

//-----------------------------------------------------------------------------
SegmentFit::SegmentFit(TrkSegment* Seg) {
  fSegment = Seg;
  Init();
}

//-----------------------------------------------------------------------------
SegmentFit::~SegmentFit() {
}

//-----------------------------------------------------------------------------
int SegmentFit::Init() {
  CalculateLsqSums();
  return 0;
}

//-----------------------------------------------------------------------------
// find line tangent to two circles defined by the two hits with known drift signs
//-----------------------------------------------------------------------------
double SegmentFit::F(double A) {
  double vd  = SegmentHit::fgVDrift;

  int npt = fSegment->nHits();

  double sum = 0;

  double b   = B(A);
  double tau = Tau(A);

  for (int i=0; i<npt; i++) {
    SegmentHit* pt = fSegment->Hit(i);
    if (pt->IsGood() != 1) continue;

    //    double ds = ((pt->y-A*pt->x-b)/sqrt(1+A*A)-(pt->t-tau)*pt->drs*vd)*(pt->x+A*pt->y);
    double ds = ((A*pt->x+b-pt->y)/sqrt(1+A*A) + (tau-pt->t)*pt->drs*vd)*(pt->x+A*pt->y);
    sum += ds;
  }

  return sum;
}

//-----------------------------------------------------------------------------
double SegmentFit::DfDa(double A) {
  double vd  = SegmentHit::fgVDrift;

  double sum = 0;

  double b       = B(A);
  double db_da   = DbDa(A);
  double tau     = Tau(A);
  double dtau_da = DbDa(A);

  int npt = fSegment->nHits();
  for (int i=0; i<npt; i++) {
    SegmentHit* pt = fSegment->Hit(i);
    if (pt->IsGood() != 1) continue;
    int    si = pt->drs;
    //    double x1 = ((A*b-db_da*(A*A+1)-pt->x-A*pt->y)/pow(1+A*A,3/2.)+si*vd*dtau_da)*(pt->x+A*pt->y);
    double x1 = ((pt->x + A*pt->y - A*b + db_da*(A*A+1))/pow(1+A*A,3/2.) + si*vd*dtau_da)*(pt->x+A*pt->y);
    //    double x2 = ((pt->y-A*pt->x-b)/sqrt(1+A*A)-pt->drs*vd*(pt->t-tau))*pt->y;
    double x2 = ((A*pt->x+b-pt->y)/sqrt(1+A*A) + si*vd*(tau-pt->t))*pt->y;
    sum += x1+x2;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// call in the beginning, and them - only if the point content chages
// calculate the sums using only points with defined drift directions
//-----------------------------------------------------------------------------
int SegmentFit::CalculateLsqSums() {
  int rc(0);

  fSxs.clear();
  fSys.clear();
  fSts.clear();

  int nhits = fSegment->nHits();
  for (int i=0; i<nhits; i++) {
    SegmentHit* pt = fSegment->Hit(i);
    if (pt->IsGood() != 1)                      continue;
    if (pt->drs != 0) {
      fSxs.addPoint(pt->x,pt->drs);
      fSys.addPoint(pt->y,pt->drs);
      fSts.addPoint(pt->t,pt->drs);
    }
  }

  return rc;
}

//-----------------------------------------------------------------------------
int SegmentFit::DefineDriftDirections(const Par_t* Pin) {
//-----------------------------------------------------------------------------
// use two hits with known drift directions and also the first ant the last hits
//-----------------------------------------------------------------------------
  int sgn[4][2] = {{-1,-1}, {-1,1}, {1,-1}, {1,1}};
  int ifirst(-1), ilast(-1);

  if (SegmentFit::fgDebugMode != 0) {
    std::cout << std::format("-- SegmentFit::{}:{} START\n",__func__,__LINE__);
  }

  int nhits = fSegment->nHits();
  for (int i=0; i<nhits; ++i) {
    SegmentHit* pt  = fSegment->Hit(i);
    if (pt->IsGood() == 0) continue;
                                        // first good point
    ifirst = i;
    break;
  }

  for (int i=nhits-1; i>=0; --i) {
    SegmentHit* pt  = fSegment->Hit(i);
    if (pt->IsGood() == 0) continue;
                                        // last good point
    ilast = i;
    break;
  }

  SegmentHit* pt[2];
  pt[0] = fSegment->Hit(ifirst);
  pt[1] = fSegment->Hit(ilast );

  int    ibest(-1);
  double chi2_best(1.e12);
                                        // initial parameters
  const Par_t* pin = Pin;
  if (pin == nullptr) pin = &fSegment->fTangentLine;
  
                                        // parameters: a,b,T0
  Par_t   fit_res[4];
//-----------------------------------------------------------------------------
// at this moment, only four points have their drift signs defined
//-----------------------------------------------------------------------------
  for (int i=0; i<4; i++) {
    pt[0]->drs = sgn[i][0];
    pt[1]->drs = sgn[i][1];

    if (SegmentFit::fgDebugMode != 0) {
      std::cout << std::format("-- SegmentFit::{}:{} signs: {}:{}\n",__func__,__LINE__,pt[0]->drs,pt[1]->drs);
    }
                                        // make sure always start from the same parameter values
    Fit(1,0,pin,&fit_res[i]);
//-----------------------------------------------------------------------------
// four drift signs defined, fit , do not exclude the hits
//-----------------------------------------------------------------------------
    if (fit_res[i].chi2dof < chi2_best) {
      ibest     = i;
      chi2_best = fit_res[i].chi2dof;
    }
  }

  if (SegmentFit::fgDebugMode != 0) {
    std::cout << std::format("-- SegmentFit::{}:{} ibest:{} signs: {}:{} chi2_best:{:8.2f}\n",__func__,__LINE__,ibest,pt[0]->drs,pt[1]->drs,chi2_best);
  }

  Par_t* pbest = &fit_res[ibest];
//-----------------------------------------------------------------------------
// update segment line parameters
//-----------------------------------------------------------------------------
  fSegment->fPar4 = *pbest;
  fSegment->fPar  = *pbest;
//-----------------------------------------------------------------------------
// set drift directions of the first and the last hits
//-----------------------------------------------------------------------------
  pt[0]->drs = sgn[ibest][0];
  pt[1]->drs = sgn[ibest][1];
//-----------------------------------------------------------------------------
// define drift directions of the rest hits wrt the best parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; ++i) {
    SegmentHit* p  = fSegment->Hit(i);
    if (p->IsGood() == 0) continue;
                                        // skip four points with [supposedly] defined drift directions
    if (p->drs != 0) continue;
    double dist   = (pbest->X0()-p->x)*pbest->Nux()+(pbest->Y0()-p->y)*pbest->Nuy();
    double rdrift = fSegment->R(p,pbest->T0());
    double dr1    = dist-rdrift;
    double dr2    = dist+rdrift;
    double min_dist(1.e6);
    if (fabs(dr1) < fabs(dr2)) {
      p->drs   =  1;
      min_dist = fabs(dr1);
    }
    else {
      p->drs   = -1;
      min_dist = fabs(dr2);
    }
    if (min_dist > 1.) p->fMask |= TrkSegment::kLargeDrhoBit;
  }
  if (SegmentFit::fgDebugMode != 0) {
    fSegment->print();
    std::cout << std::format("-- SegmentFit::{}:{} END\n",__func__,__LINE__);
  }

  return 0;
}

//-----------------------------------------------------------------------------
// which parameters to use ?
//-----------------------------------------------------------------------------
int SegmentFit::DisplaySegment() {
//-----------------------------------------------------------------------------
// plot the results after the last iteration, but before updating the points
//-----------------------------------------------------------------------------
  if (fgDebugMode) std::cout << std::format("-- SegmentFit::{}:{} START\n",__func__,__LINE__);
      
  int nhits = fSegment->nHits();

  double xm = (fSegment->Hit(nhits-1)->x+fSegment->Hit(0)->x)/2;
  double dx = (fSegment->Hit(nhits-1)->x-fSegment->Hit(0)->x);
  if (dx < 10) dx = 10;

  // draw the line
  double xmin = xm-dx/2-10.;  // mm
  double xmax = xm+dx/2+10.;
  double ymin = fSegment->DyDx()*xmin+fSegment->Y0();
  double ymax = fSegment->DyDx()*xmax+fSegment->Y0();
  TLine* l = new TLine(xmin,ymin,xmax,ymax);
  l->SetLineColor(kRed+2);
  l->Draw();

  // draw hits

  double t0 = fSegment->T0();
  for (int i=0; i<nhits; i++) {
    SegmentHit* pt = fSegment->Hit(i);
    if (fgDebugMode and fgDebugBits[0]) pt->print();
    // straw
    TEllipse* e = new TEllipse(pt->x,pt->y,2.5,0,0,360);
    e->SetLineColor(kRed);
    e->SetFillStyle(0);
    e->Draw();
    // hit
    double r = fSegment->R(pt,t0);
      TEllipse* e2 = new TEllipse(pt->x,pt->y,r,0,0,360);
    if (pt->fMask == 0) {
      e2->SetFillColor(kBlue-10);
      e2->SetFillStyle(3003);
    }
    else {
      e2->SetFillColor(kRed-10);
      e2->SetFillStyle(3001);
    }
    e2->Draw();
  }
  if (fgDebugMode) std::cout << std::format("-- SegmentFit::{}:{} END\n",__func__,__LINE__);
  return 0;
}

//-----------------------------------------------------------------------------
// rc = 1: fit converged ???
// assume that the initial tangent line is defined
//-----------------------------------------------------------------------------
int SegmentFit::Fit(int NItMax, int DoCleanup, const Par_t* Pin, Par_t* Par) {
  int  converged(0);

  if (fgDebugMode) std::cout << std::format("-- SegmentFit::{}:{} START\n",__func__,__LINE__);

                                        // pin: initial parameter values
  const Par_t* pin = Pin;
  if (pin == nullptr) pin = &fSegment->fPar4;

  double a0   = pin->a;
  double tau0 = pin->tau;

  for (int iter=0; iter<NItMax; iter++) {
    if (fgDebugMode) {
      std::cout << std::format("SegmentFit::{} -- fit iteration:{} a0:{:10.5}\n", __func__,iter,a0);
    }
                                        // all drift directions have to be defined
    CalculateLsqSums();

    double f0    = F(a0);
    double dfda0 = DfDa(a0);
                                        // parameters of the line at the next iteration
    double a1    = a0-f0/dfda0;
    double b1    = B  (a1);
    double tau   = Tau(a1);

    if (fgDebugMode) {
      std::cout << std::format(" f0:{:10.3f} dfda0:{:10.5f} a1:{:10.5f} b1:{} tau:{:10.5f}\n",
                               f0,dfda0,a1,b1,tau);
    }
//-----------------------------------------------------------------------------
// update parameters, recalculate the chi2
//-----------------------------------------------------------------------------
    fSegment->UpdateParameters(a1,b1,tau);
//-----------------------------------------------------------------------------
// the rest is an attempt to remove bad hits if any, ignore for now
//-----------------------------------------------------------------------------
//     double x0    = 0;
//     double y0    = b1;
//     double ny    = a1/sqrt(1+a1*a1); // keep nx > 0;
//     double nx    =  1/sqrt(1+a1*a1);
//     double nux   = -ny;
//     double nuy   =  nx;

//     std::cout << std::format("-- {} iteration:{:2} a0:{:12.7f}  a1:{:12.7f} b1:{:12.7f} tau:{:10.3f}\n",
//                              __func__,iter,a0,a1,b1,tau);

//                                         // at this point need to recalculate the radii
//     double   chi2_tot   = 0;
//     double   chi2_dof   = 0;
//     double   drho       = 0;
//     double   chi2_worst (-1);
//     double   rho_worst(0);
//     Point2D* worst_hit(nullptr);

//     int nhits = fSegment->points.size();
//     int n_good_hits = 0;
//     if (nhits > 3) {
//       if (fgDebugMode) {
//         std::cout << std::format(" i     xi       yi      ri       si      rho     chi2\n");
//         std::cout << std::format("-----------------------------------------------------\n");
//       }
//       for (int i=0; i<nhits; i++) {
//         Point2D* pt = &fSegment->points[i];
//         if (pt->IsGood() != 1)                           continue;
//         n_good_hits += 1;
//         int    si  = pt->drs;
//         double ri  = fSegment->R(i)*si;
//         double xi  = pt->x+nux*ri;
//         double yi  = pt->y+nuy*ri;
//         double rho = ((x0-xi)*nux+(y0-yi)*nuy)*si;
//         double chi2_hit = rho*rho/(0.2*0.2);
//         if (fgDebugMode) {
//           std::cout << std::format("{:2} {:10.4f} {:10.4f} {:10.4f} {:2} {:10.4f} {:10.3f}",
//                                    i,xi,yi,ri,si,rho,chi2_hit)
//                     << std::endl;
//         }
//         drho     += rho;
//         chi2_tot += chi2_hit;

//         if (chi2_hit > chi2_worst) {
//           worst_hit  = pt;
//           chi2_worst = chi2_hit;
//           rho_worst  = rho;
//         }
//       }

//       chi2_dof = chi2_tot/(n_good_hits-2.99999);
//       drho    /= n_good_hits;

//       Par->chi2dof  = chi2_dof;

//       if (chi2_dof > 5.) {
//         if (DoCleanup) {
// //-----------------------------------------------------------------------------
// // check if the chi2 is dominated by a single hit and remove it if needed
// // one could do it better, but should be good enough for the purpose
// // in case of a cleanup, do we need to recalculate the fit parameters ?
// //-----------------------------------------------------------------------------
//           if (chi2_tot - chi2_worst < n_good_hits) {
//             worst_hit->fMask = 0x1;
//             chi2_tot    -= chi2_worst;
//             drho        -= rho_worst;
//             n_good_hits -= 1;
//                                         // and recalculate the chi2
//             chi2_dof     = chi2_tot/(n_good_hits-2.99999);
//             drho        /= n_good_hits;
//           }
//         }
//       }
//       std::cout << std::format("SegmentFit::{}:{} chi2:{:8.3f} drho:{:8.5f}",__func__,__LINE__,chi2_dof,drho) << std::endl;
//     }
    if (fabs(tau - tau0) < 0.1) {
      converged = 1;
      break;
    }
//-----------------------------------------------------------------------------
// continue iterations
//-----------------------------------------------------------------------------
    a0 = a1;
  }
                                        // return best parameters
  *Par = fSegment->fPar;

  if (fgDebugMode)  std::cout << std::format("-- SegmentFit::{} END  Par->chi2_dof:{:8.2f} converged:{}\n",__func__,Par->chi2dof,converged);

  return converged;
}

//-----------------------------------------------------------------------------
double SegmentFit::DChi2Db  (double A, double B, double Tau) {
  double sum = 0;
  int nhits = fSegment->nHits();
  for (int i=0; i<nhits; ++i) {
    SegmentHit* p  = fSegment->Hit(i);
    if (not p->IsGood())                                    continue;
                                        // do we need to skip bad hits ? - I think, we do! 
    double tdrift = p->t+fSegment->fTMean-Tau;
    double dist   = (A*p->x + B - p->y)/sqrt(1+A*A) - p->drs*tdrift*SegmentHit::fgVDrift;
    sum += dist;
  }

  return sum;
}

//-----------------------------------------------------------------------------
double SegmentFit::DChi2Dtau  (double A, double B, double Tau) {
  double sum = 0;
  int nhits = fSegment->nHits();
  for (int i=0; i<nhits; ++i) {
    SegmentHit* p = fSegment->Hit(i);
    if (not p->IsGood())                                    continue; 
    double tdrift = p->t+fSegment->fTMean-Tau;
    double dist   = (A*p->x + B - p->y)/sqrt(1+A*A) - p->drs*tdrift*SegmentHit::fgVDrift;
    sum          += dist*p->drs;
  }

  return sum;
}

//-----------------------------------------------------------------------------
double SegmentFit::DChi2Da  (double A, double B, double Tau) {
  double sum = 0;
  int nhits = fSegment->nHits();
  for (int i=0; i<nhits; ++i) {
    SegmentHit* p = fSegment->Hit(i);
    if (not p->IsGood())                                    continue; 
    double tdrift = p->t+fSegment->fTMean-Tau;
    double rdrift = tdrift*SegmentHit::fgVDrift;
    double d1     = (p->y-A*p->x-B)/sqrt(1+A*A);
    double dist   = d1+rdrift*p->drs;
    printf("i:%2i  d1:%12.5f dist:%12.5f  rdrift:%12.5f\n",i,d1,dist,rdrift);
    double d2     = p->x*(1+A*A)+A*(p->y-A*p->x-B);
    sum          += dist*d2;
  }

  return sum;
}
