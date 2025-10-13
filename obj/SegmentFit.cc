///////////////////////////////////////////////////////////////////////////////
#include "TLine.h"
#include "TEllipse.h"

#include "daqana/obj/SegmentFit.hh"
#include <memory>


int SegmentFit::fgDebugMode(0);

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
int SegmentFit::FindTangentLine(TrkSegment* Seg) {
  int rc(0);

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  // Point2D* p1  = &points[fIhit[0]];
  // Point2D* p2  = &points[fIhit[1]];

  // double dx    = p2->x-p1->x;
  // double dy    = p2->y-p1->y;
  // double alp   = -p2->drs*p2->r()+p1->drs*p1->r();
  
  // double dr    = sqrt(dx*dx+dy*dy);
  // double alpdr = alp/dr;
  // double dxr   = dx/dr;    
  // double dyr   = dy/dr;    
  // double nx    = -alpdr*dyr-dxr*sqrt(1-alpdr*alpdr);
  // double ny    =  alpdr*dxr-dyr*sqrt(1-alpdr*alpdr);

  return rc;
}

//-----------------------------------------------------------------------------
// find line tangent to two circles defined by the two hits with known drift signs
//-----------------------------------------------------------------------------
double SegmentFit::F(double A) {
  // double smn = fSts.yMean();
  double vd  = Point2D::fgVDrift;

  int npt = fSegment->points.size();

  double sum = 0;

  double b   = B(A);
  double tau = Tau(A);

  for (int i=0; i<npt; i++) {
    Point2D* pt = &fSegment->points[i];
    if (pt->IsGood() != 1) continue;

    double ds = ((pt->y-A*pt->x-b)/sqrt(1+A*A)-(pt->t-tau)*pt->drs*vd)*(pt->x+A*pt->y);
    sum += ds;
  }

  return sum;
}

//-----------------------------------------------------------------------------
double SegmentFit::DfDa(double A) {
  // double smn = fSts.yMean();
  double vd  = Point2D::fgVDrift;

  double sum = 0;

  double b       = B(A);
  double db_da   = DbDa(A);
  double tau     = Tau(A);
  double dtau_da = DbDa(A);

  int npt = fSegment->points.size();
  for (int i=0; i<npt; i++) {
    Point2D* pt = &fSegment->points[i];
    if (pt->IsGood() != 1) continue;

    double x1 = ((A*b-db_da*(A*A+1)-pt->x-A*pt->y)/pow(1+A*A,3/2.)+pt->drs*vd*dtau_da)*(pt->x+A*pt->y);
    double x2 = ((pt->y-A*pt->x-b)/sqrt(1+A*A)-pt->drs*vd*(pt->t-tau))*pt->y;
    sum += x1+x2;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// call in the beginning, and them - only if the point content chages
// calculate the sums using only points with defined drift directions
//-----------------------------------------------------------------------------
int SegmentFit::DefineDriftDirections() {
//-----------------------------------------------------------------------------
// use two hits with known drift directions and also the first ant the last hits
//-----------------------------------------------------------------------------
  int sgn[4][2] = {{-1,-1}, {-1,1}, {1,-1}, {1,1}};
  int ifirst(-1), ilast(-1);

  int nhits = fSegment->points.size();
  for (int i=0; i<nhits; ++i) {
    Point2D* pt  = &fSegment->points[i];
    if (pt->IsGood() == 0) continue;
                                        // first good point
    ifirst = i;
    break;
  }

  for (int i=nhits-1; i>=0; --i) {
    Point2D* pt  = &fSegment->points[i];
    if (pt->IsGood() == 0) continue;
                                        // last good point
    ilast = i;
    break;
  }

  Point2D* p0 = &fSegment->points[ifirst];
  Point2D* p1 = &fSegment->points[ilast ];

  int    ibest(-1);
  double chi2_best(1.e12);
                                        // parameters: a,b,T0
  Par_t   fit_res[4];
//-----------------------------------------------------------------------------
// at this moment, only four points have their drift signs defined
//-----------------------------------------------------------------------------
  for (int i=0; i<4; i++) {
    p0->drs = sgn[i][0];
    p1->drs = sgn[i][1];
                                        // four drift signs defined
    Fit(1,&fit_res[i]);
    if (fit_res[i].chi2dof < chi2_best) {
      ibest     = i;
      chi2_best = fit_res[i].chi2dof;
    }
  }
//-----------------------------------------------------------------------------
// choose the fit returning the best chi2
//-----------------------------------------------------------------------------
  double a  = fit_res[ibest].a;
  // double b  = fit_res[ibest].b;
  // double t0 = fit_res[ibest].tau;

  fSegment->fNx      = 1/sqrt(a*a+1);                  // choose positive sign
  fSegment->fNy      = a/sqrt(a*a+1);
  fSegment->fDyDx    = a;
  double t = -(p1->x+p1->drs*p1->r()*fSegment->Nux())/fSegment->Nx();
  fSegment->fY0      = p1->y+p1->drs*p1->r()*fSegment->Nuy()+fSegment->Ny()*t;
  if (fgDebugMode != 0) {
    // double n2 = sqrt(nx*nx+ny*ny);
    std::cout << std::format("-- TrkSegment::{}:{} nx,ny, n2 = : {:12.5e} {:12.5e}\n",
                             __func__,__LINE__,fSegment->fNx,fSegment->fNy);
  }
//-----------------------------------------------------------------------------
// after that, define drift directions of all other hits
//-----------------------------------------------------------------------------
  p0->drs = sgn[ibest][0];
  p1->drs = sgn[ibest][1];

  // int nhits = fSegment->points.size();
  for (int i=0; i<nhits; ++i) {
    Point2D* p  = &fSegment->points[i];
    if (p->IsGood() == 0) continue;
                                        // skip four points with [supposedly] defined drift directions
    if (p->drs != 0) continue;
    double dr1 = (0-p->x)*fSegment->Nux()+(fSegment->fY0-p->y)*fSegment->Nuy()+p->r();
    double dr2 = (0-p->x)*fSegment->Nux()+(fSegment->fY0-p->y)*fSegment->Nuy()-p->r();
    if (fabs(dr1) < fabs(dr2)) p->drs =  1;
    else                       p->drs = -1;
  }

  return 0;
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

  int nhits = fSegment->points.size();
  for (int i=0; i<nhits; i++) {
    Point2D* pt = &fSegment->points[i];
    pt->print();
    if (pt->IsGood() != 1)                      continue;
    if (pt->drs != 0) {
      fSxs.addPoint(pt->x,pt->drs);
      fSys.addPoint(pt->x,pt->drs);
      fSts.addPoint(pt->t,pt->drs);
    }
  }

  return rc;
}

//-----------------------------------------------------------------------------
// rc = 0: fit converged ???
//-----------------------------------------------------------------------------
int SegmentFit::Fit(int NItMax, Par_t* Par) {
  int rc(0);
  bool converged(false);
  for (int iter=0; iter<NItMax; iter++) {
    std::cout << std::format("-- SegmentFit::{}:{} ----------------------------------- iteration {}\n",
                             __func__,__LINE__,iter);

                                        // 1. update drift radii - is that needed ? - perhaps not, the hit radii are calculable
    CalculateLsqSums();

    double a0    = fSegment->DyDx();
    double f0    = F(a0);
    double dfda0 = DfDa(a0);
                                        // parameters of the line at the next iteration
    double a1    = a0-f0/dfda0;
    double b1    = B  (a1);
    double tau   = Tau(a1);
//-----------------------------------------------------------------------------
// update parameters, recalculate the chi2
//-----------------------------------------------------------------------------
    fSegment->UpdateLineParameters(a1,b1,tau);
    fSegment->CalculateChi2();

    double x0    = 0;
    double y0    = b1;
    double ny    = a1/sqrt(1+a1*a1); // keep nx > 0;
    double nx    =  1/sqrt(1+a1*a1);
    double nux   = -ny;
    double nuy   =  nx;

    std::cout << std::format("-- {} iteration:{:2} a0:{:12.7f}  a1:{:12.7f} b1:{:12.7f}",
                             __func__,iter,a0,a1,b1)
              << std::endl;
                                        // at this point need to recalculate the radii
    double   chi2_tot   = 0;
    double   chi2_dof   = 0;
    double   drho       = 0;
    double   chi2_worst (-1);
    double   rho_worst(0);
    Point2D* worst(nullptr);

    int nhits = fSegment->points.size();
    int n_good_hits = 0;
    if (nhits > 3) {
      for (int i=0; i<nhits; i++) {
        Point2D* pt = &fSegment->points[i];
        if (pt->IsGood() != 1)                           continue;
        n_good_hits += 1;
        int    si  = pt->drs;
        double ri  = pt->r()*si;
        double xi  = pt->x+nux*ri;
        double yi  = pt->y+nuy*ri;
        double rho = ((x0-xi)*nux+(y0-yi)*nuy)*si;
        std::cout << std::format("i:{:2} xi,yi,ri,si,rho: {:10.4f} {:10.4f} {:10.4f} {:2} {:10.4f}",
                                 i,xi,yi,ri,si,rho)
                  << std::endl;
        
        drho     += rho;
        double chi2_hit = rho*rho/(0.2*0.2);
        chi2_tot += chi2_hit;

        if (chi2_hit > chi2_worst) {
          worst      = pt;
          chi2_worst = chi2_hit;
          rho_worst  = rho;
        }
      }

      chi2_dof = chi2_tot/(n_good_hits-2.99999);
      drho    /= n_good_hits;

      if (chi2_dof > 5.) {
                                        // check if everything is coming from one hit
        if (chi2_tot - chi2_worst < n_good_hits) {
          worst->fMask = 0x1;
          chi2_tot    -= chi2_worst;
          drho        -= rho_worst;
          n_good_hits -= 1;
                                        // and recalculate the chi2
          chi2_dof     = chi2_tot/(n_good_hits-2.99999);
          drho        /= n_good_hits;
        }
      }
    }

    std::cout << std::format("{}:{} chi2:{:8.3f} drho:{:8.5f}",__func__,__LINE__,chi2_dof,drho) << std::endl;

    if (fabs(drho) < 0.001) {
      converged = true;
    }

    if ((iter == NItMax-1) or converged) { // one micron 
      DisplaySegment();
    }

    if (converged) break;
  }
  
  return rc;
}


//-----------------------------------------------------------------------------
int SegmentFit::DisplaySegment() {
//-----------------------------------------------------------------------------
// plot the results after the last iteration, but before updating the points
//-----------------------------------------------------------------------------
  // std::cout << std::format(" plotting: Seg->drho:{:10.5f}",Seg->drho)
  //           << std::endl;
      
  int nhits = fSegment->nHits();

  double xm = (fSegment->points[nhits-1].x+fSegment->points[0].x)/2;
  double dx = (fSegment->points[nhits-1].x-fSegment->points[0].x);
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
  for (int i=0; i<nhits; i++) {
    Point2D* pt = &fSegment->points[i];
    pt->print();
    // straw
    TEllipse* e = new TEllipse(pt->x,pt->y,2.5,0,0,360);
    e->SetLineColor(kRed);
    e->SetFillStyle(0);
    e->Draw();
    // hit
    TEllipse* e2 = new TEllipse(pt->x,pt->y,fSegment->R(i),0,0,360);
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
  return 0;
}

