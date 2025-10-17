//
#include <format>
#include <ostream>
#include "daqana/obj/TrkSegment.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"

int TrkSegment::_debugMode(0);

//
double const Point2D::fgVDrift    = 65.e-3;   // mm/ns
double const TrkSegment::fgRStraw =  2.5;     // mm
//-----------------------------------------------------------------------------
TrkSegment::TrkSegment(int Plane, int Panel) {
  plane       = Plane;
  panel       = Panel;
  trkPanel    = nullptr;
  fCombiTrans = nullptr;
  //  fChi2Dof    = -1;
  fMask       = 0;
  fIhit[0]    = -1;
  fIhit[1]    = -1;
}

//-----------------------------------------------------------------------------
int TrkSegment::InitHits() {
  int rc(0);
  std::cout << "ERROR: TrkSegment::InitHits() yet to be written! \n" ;
  return rc;
}

//-----------------------------------------------------------------------------
int TrkSegment::InitHits(std::vector<ComboHitData_t>& Hits, int UniquePlane, int Panel) {
  int rc(0);

  int nhits = Hits.size();
  std::cout << std::format("TrkSegment::{}: BEGIN nhits:{}\n",__func__,nhits);

  //  int all_hits(1);
  // if ((UniquePlane != -1) && (Panel != -1)) all_hits = 0;

//   fT0           = 0;
//   fNGoodHits    = 0;
//   for (int i=0; i<nhits; i++) {
//     ComboHitData_t* hit = &Hits[i];
//     if (hit->IsGood() != 1)              continue;
//     fT0        += hit->tcorr;
//     fNGoodHits += 1;
//   }
//                                         // in principle, this approximation for T0 should be good enough
//   fT0 = fT0/(fNGoodHits+1.e-12);
//   if (fNGoodHits < 4) {
//     fMask = kNHitsBit;
//   }
// //-----------------------------------------------------------------------------
// // make sure initial drift radii do not go above the straw radius - within the assumptioins used,
// // R > Rstraw means too small value of T0
// //-----------------------------------------------------------------------------
//   double rstraw(2.5);
//   double rmax(rstraw);
//   for (int i=0; i<nhits; i++) {
//     ComboHitData_t* hit = &Hits[i];
//     double td = hit->time-hit->prtime-fT0;
//     double r  = td*Point2D::fgVDrift;
//     if (r > rstraw) rmax = r;
//   }
//   if (rmax > rstraw) fT0 = fT0+(rmax-rstraw)/Point2D::fgVDrift;
//-----------------------------------------------------------------------------
// find seed hits, ignore the very first and the very last pairs
//-----------------------------------------------------------------------------
  fNTransitions = 0;
  ComboHitData_t* last_hit = nullptr;

  for (int i=0; i<nhits; i++) {
    ComboHitData_t* hit = &Hits[i];
                                        // figure out the drift signs
    int layer = mu2e::StrawId(hit->sid).getLayer();
    if (last_hit != nullptr) {
      int last_hit_layer = mu2e::StrawId(last_hit->sid).getLayer();
                                        // layer is either 0 or 1
      if (last_hit_layer != layer) {
//-----------------------------------------------------------------------------
// transition between the layers found.
// ignore back-and-forth transitions - instead, mark hit as bad
//-----------------------------------------------------------------------------
        if (i < nhits-1) {
          ComboHitData_t* next_hit = &Hits[i+1];
          int next_hit_layer = mu2e::StrawId(next_hit->sid).getLayer();
          if (next_hit_layer == last_hit_layer) {
            hit->flag = 0x1;
                                        // start from scratch
            last_hit  = nullptr;
            continue;
          }
        }
//-----------------------------------------------------------------------------
// also ignore transitions involving the first and the last hits
//-----------------------------------------------------------------------------
        if ((i == 1) or (i == nhits-1)) {
          last_hit = hit;
          continue;
        }
//-----------------------------------------------------------------------------
// get here only if hits is initially good (not 'back-n-forth')
// rotations don't change the Z-coordinate
//-----------------------------------------------------------------------------
        fNTransitions += 1;
                                        // z becomes Y in points
        if (last_hit->z > hit->z) {
          last_hit->drs      = -1;
          hit->drs           =  1;
        }
        else if (last_hit->z < hit->z) {
          last_hit->drs      =  1;
          hit->drs           = -1;
        }
        fIhit[0] = i-1;
        fIhit[1] = i;
      }
    }
    last_hit    = hit;
  }
  if (fNTransitions != 1) fMask |= kNTransitionsBit;
//-----------------------------------------------------------------------------
// add Point2D's :
//-----------------------------------------------------------------------------
  //  int    imax      (-1);
  //  double max_drtime(-1);

  fXMean = 0;
  fYMean = 0;
  fTMean = 0;

  if (_debugMode != 0) {
    std::cout << std::format(" i   sid                  XYZ(M)                              XYZ(L)                 time    prtime  rdrift drs\n");
    std::cout << std::format("--------------------------------------------------------------------------------------\n");
  }

  fNGoodHits = 0;
  for (int i=0; i<nhits; i++) {
    ComboHitData_t* ch = &Hits[i];
    if (ch->IsGood() == 0) continue;
    fNGoodHits += 1;

    double posm[3], posl[3];
    posm[0] = ch->x;
    posm[1] = ch->y;
    posm[2] = ch->z;
    fCombiTrans->MasterToLocalVect(posm, posl);
//-----------------------------------------------------------------------------
// points have coordinates in the local coordinate system of the panel
//-----------------------------------------------------------------------------
    addPoint(ch->sid,ch->flag, posl[0],posl[2],ch->time,ch->prtime, ch->tcorr, ch->drs);
//-----------------------------------------------------------------------------
// to avoid mistakes, perform the subtraction just once
// hits flagged at this stage will never be unflagged
// and don't assume that <X>, <Y> etc are exactly equal to zero
//-----------------------------------------------------------------------------
    fXMean += posl[0];
    fYMean += posl[2];
    fTMean += ch->time-ch->prtime;
    // double drtime = ch->time-ch->prtime-fT0;
    // if (drtime > max_drtime) {
    //   // imax       = i;
    //   max_drtime = drtime;
    // }
    if (_debugMode != 0) {
      std::cout << std::format("{:2} 0x{:05x} ({:10.3f} {:10.3f} {:10.3f})",i,ch->flags,ch->x, ch->y, ch->z)
                << std::format(" ({:10.3f} {:10.3f} {:10.3f}) {:10.3f} {:7.3f} {:2}",
                               posl[0], posl[1], posl[2], ch->time, ch->prtime, ch->drs)
                << std::endl;
    }
  }

  fXMean = fXMean/fNGoodHits;
  fYMean = fYMean/fNGoodHits;
  fTMean = fTMean/fNGoodHits;
//-----------------------------------------------------------------------------
// it is enough to know that there is no large offsets which could affect the
// calculations; x,y, and t to be used in calculations, their values should be close to zero
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; i++) {
    Point2D* p  = &points[i];
    p->x       -= fXMean;
    p->y       -= fYMean;
    p->t        = p->time-p->tprop-fTMean;
  }
//-----------------------------------------------------------------------------
// finally, sort points and define tangent line
//-----------------------------------------------------------------------------
  std::sort(points.begin(),points.end(),
            [] (const Point2D& a, const Point2D& b) {
              return a.strawId().straw() < b.strawId().straw();
            });

  DefineTangentLine();

  if (_debugMode != 0) {
    print("InitHits after DefineTangentLine");
    std::cout << std::format("-- TrkSegment.cc:{} END rc:{}\n",__LINE__,rc);
  }
  return rc;
}

//-----------------------------------------------------------------------------
void TrkSegment::print(const char* Message, std::ostream& Stream) {

  if (Message) Stream << Message << std::endl;

  int nhits = points.size();

  Stream << std::format(" nhits:{} fNGoodHits:{} n_transitions:{} transition hits: {}:{} ",
                        nhits,fNGoodHits,fNTransitions,fIhit[0],fIhit[1])
         << std::format(" Parameters: DyDx:{:10.5f} Y0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f} T0:{:10.3f} chi2/DOF:{:8.3f}",
                        fPar.a,fPar.b,fPar.nx,fPar.ny,fPar.tau+fTMean,fPar.chi2dof)
         << std::endl;

  Stream << std::format(" i   sid   mask       x      y        t        time   tprop tdrift   radius drs      rho  fixed\n");
  Stream <<             "-----------------------------------------------------------------------------------------------\n" ;
  double t0 = fPar.T0();
  for (int i=0; i<nhits; i++) {
    Point2D* pt = &points[i];
    Stream << std::format("{:2d} 0x{:04x} 0x{:04x} {:8.3f} {:6.3f} {:7.3f} {:10.3f} {:7.3f} {:7.3f} {:7.3f} {:2d} {:10.4f} {:5d}",
                          i,pt->sid,pt->fMask,pt->x,pt->y,pt->t,pt->time,pt->tprop,
                          DriftTime(pt,t0),R(pt,t0),pt->drs,Rho(pt,&fPar),pt->sign_fixed)
           << std::endl;
  }
}


//------------------------------------------------------------------------------
// in general, the worst hit could be not the one which has the largest chi2
//-----------------------------------------------------------------------------
double TrkSegment::Chi2() {
  double const res2(0.2*0.2);            // assume 200 um

  double chi2 = 0;
  int nhits = points.size();
  int n_good_hits(0);
  double t0 = T0();
  double a  = DyDx();
  double b  = Y0();
  for (int i=0; i<nhits; ++i) {
    Point2D* pt = &points[i];
    if (pt->IsGood() == 0) continue;
    double rd   = R(pt,t0);
    if (rd < 0) rd = 0;
    if (rd > fgRStraw) rd = fgRStraw;
    double dist = (a*pt->x + b - pt->y)/sqrt(1+a*a)-rd*pt->drs;
    chi2        += dist*dist/res2;
    n_good_hits += 1;
  }

  chi2 = chi2/(n_good_hits - 2.999999);
  return chi2;
}

//-----------------------------------------------------------------------------
void  TrkSegment::UpdateParameters(double DyDx, double Y0, double T0, double Chi2Dof) {
  fPar.a   = DyDx;
  fPar.b   = Y0;
  fPar.nx  =  1.  /sqrt(1+DyDx*DyDx);
  fPar.ny  =  DyDx/sqrt(1+DyDx*DyDx);
  fPar.tau = T0;

  if (Chi2Dof >= 0) {
    fPar.chi2dof = Chi2Dof;
  }
  else {
    fPar.chi2dof = Chi2();
  }
}


//-----------------------------------------------------------------------------
// find line tangent to two circles defined by the two hits with known drift signs
// the segment fT0 is expected to be defined
//-----------------------------------------------------------------------------
int TrkSegment::DefineTangentLine() {
  int rc(0);
  if (_debugMode) std::cout << std::format("-- TrkSegment::{}: BEGIN\n",__func__);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  int i1       = fIhit[0];
  int i2       = fIhit[1];

  Point2D* p1  = &points[i1];
  Point2D* p2  = &points[i2];
//-----------------------------------------------------------------------------
// define T0
//-----------------------------------------------------------------------------
  double t0(0);
  int    n_good_hits(0);

  int nhits = points.size();
  for (int i=0; i<nhits; i++) {
    Point2D* pt = &points[i];
    if (pt->IsGood() != 1)              continue;
    t0          += pt->tcorr;
    n_good_hits += 1;
  }
                                         // in principle, this approximation for T0 should be good enough
  t0 = t0/(n_good_hits+1.e-12) - fTMean;
//-----------------------------------------------------------------------------
// make sure t0 is not larger than the smallest measured time
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; i++) {
    Point2D* pt = &points[i];
    if (pt->IsGood() != 1)              continue;
    if (t0 > pt->t) t0 = pt->t;
  }


  double r1    = R(p1,t0);
  double r2    = R(p2,t0);

  double dx    = p2->x-p1->x;
  double dy    = p2->y-p1->y;
  double alp   = -p2->drs*r2+p1->drs*r1;
  
  double dr    = sqrt(dx*dx+dy*dy);
  double alpdr = alp/dr;
  double dxr   = dx/dr;    
  double dyr   = dy/dr;    
  double nx    = -alpdr*dyr-dxr*sqrt(1-alpdr*alpdr);
  double ny    =  alpdr*dxr-dyr*sqrt(1-alpdr*alpdr);

  if (nx < 0) {
                                        // force nx > 0
    nx = -nx;
    ny = -ny;
  }
                                        // direction of the normal to the line: rotated by +90 deg
  double nux   = -ny;
  double nuy   =  nx;

  double x1    = p1->x+p1->drs*r1*nux;
  double y1    = p1->y+p1->drs*r1*nuy;

  double t     = -x1/nx;
                                        // parameters of the line, don't update T0
  fTangentLine.a       = ny/nx;
  fTangentLine.b       = y1 + t*ny;  // Y0
  fTangentLine.nx      = nx;
  fTangentLine.ny      = ny;
  fTangentLine.tau     = t0;
  fTangentLine.chi2dof = 0;
                                        // at the moment, this is the best...
  fPar = fTangentLine;

  if (_debugMode) {
    std::cout << std::format("-- TrkSegment::{}: END rc:{} fDyDx:{:10.5f} fY0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f}\n",
                             __func__,rc,fTangentLine.a,fTangentLine.b,fTangentLine.nx,fTangentLine.ny);
  }
  return rc;
}
