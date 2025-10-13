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
  //  fDRho       = 0;
  fChi2Dof    = -1;
  fMask       = 0;
  fIhit[0]    = -1;
  fIhit[1]    = -1;
}

//-----------------------------------------------------------------------------
int TrkSegment::InitHits(std::vector<ComboHitData_t>& Hits, int UniquePlane, int Panel) {
  int rc(0);

  int nhits = Hits.size();
  std::cout << std::format("TrkSegment::{}: BEGIN nhits:{}\n",__func__,nhits);

  //  int all_hits(1);
  // if ((UniquePlane != -1) && (Panel != -1)) all_hits = 0;

  fT0           = 0;
  fNGoodHits    = 0;
  for (int i=0; i<nhits; i++) {
    ComboHitData_t* hit = &Hits[i];
    if (hit->IsGood() != 1)              continue;
    fT0        += hit->tcorr;
    fNGoodHits += 1;
  }
                                        // in principle, this approximation for T0 should be good enough
  fT0 = fT0/(fNGoodHits+1.e-12);
  if (fNGoodHits < 4) {
    fMask = kNHitsBit;
  }
//-----------------------------------------------------------------------------
// make sure initial drift radii do not go above the straw radius - that reflects
// too small value of T0
//-----------------------------------------------------------------------------
  double rstraw(2.5);
  double drmax(0);
  for (int i=0; i<nhits; i++) {
    ComboHitData_t* hit = &Hits[i];
    double td = hit->time-fT0-hit->prtime;
    double dr = td*Point2D::fgVDrift-rstraw;
    if (dr > drmax) drmax = dr;
  }
  if (drmax > 0) fT0 = fT0+drmax/Point2D::fgVDrift;
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
// ignore transitions involving the first and the last hits
//-----------------------------------------------------------------------------
        if ((i == 1) or (i == nhits-1))                 continue;
//-----------------------------------------------------------------------------
// ignore back-and-forth transitions - instead, mark hit as bad
//-----------------------------------------------------------------------------
        if (i < nhits-1) {
          ComboHitData_t* next_hit = &Hits[i+1];
          int next_hit_layer = mu2e::StrawId(next_hit->sid).getLayer();
          if (next_hit_layer == last_hit_layer) {
            hit->flag = 0x1;
            continue;
          }
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
  std::cout << std::format("-- TrkSegment::{}:{}: fT0:{:10.3f} fNTransitions:{} fNGoodHits:{}\n",
                           __func__,__LINE__,fT0,fNTransitions,fNGoodHits);

  //  int    imax      (-1);
  double max_drtime(-1);

  fXMean = 0;
  fYMean = 0;
  fTMean = 0;
  
  for (int i=0; i<nhits; i++) {
    ComboHitData_t* hit = &Hits[i];
    if (hit->IsGood() == 0) continue;

    double posm[3], posl[3];
    posm[0] = hit->x;
    posm[1] = hit->y;
    posm[2] = hit->z;
    fCombiTrans->MasterToLocalVect(posm, posl);

    if (_debugMode != 0) {
      std::cout << std::format("i:{} sid:0x{:05x} xyz(M):({:10.3f} {:10.3f} {:10.3f})", i, hit->sid, hit->x, hit->y, hit->z)
                << std::format(" xyz(L):({:10.3f} {:10.3f} {:10.3f}) time:{:10.3f} prtime:{:10.3f} r:{:10.3f} drs:{:2d}\n",
                               posl[0], posl[1], posl[2], hit->time, hit->prtime, hit->r,hit->drs);
    }
//-----------------------------------------------------------------------------
// points have coordinates in the local coordinate system of the panel
//-----------------------------------------------------------------------------
    addPoint(hit->sid,hit->flag, posl[0],posl[2],hit->time,hit->prtime, hit->drs);
//-----------------------------------------------------------------------------
// to avoid mistakes, perform the subtraction just once
// hits flagged at this stage will never be unflagged
// and don't assume that <X>, <Y> etc are exactly equal to zero
//-----------------------------------------------------------------------------
    fXMean += posl[0];
    fYMean += posl[2];
    fTMean += hit->time-hit->prtime;
    double drtime = hit->time-hit->prtime-fT0;
    if (drtime > max_drtime) {
      // imax       = i;
      max_drtime = drtime;
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
// sort points
//-----------------------------------------------------------------------------
  std::sort(points.begin(),points.end(),
            [] (const Point2D& a, const Point2D& b) {
              return a.strawId().straw() < b.strawId().straw();
            });

  if (_debugMode != 0) {
    std::cout << std::format("-- TrkSegment.cc:{}",__LINE__) << std::endl;
    print();
  }

  std::cout << std::format("TrkSegment::{}: END rc={}\n",__func__,rc);
  return rc;
}

//------------------------------------------------------------------------------
// in general, the worst hit could be not the one which has the largest chi2
//-----------------------------------------------------------------------------
int TrkSegment::CalculateChi2() {
  double const res2(0.2*0.2);            // assume 200 um

  double chi2 = 0;
  int nhits = points.size();
  int n_good_hits(0);
  for (int i=0; i<nhits; ++i) {
    Point2D* pt = &points[i];
    if (pt->IsGood() == 0) continue;
    double rd   = Point2D::fgVDrift*(pt->time-pt->tprop-fT0);
    if (rd < 0) rd = 0;
    if (rd > fgRStraw) rd = fgRStraw;
    double dist = (pt->y-fDyDx*pt->x-fY0)/sqrt(1+fDyDx*fDyDx)-rd*pt->drs;
    chi2        += dist*dist/res2;
    n_good_hits += 1;
  }

  fChi2Dof = chi2/(n_good_hits - 2.999999);
  return 0;
}

//-----------------------------------------------------------------------------
void TrkSegment::print(std::ostream& Stream) {

  int nhits = points.size();

  Stream << std::format("DyDx:{:10.5f} Y0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f}",fDyDx,fY0,fNx,fNy)
         << std::format("-- TrkSegment::{}: nhits:{} fIhit[0]:{} fIhit[1]:{} n_transitions:{} chi2/DOF:{:8.3f}",
                        __func__,nhits,fIhit[0],fIhit[1],fNTransitions,fChi2Dof)
         << std::endl;

  for (int i=0; i<nhits; i++) {
    Point2D* pt = &points[i];
    pt->print();
  }
}
