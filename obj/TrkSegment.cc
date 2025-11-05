//
#include <format>
#include <ostream>
#include "daqana/obj/TrkSegment.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

int TrkSegment::fgDebugMode(0);

//
double       Point2D::fgTOffset   =  0.0;      // 1.6;       // ns
double       Point2D::fgVDrift    = 61.3e-3;   // mm/ns
double const TrkSegment::fgRStraw =  2.5;      // mm
//-----------------------------------------------------------------------------
TrkSegment::TrkSegment(int Plane, int Panel) {
  Clear();
  fPlane      = Plane;
  fPanel      = Panel;
  //  fCombiTrans = nullptr; // geometry - to be initialized just once
  fTrkPanel    = nullptr;
}

//-----------------------------------------------------------------------------
void TrkSegment::Clear(){
  fMask       = 0;
  fIhit[0]    = -1;
  fIhit[1]    = -1;
  hits.clear();
  points.clear();
  fNghl[0] = 0;
  fNghl[1] = 0;
  fNmhl [0] = 0;
  fNmhl [1] = 0;
}

//-----------------------------------------------------------------------------
int TrkSegment::InitHits(std::vector<const mu2e::ComboHit*>* ListOfHits, int UniquePlane, int Panel) {
  int rc(0);

  
  int nhits(-1);
  if (ListOfHits) {
                                        // standalone test
    nhits = ListOfHits->size();
    std::cout << std::format("TrkSegment::{}:START UniquePlane:{} Panel:{} plane:{} panel:{} nhits:{}\n",
                           __func__,UniquePlane,Panel,fPlane,fPanel,nhits);

                                        // copy pointers to hits and dont worry about the memory management
                                        // hits are owned by the caller
    for (int i=0; i<nhits; i++) {
      const mu2e::ComboHit* ch = ListOfHits->at(i);
      hits.emplace_back(ch);
    }
  }
  else {
                                        // MakeDigiNtuple
    nhits = hits.size();
  }
  std::vector<int> drs (nhits);
  std::vector<int> flag(nhits);
//-----------------------------------------------------------------------------
// sort combohits
//-----------------------------------------------------------------------------
  std::sort(hits.begin(),hits.end(),
            [] (const mu2e::ComboHit* a, const mu2e::ComboHit* b) {
              return a->strawId().straw() < b->strawId().straw();
            });

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
// in the befinning, set all hit flag bits to zero
/// mu2e::StrawHitFlagDetail::dead : back-and-forth transition 
//-----------------------------------------------------------------------------
  fNTransitions = 0;
  const mu2e::ComboHit* last_hit = nullptr;
  int                   last_i(-1);

  for (int i=0; i<nhits; i++) {
    const mu2e::ComboHit* hit = hits[i];
                                        // figure out the drift signs
    int layer = hit->strawId().getLayer();
    if (last_hit != nullptr) {
      int last_hit_layer = last_hit->strawId().getLayer();
                                        // layer is either 0 or 1
      if (last_hit_layer != layer) {
//-----------------------------------------------------------------------------
// transition between the layers found.
// ignore back-and-forth transitions - instead, mark hit as bad
//-----------------------------------------------------------------------------
        if (i < nhits-1) {
          const mu2e::ComboHit* next_hit = hits[i+1];
          int next_hit_layer = next_hit->strawId().getLayer();
          if (next_hit_layer == last_hit_layer) {
            flag[i] |= 0x1;
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
        double dz = (last_hit->pos().z()-hit->pos().z())*fTrkPanel->dsToPanel().rotation().zz();
        //        if (last_hit->pos().z() > hit->pos().z()) {
        if (dz > 0) {
          drs[last_i]      = -1; // those are uint16_t, so use sign+2
          drs[i]           =  1;
        }
        else {
          drs[last_i]      =  1;
          drs[i]           = -1;
        }
        fIhit[0] = i-1;
        fIhit[1] = i;
      }
    }
    last_hit    = hit;
    last_i      = i;
  }

  int last_straw[2] = {-1, -1};
  if (fNTransitions != 1) fMask |= kNTransitionsBit;
  for (int i=0; i<nhits; i++) {
    const mu2e::ComboHit* hit = hits[i];
    if (flag[i] != 0) continue;
    int layer = hit->strawId().getLayer();
    int straw = hit->strawId().straw();
    fNghl[layer] += 1;
    if (last_straw[layer] != -1) {
      int miss = (straw-last_straw[layer]) > 2;
      fNmhl[layer]  += miss;
    }
    last_straw[layer] = straw;
  }

  // for (int i=fIhit[1]; i<nhits; i++) {
  //   const mu2e::ComboHit* hit = hits[i];
  //   if (flag[i] != 0) continue;
  //   int layer = hit->strawId().getLayer();
  //   int straw = hit->strawId().straw();
  //   fNghl[layer] += 1;
  //   int miss = (straw-last_straw[layer]) > 2;
  //   fNmhl[layer]  += miss;
  // }

  if (fgDebugMode != 0) {
    std::cout << std::format("fIhit[0]:{:2d} fIhit[1]:{:2d} ngh: {} {} miss: {} {}\n",
                             fIhit[0],fIhit[1],fNghl[0],fNghl[1],fNmhl[0],fNmhl[1]);
    //    fCombiTrans->Print();
  }
//-----------------------------------------------------------------------------
// add Point2D's :
//-----------------------------------------------------------------------------
  //  int    imax      (-1);
  //  double max_drtime(-1);

  fXMean = 0;
  fYMean = 0;
  fTMean = 0;

  if (fgDebugMode != 0) {
    std::cout << std::format(" i   sid   flag                 XYZ(M)                              XYZ(L)                 time    prtime  rdrift drs\n");
    std::cout << std::format("--------------------------------------------------------------------------------------\n");
  }

  fNGoodHits = 0;
  for (int i=0; i<nhits; i++) {
    const mu2e::ComboHit* ch = hits[i];

    CLHEP::Hep3Vector posl = fTrkPanel->dsToPanel()*ch->posCLHEP();
//-----------------------------------------------------------------------------
// points have coordinates in the local coordinate system of the panel
// add all hits, including flagged ones - those will not be used in the fit
//-----------------------------------------------------------------------------
    AddPoint(ch->_sid.asUint16(),flag[i], posl[1],posl[2],ch->_etime[ch->_eend],ch->_ptime, ch->correctedTime(), drs[i]);
//-----------------------------------------------------------------------------
// to avoid mistakes, perform the subtraction just once
// hits flagged at this stage will never be unflagged
// and don't assume that <X>, <Y> etc are exactly equal to zero
//-----------------------------------------------------------------------------
    if (flag[i] == 0) {
      fNGoodHits += 1;
      fXMean     += posl[1];
      fYMean     += posl[2];
                                        // time() is the early end time, assum ptime is the propagation time for that end
      fTMean     += ch->time()-ch->_ptime;
    }
    // double drtime = ch->time-ch->prtime-fT0;
    // if (drtime > max_drtime) {
    //   // imax       = i;
    //   max_drtime = drtime;
    // }
    if (fgDebugMode != 0) {
      int iflag = std::stoi(ch->flag().hex(),nullptr,0);
      std::cout << std::format("{:2} 0x{:04x} 0x{:05x} ({:10.3f} {:10.3f} {:10.3f})",i,ch->strawId().asUint16(),
                               iflag,ch->pos().x(), ch->pos().y(), ch->pos().z())
                << std::format(" ({:10.3f} {:10.3f} {:10.3f}) {:10.3f} {:7.3f} {:2}",
                               posl[0], posl[1], posl[2], ch->correctedTime(), ch->_ptime, drs[i])
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
  if (fgDebugMode != 0) {
    std::cout << " --- after subtracting mean values\n";
    std::cout << std::format(" i      Xloc      Yloc    T\n");
    std::cout << std::format("----------------------------\n");
  }

  for (int i=0; i<nhits; i++) {
    Point2D* p  = &points[i];
    p->x       -= fXMean;
    p->y       -= fYMean;
    p->t        = p->time-p->tprop-fTMean;

    if (fgDebugMode != 0) {
      std::cout << std::format(" {:2d} {:10.3f} {:10.3f} {:10.3f}\n",i,p->x,p->y,p->t);
    }
  }
//-----------------------------------------------------------------------------
// finally, sort points and define tangent line
// the hits are already sorted, the points should be sorted by construction
//-----------------------------------------------------------------------------
  // std::sort(points.begin(),points.end(),
  //           [] (const Point2D& a, const Point2D& b) {
  //             return a.strawId().straw() < b.strawId().straw();
  //           });

  DefineTangentLine();

  if (fgDebugMode != 0) {
    print("InitHits after DefineTangentLine");
    std::cout << std::format("-- TrkSegment.cc:{} END rc:{}\n",__LINE__,rc);
  }
  return rc;
}

//------------------------------------------------------------------------------
// in general, the worst hit could be not the one which has the largest chi2
//-----------------------------------------------------------------------------
double TrkSegment::Chi2() {
  double const res2(0.2*0.2);            // assume 200 um

  double chi2 = 0;
  int nhits = hits.size();
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
  if (fgDebugMode) std::cout << std::format("-- TrkSegment::{}: BEGIN\n",__func__);
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

  int nhits = hits.size();
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
  // double nx    =  alpdr*dxr-dyr*sqrt(1-alpdr*alpdr);
  // double ny    =  alpdr*dyr+dxr*sqrt(1-alpdr*alpdr);

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

  if (fgDebugMode) {
    std::cout << std::format("-- TrkSegment::{}: END rc:{} fDyDx:{:10.5f} fY0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f}\n",
                             __func__,rc,fTangentLine.a,fTangentLine.b,fTangentLine.nx,fTangentLine.ny);
  }
  return rc;
}

//-----------------------------------------------------------------------------
void TrkSegment::print(const char* Message, std::ostream& Stream) {

  Stream << "--------------------------------------------- Segment parameters:";
  if (Message) Stream << " " << Message;
  Stream << std::endl;
                                        // no points, if N(hits) < 4...
  int nhits = hits.size();
  int npt   = points.size();

  Stream << std::format(" plane:{:2} panel:{} nhits:{} fNGoodHits:{} n_transitions:{} seed hits: {}:{} good hits/layer: {:2}:{:2} misses/layer: {:2}{:2} Npoints:{}\n",
                        fPlane,fPanel,nhits,fNGoodHits,fNTransitions,fIhit[0],fIhit[1],fNghl[0],fNghl[1],fNmhl[0],fNmhl[1],npt)
         << std::format(" Parameters: DyDx:{:10.5f} Y0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f} T0:{:10.3f} chi2/DOF:{:8.3f}\n",
                        fPar.a,fPar.b,fPar.nx,fPar.ny,fPar.tau+fTMean,fPar.chi2dof);

  Stream << std::format(" i   sid   mask       x      y        t     time      tprop  tdrift  radius drs      rho  fixed\n");
  Stream <<             "-----------------------------------------------------------------------------------------------\n" ;
  double t0 = fPar.T0();
  for (int i=0; i<npt; i++) {
    Point2D* pt = &points[i];
    Stream << std::format("{:2d} 0x{:04x} 0x{:04x} {:8.3f} {:6.3f} {:7.3f} {:10.3f} {:7.3f} {:7.3f} {:7.3f} {:2d} {:10.4f} {:5d}",
                          i,pt->sid,pt->fMask,pt->x,pt->y,pt->t,pt->time,pt->tprop,
                          DriftTime(pt,t0),R(pt,t0),pt->drs,Rho(pt,&fPar),pt->sign_fixed)
           << std::endl;
  }
}
