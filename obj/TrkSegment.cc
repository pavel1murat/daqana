//
#include <format>
#include <ostream>
#include "daqana/obj/TrkSegment.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

int TrkSegment::fgDebugMode(0);

//
double       SegmentHit::fgTOffset =  0.0;      // 1.6;       // ns
double       SegmentHit::fgVDrift  = 61.3e-3;   // mm/ns
double const TrkSegment::fgRStraw  =  2.5;      // mm
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
  fListOfHits.clear();
  // points.clear();
  fNghl[0]    = 0;
  fNghl[1]    = 0;
  fNmhl [0]   = 0;
  fNmhl [1]   = 0;
}

//-----------------------------------------------------------------------------
// expect 'ListOfHits' to be already ordered in the ascending straw order
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
      SegmentHit sgs(ch);
      fListOfHits.emplace_back(sgs);
    }
  }
  else {
                                        // MakeDigiNtuple
    nhits = nHits();
  }

  std::vector<int> drs (nhits);
//-----------------------------------------------------------------------------
// find seed hits, ignore the very first and the very last pairs
// in the befinning, set all hit flag bits to zero
/// mu2e::StrawHitFlagDetail::dead : back-and-forth transition 
//-----------------------------------------------------------------------------
  fNTransitions = 0;
  const mu2e::ComboHit* last_hit = nullptr;
  int                   last_i(-1);

  for (int i=0; i<nhits; i++) {
    SegmentHit* sgh = Hit(i);
    if (not sgh->IsGood())                                  continue;
    const mu2e::ComboHit* hit = sgh->ComboHit();
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
          const mu2e::ComboHit* next_hit = Hit(i+1)->ComboHit();
          int next_hit_layer = next_hit->strawId().getLayer();
          if (next_hit_layer == last_hit_layer) {
            sgh->fMask |= TrkSegment::kNTransitionsBit;
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
    SegmentHit* sgh = Hit(i);
    if (not sgh->IsGood())                                  continue;

    const mu2e::ComboHit* hit = sgh->ComboHit();
    int layer = hit->strawId().getLayer();
    int straw = hit->strawId().straw();
    fNghl[layer] += 1;
    if (last_straw[layer] != -1) {
      int miss = (straw-last_straw[layer]) > 2;
      fNmhl[layer]  += miss;
    }
    last_straw[layer] = straw;
  }

  if (fgDebugMode != 0) {
    std::cout << std::format("fIhit[0]:{:2d} fIhit[1]:{:2d} ngh: {} {} miss: {} {}\n",
                             fIhit[0],fIhit[1],fNghl[0],fNghl[1],fNmhl[0],fNmhl[1]);
  }

  fXMean = 0;
  fYMean = 0;
  fTMean = 0;

  if (fgDebugMode != 0) {
    std::cout << std::format(" i   sid   mask pnl:str ch_flag               XYZ(M)                              XYZ(L)                 time    prtime  rdrift drs\n");
    std::cout << std::format("--------------------------------------------------------------------------------------\n");
  }

  fNGoodHits = 0;
  for (int i=0; i<nhits; i++) {
    SegmentHit*           sgh  = Hit(i);
    const mu2e::ComboHit* ch   = sgh->ComboHit();
    CLHEP::Hep3Vector     posl = fTrkPanel->dsToPanel()*ch->posCLHEP();
//-----------------------------------------------------------------------------
// points have coordinates in the local coordinate system of the panel
// add all hits, including flagged ones - those will not be used in the fit
//-----------------------------------------------------------------------------
    sgh->x          = posl.y();
    sgh->y          = posl.z();
    sgh->drs        = drs[i];
    sgh->sign_fixed = 0;
    sgh->fChi2      = -1.;
//-----------------------------------------------------------------------------
// to avoid mistakes, perform the subtraction just once
// hits flagged at this stage will never be unflagged
// and don't assume that <X>, <Y> etc are exactly equal to zero
//-----------------------------------------------------------------------------
    if (sgh->IsGood()) {
      fNGoodHits += 1;
      fXMean     += posl.y();
      fYMean     += posl.z();
                                        // time() is the early end time, assum ptime is the propagation time for that end
      fTMean     += ch->time()-ch->propTime();
    }
    if (fgDebugMode != 0) {
      int ch_flag = std::stoi(ch->flag().hex(),nullptr,0);
      mu2e::StrawId sid = ch->strawId();
      std::cout << std::format("{:2d} 0x{:04x} 0x{:04x} 0x{:08x} {:2}:{}:{:2}  ({:10.3f} {:10.3f} {:10.3f})",
                               i,sid.asUint16(),sgh->fMask,ch_flag, sid.plane(),sid.panel(),sid.straw(),
                               ch->pos().x(), ch->pos().y(), ch->pos().z())
                << std::format(" ({:10.3f} {:10.3f} {:10.3f}) {:10.3f} {:7.3f} {:2}",
                               posl[0], posl[1], posl[2], ch->correctedTime(), ch->propTime(), drs[i])
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
    std::cout << std::format("  i       Xloc       Yloc        T\n");
    std::cout << std::format("------------------------------------\n");
  }

  for (int i=0; i<nhits; i++) {
    // Point2D* p  = &points[i];
    SegmentHit* p  = &fListOfHits[i];
    p->x          -= fXMean;
    p->y          -= fYMean;
    p->t           = p->ComboHit()->time()-p->ComboHit()->propTime()-fTMean;

    if (fgDebugMode != 0) {
      std::cout << std::format(" {:2d} {:10.3f} {:10.3f} {:10.3f}\n",i,p->x,p->y,p->t);
    }
  }
//-----------------------------------------------------------------------------
// finally, define tangent line. The hits are already sorted
//-----------------------------------------------------------------------------
  rc = DefineTangentLine();

  if (fgDebugMode != 0) {
    print("InitHits after DefineTangentLine");
    std::cout << std::format("-- TrkSegment.cc:{} END rc:{}\n",__LINE__,rc);
  }
  return rc;
}

//------------------------------------------------------------------------------
// in general, the worst hit could be not the one which has the largest
// updates hit contribuitons to the total chi2 and the numbers of good hits
//-----------------------------------------------------------------------------
double TrkSegment::Chi2() {
  double const res2(0.2*0.2);            // assume 200 um

  double chi2 = 0;
  int nhits = nHits();
  int n_good_hits(0);
  int nghl[2] = {};
  double t0 = T0();
  double a  = DyDx();
  double b  = Y0();

  for (int i=0; i<nhits; ++i) {
    SegmentHit* sgh = Hit(i);
    sgh->fChi2 = -1;
    if (sgh->IsGood() == 0)                                 continue;
    double rd   = R(sgh,t0);
    if (rd < 0) rd = 0;
    if (rd > fgRStraw) rd = fgRStraw;
    double dist = (a*sgh->x + b - sgh->y)/sqrt(1+a*a)-rd*sgh->drs;
    sgh->fChi2   = dist*dist/res2; 
    chi2        += sgh->fChi2;
    n_good_hits += 1;
    nghl[sgh->fComboHit->strawId().layer()] += 1;
  }

  fNGoodHits = n_good_hits;
  fNghl[0]   = nghl[0];
  fNghl[1]   = nghl[1];
  chi2       = chi2/(n_good_hits - 2.999999);
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

  if ((i1 < 0) or (i2 < 0)) {
                                        // don't do anything for 1-layer segments
    fTangentLine.a       = -999.;
    fTangentLine.b       = -999.;
    fTangentLine.nx      = -999;
    fTangentLine.ny      = -999;
    fTangentLine.tau     = -999;
    fTangentLine.chi2dof = -1;
    fPar = fTangentLine;
    return -1;
  }

  SegmentHit* p1 = Hit(i1);
  SegmentHit* p2 = Hit(i2);
//-----------------------------------------------------------------------------
// define T0
//-----------------------------------------------------------------------------
  double t0(0);
  int    n_good_hits(0);

  int nhits = nHits();
  for (int i=0; i<nhits; i++) {
    SegmentHit* pt = Hit(i);
    if (pt->IsGood() != 1)              continue;
    t0          += pt->fComboHit->correctedTime();
    n_good_hits += 1;
  }
                                         // in principle, this approximation for T0 should be good enough
  t0 = t0/(n_good_hits+1.e-12) - fTMean;
//-----------------------------------------------------------------------------
// make sure t0 is not larger than the smallest measured time
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; i++) {
    SegmentHit* pt = Hit(i);
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

  Stream << "-----------------------------------------------------------------------------------------------";
  Stream << std::endl
         << "Segment parameters:";
  if (Message) Stream << " " << Message;
  Stream << std::endl;
                                        // no points, if N(hits) < 4...
  int nhits = nHits();

  Stream << std::format(" panel:{:02}:{} nhits:{} ngood:{} n_transitions:{} seed hits: {}:{} ngood/layer: {:2}:{:2} misses/layer: {:2}{:2}\n",
                        fPlane,fPanel,nhits,fNGoodHits,fNTransitions,fIhit[0],fIhit[1],fNghl[0],fNghl[1],fNmhl[0],fNmhl[1])
         << std::format("DyDx:{:10.5f} Y0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f} T0:{:10.3f} chi2/DOF:{:8.3f}\n",
                        fPar.a,fPar.b,fPar.nx,fPar.ny,fPar.tau+fTMean,fPar.chi2dof);

  Stream << std::format(" i   sid  pnl:str  mask   ch_mask       x      y        t      time      tprop  tdrift  radius drs      rho     chi2   fixed\n");
  Stream <<             "----------------------------------------------------------------------------------------------------------------------------\n" ;
  double t0 = fPar.T0();
  for (int i=0; i<nhits; i++) {
    SegmentHit* sgh = Hit(i);
    mu2e::StrawId sid = sgh->ComboHit()->strawId();
    int ch_flag = std::stoi(sgh->ComboHit()->flag().hex());
    Stream << std::format("{:2d} 0x{:04x} {:2}:{}:{:2} 0x{:04x} 0x{:08x} {:8.3f} {:6.3f} {:7.3f} {:10.3f} {:7.3f} {:7.3f} {:7.3f} {:2d} {:10.4f} {:8.2f} {:5d}",
                          i,sid.asUint16(),sid.plane(),sid.panel(),sid.straw(),sgh->fMask, ch_flag,
                          sgh->x,sgh->y,sgh->t,sgh->fComboHit->time(),sgh->fComboHit->propTime(),
                          DriftTime(sgh,t0),R(sgh,t0),sgh->drs,Rho(sgh,&fPar),sgh->fChi2,sgh->sign_fixed)
           << std::endl;
  }
}
