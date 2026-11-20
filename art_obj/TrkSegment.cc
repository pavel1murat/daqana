//
#include <format>
#include <ostream>

#include "TROOT.h"

#include "cetlib_except/exception.h"
#include "daqana/art_obj/TrkSegment.hh"

int TrkSegment::fgDebugMode(0);

double const TrkSegment::fgRStraw  =  2.5;      // mm

//-----------------------------------------------------------------------------
TrkSegment::TrkSegment(int Plane, int Panel) {
  Clear();
  fPlane      = Plane;
  fPanel      = Panel;
  fTrkPanel    = nullptr;
}

//-----------------------------------------------------------------------------
// the hits keep the index of their ComboHit; the pointers are set here, once the
// collection is known (after reading the segments from a file, or after making them)
//-----------------------------------------------------------------------------
void TrkSegment::setComboHits(mu2e::ComboHitCollection const* Hits) {
  for (auto& sgh : fListOfHits) {
    if ((Hits == nullptr) or (sgh.fIndex < 0)) {
      sgh.fComboHit = nullptr;
      continue;
    }
    sgh.fComboHit = &Hits->at(sgh.fIndex);
    sgh.fStraw    = sgh.fComboHit->strawId().straw();
  }
}

//-----------------------------------------------------------------------------
int TrkSegment::NHitsWithComboHit() const {
  int n(0);
  for (auto const& sgh : fListOfHits) {
    if (sgh.fIndex >= 0) n += 1;
  }
  return n;
}

//-----------------------------------------------------------------------------
SegmentHit* TrkSegment::NextHit(int I) const {
  SegmentHit* next{nullptr};

  int ns = NStraws();
  for (int i=I+1; i<ns; i++) {
    SegmentHit* hit = Hit(i);
    if (hit->ComboHit() != nullptr) {
      next = hit;
      break;
    }
  }
  return next;
}

//-----------------------------------------------------------------------------
// the hits are expected to be sorted by straw, and to have their ComboHits set
// TrkSegments have a well defined pattern with just one transition between the
// layers
//-----------------------------------------------------------------------------
void TrkSegment::fillMissingStraws() {

  int nstraws = NStraws();

  std::vector<SegmentHit> filled;
  filled.push_back(*Hit(0));

  for (int i=1; i<nstraws; i++) {
    SegmentHit* sh_im1 = Hit(i-1);
    SegmentHit*  sh_i   = Hit(i);
    int is_im1 = Hit(i-1)->Straw();
    int is_i   = Hit(i)->Straw();
    if (sh_i->Layer() == sh_im1->Layer()) {
      // hits in the same layer, add potentially missing straws in that layer
      for (int is=is_im1+2; is<is_i; is+=2) {
        SegmentHit gap;
        gap.fIndex = -1;
        gap.fStraw = is;
        gap.fMask  = kNoComboHitBit;
        filled.push_back(gap);
      }
    }
    else {
      // dont' know where the transition had occured, add all straws in the gap
      for (int is=is_im1+1; is<is_i; is++) {
        mu2e::StrawId sid(fPlane,fPanel,is);
        SegmentHit gap(-1,sid);
        gap.fMask  = kNoComboHitBit;
        filled.push_back(gap);
      }
    }
    filled.push_back(*sh_i);
  }

  fListOfHits = filled;
}

//-----------------------------------------------------------------------------
void TrkSegment::UpdateHitParameters() {
  for (auto& sgh : fListOfHits) {
    sgh.fR    = R        (&sgh);
    sgh.t     = DriftTime(&sgh);
    sgh.fDoca = Doca     (&sgh);
  }
}

//-----------------------------------------------------------------------------
void TrkSegment::Clear(){
  fIndex     = -1;
  fMask      =  0;
  fIhit[0]   = -1;
  fIhit[1]   = -1;
  fListOfHits.clear();
  fNghl[0]   =  0;
  fNghl[1]   =  0;
  fNmhl[0]   =  0;
  fNmhl[1]   =  0;
}

//-----------------------------------------------------------------------------
// the hits are expected to be already ordered in the ascending straw order
// when InitHits is called, empty straws are already added
//-----------------------------------------------------------------------------
int TrkSegment::InitHits(int DoTransform) {
  int rc(0);

  int nstraws = NStraws();
  //  int n_hits  = NHitsWithComboHit();

  std::vector<int> drs (nstraws);         // drift signs
//-----------------------------------------------------------------------------
// find seed hits, ignore the very first and the very last pairs
// in the beginning, set all hit flag bits to zero
/// mu2e::StrawHitFlagDetail::dead : back-and-forth transition 
//-----------------------------------------------------------------------------
  fNTransitions = 0;
  const mu2e::ComboHit* last_hit = nullptr;
  int                   last_i(-1);

  for (int i=0; i<nstraws; i++) {
    SegmentHit* sgh = Hit(i);
    if (not sgh->IsGood())                                  continue;

    // figure out the drift signs
    const mu2e::ComboHit* hit = sgh->ComboHit();
    int layer = hit->strawId().getLayer();
    if (last_hit != nullptr) {
      int last_hit_layer = last_hit->strawId().getLayer();

      // layer is either 0 or 1
      if (last_hit_layer != layer) {
//-----------------------------------------------------------------------------
// transition between the layers found.
// ignore back-and-forth transitions - instead, mark hit as bad
//-----------------------------------------------------------------------------
        if (i < nstraws-1) {
          SegmentHit* next = NextHit(i);
          // figure next hit
          
          const mu2e::ComboHit* next_hit = next->ComboHit();
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
        if ((i == 1) or (i == nstraws-1)) {
          last_hit = hit;
          continue;
        }
//-----------------------------------------------------------------------------
// get here only if hits is initially good (not 'back-n-forth')
// rotations don't change the Z-coordinate
//-----------------------------------------------------------------------------
        fNTransitions += 1;
                                        // z becomes Y in points

        double dz = last_hit->pos().z()-hit->pos().z();
        if (DoTransform) dz = dz*fTrkPanel->dsToPanel().rotation().zz();

        if (dz > 0) {
          drs[last_i]      = -1; // those are uint16_t, so use sign+2
          drs[i]           =  1;
        }
        else {
          drs[last_i]      =  1;
          drs[i]           = -1;
        }
        // 
        fIhit[0] = last_i; // this doesnt look right: i-1;
        fIhit[1] = i;
      }
    }
    last_hit    = hit;
    last_i      = i;
  }

  int last_straw[2] = {-1, -1};
  if (fNTransitions != 1) fMask |= kNTransitionsBit;
  for (int i=0; i<nstraws; i++) {
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
  for (int i=0; i<nstraws; i++) {
    SegmentHit*           sgh = Hit(i);
    const mu2e::ComboHit* ch  = sgh->ComboHit();
    if (ch == nullptr) {
//-----------------------------------------------------------------------------
// a straw without a ComboHit: its wire position, the hit is not used in the fit
//-----------------------------------------------------------------------------
      if (fTrkPanel == nullptr) {
        throw cet::exception("TrkSegment") << "TrkSegment::InitHits: the panel is needed for the straws without a ComboHit\n";
      }
      CLHEP::Hep3Vector wire = fTrkPanel->getStraw(sgh->fStraw).getMidPoint();
      if (DoTransform) wire  = fTrkPanel->dsToPanel()*wire;
      sgh->x          = wire.y();
      sgh->y          = wire.z();
      sgh->drs        = 0;
      sgh->sign_fixed = 0;
      sgh->fChi2      = -1.;
      continue;
    }
//-----------------------------------------------------------------------------
// points have coordinates in the local coordinate system of the panel
// add all hits, including flagged ones - those will not be used in the fit
//-----------------------------------------------------------------------------
    CLHEP::Hep3Vector     posl = ch->posCLHEP();
    if (DoTransform) {
      posl            = fTrkPanel->dsToPanel()*posl;
    }
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
//-----------------------------------------------------------------------------
// time() is the early end time, assum ptime is the propagation time for that end
//-----------------------------------------------------------------------------
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
    std::cout << std::format("  i ChIndex      Xloc        Yloc         T\n");
    std::cout << std::format("-------------------------------------------\n");
  }

  for (int i=0; i<nstraws; i++) {
    SegmentHit* sgh  = Hit(i);
    sgh->x -= fXMean;
    sgh->y -= fYMean;
    const mu2e::ComboHit* ch = sgh->ComboHit();
    sgh->t  = (ch == nullptr) ? 0. : ch->time()-ch->propTime()-fTMean;

    if (fgDebugMode != 0) {
      std::cout << std::format(" {:2d} {:5d} {:10.3f} {:10.3f} {:10.3f}\n",
                               i,sgh->fIndex,sgh->x,sgh->y,sgh->t);
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

  int    n_good_hits{0};              // total number of good hits
  int    nghl[2] {};                  // number of good hits per layer
  double chi2{0};

  double t0    = T0();
  double a     = DyDx();
  double b     = Y0();

  int    nstraws = NStraws();
  for (int i=0; i<nstraws; ++i) {
    SegmentHit* sgh = Hit(i);
    sgh->fChi2 = -1;
    if (sgh->IsGood() == 0)                                 continue;
    const mu2e::ComboHit* ch = sgh->ComboHit();
    double rd   = R(sgh,t0);
    if (rd < 0) rd = 0;
    if (rd > fgRStraw) rd = fgRStraw;
    double dist  = (a*sgh->x + b - sgh->y)/sqrt(1+a*a)-rd*sgh->drs;
    double sigr2 = sgh->fSigmaR*sgh->fSigmaR;
    sgh->fChi2   = dist*dist/sigr2; 
    chi2        += sgh->fChi2;
    n_good_hits += 1;
    nghl[ch->strawId().layer()] += 1;
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

  int nstraws = NStraws();
  for (int i=0; i<nstraws; i++) {
    SegmentHit* sgh = Hit(i);
    if (sgh->IsGood() != 1)              continue;
    t0          += sgh->ComboHit()->correctedTime();
    n_good_hits += 1;
  }
                                         // in principle, this approximation for T0 should be good enough
  t0 = t0/(n_good_hits+1.e-12) - fTMean;
//-----------------------------------------------------------------------------
// make sure t0 is not larger than the smallest measured time
//-----------------------------------------------------------------------------
  for (int i=0; i<nstraws; i++) {
    SegmentHit* sgh = Hit(i);
    if (sgh->IsGood() != 1)              continue;
    if (t0 > sgh->t) t0 = sgh->t;
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
         << std::format("Segment {:3d} parameters:",fIndex);
  if (Message) Stream << " " << Message;
  Stream << std::endl;
                                        // no points, if N(hits) < 4...
  int nstraws = NStraws();
  int n_hits  = NHitsWithComboHit();

  Stream << std::format(" panel:{:02}:{} nstraws:{:2} nhits:{} ngood:{} n_transitions:{} seed hits: {}:{} ngood/layer: {:2}:{:2} misses/layer: {:2}{:2}\n",
                        fPlane,fPanel,nstraws,n_hits,fNGoodHits,fNTransitions,fIhit[0],fIhit[1],fNghl[0],fNghl[1],fNmhl[0],fNmhl[1])
         << std::format("DyDx:{:10.5f} Y0:{:10.5f} Nx:{:10.5f} Ny:{:10.5f} T0:{:10.3f} chi2/DOF:{:8.3f}\n",
                        fPar.a,fPar.b,fPar.nx,fPar.ny,fPar.tau+fTMean,fPar.chi2dof);

  Stream << std::format(" i   sid  pnl:str  mask   ch_mask       x      y        t      time      tprop  tdrift  radius drs      dR       doca        sigr    chi2   fixed\n");
  Stream <<             "------------------------------------------------------------------------------------------------------------------------------------\n" ;
  double t0 = fPar.T0();
  for (int i=0; i<nstraws; i++) {
    SegmentHit* sgh = Hit(i);
    const mu2e::ComboHit* ch = sgh->ComboHit();
    if (ch == nullptr) {
      Stream << std::format("{:2d}        {:2}:{}:{:2} 0x{:04x}            {:8.3f} {:6.3f}                                                        {:12.4f}\n",
                            i,fPlane,fPanel,sgh->fStraw,sgh->fMask,sgh->x,sgh->y,sgh->fDoca);
      continue;
    }
    mu2e::StrawId sid = ch->strawId();
    int ch_flag = std::stoi(ch->flag().hex());
    Stream << std::format("{:2d} 0x{:04x} {:2}:{}:{:2} 0x{:04x} 0x{:08x}",
                          i,sid.asUint16(),sid.plane(),sid.panel(),sid.straw(),sgh->fMask, ch_flag);
    Stream << std::format(" {:8.3f} {:6.3f} {:7.3f} {:10.3f} {:7.3f} {:7.3f} {:7.3f} {:2d}",
                          sgh->x,sgh->y,sgh->t,ch->time(),ch->propTime(),
                          DriftTime(sgh,t0),R(sgh,t0),sgh->drs);
    Stream << std::format(" {:10.4f} {:10.4f} {:7.3f} {:8.2f} {:5d}",
                          Dr(sgh,&fPar),Doca(sgh,&fPar),sgh->fSigmaR,sgh->fChi2,sgh->sign_fixed)
           << std::endl;
  }
}
