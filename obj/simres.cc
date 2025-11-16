///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "daqana/obj/simres.hh"


//-----------------------------------------------------------------------------
simres::simres(int DebugMode) {
  fDebugMode   = DebugMode;
  fStrawRadius = 2.5;                   // mm
  fVDrift      = fStrawRadius/40.;      // max 'legal' drift time of 40 ns
}

//-----------------------------------------------------------------------------
// R in mm, t in ns
//-----------------------------------------------------------------------------
double simres::D2T(double R) {
  double t = R/fVDrift;
  return t;
}

//-----------------------------------------------------------------------------
// R in mm, t in ns
//-----------------------------------------------------------------------------
double simres::T2D(double T) {
  double r = T*fVDrift;
  return r;
}

//-----------------------------------------------------------------------------
int simres::init_geometry() {

  double pitch = 6.25;              // step between straws in the layer
  for (int is=0; is<96; ++is) {
    Straw_t* s = new Straw_t(is);
    s->fY      = 3.*(2*s->fLayer-1);
    s->fX      = (s->fNumber-24)*pitch + pitch/2*s->fLayer; // layers are staggered

    fPanel.fListOfStraws.emplace_back(s);

    std::cout << std::format("is, layer, y, x: {:2d} {} {:10.3f} {:10.3f}\n",
                             is, s->fLayer,s->fX,s->fY);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int simres::book_histograms() {
  fHist.fNHits = new TH1F("nhits","nhits",20,0,20);
  return 0;
}

//-----------------------------------------------------------------------------
int simres::fill_histograms(Hist_t* Hist) {
  Hist->fNHits->Fill(fPanel.fNHits);
  return 0;
}

//-----------------------------------------------------------------------------
int simres::run(int NEvents, int DebugMode) {
  init_geometry();

  book_histograms();

  // simulation per se
  TRandom3 trn3;

  // double y0, dydx;

  for (int ievent=0; ievent<NEvents; ++ievent) {
    fPanel.fNHits = 0;
    fPanel.fListOfHitStraws.clear();
//----------------------------------------------------------------------------
// 1) simulate the line parameters
//-----------------------------------------------------------------------------
    fT0   = 0.; // trn3.Gaus(0,10); // turn off the T0 smearing to begin with
    fX0   = 0;                                                 // 
    fY0   = 0; // (2*trn3.Rndm()-1)*1.0; // Y at X=0, in mm
    fDyDx = 0.4; // trn3.Gaus(0.4,0.1);

    std::cout << std::format("[{}:{}] y0:{:10.4f} dydx:{:10.4f}\n",__func__,__LINE__,fY0,fDyDx);

    fNx   = 1./sqrt(1+fDyDx*fDyDx);
    fNy   = fDyDx*fNx;
    fNux  = -fNy;
    fNuy  = fNx;

    std::cout << __func__ << ":" << __LINE__ << ": nx :" << fNx  << " ny :" << fNy  << std::endl;
    std::cout << __func__ << ":" << __LINE__ << ": nux:" << fNux << " nuy:" << fNuy << std::endl;
//-----------------------------------------------------------------------------
// 2) simulate hits
//-----------------------------------------------------------------------------
    for (int is=0; is<96; is++) {
      std::cout << __func__ << ":" << __LINE__ << ": ------------------------- is:" << is << std::endl;
      Straw_t* s = Straw(is);

      double dr = (fY0-s->fY)*fNuy+(fX0-s->fX)*fNux; // distance (signed) from the track to the wire

      std::cout << " s->fX:" << s->fX << " s->fY:" << s->fY << " dr:" << dr << std::endl;

      if (fabs(dr) < fStrawRadius) {
        s->fDoca = fabs(dr);
        fPanel.fNHits++;
        fPanel.fListOfHitStraws.push_back(s);
        std::cout << std::format("-- hit in straw:{:2d} x0:{:10.3f}  y0:{:10.3f} doca:{:10.3f}\n",
                                 is,s->fX,s->fY,s->fDoca);
      }
      std::cout <<  __func__ << ":" << __LINE__ << ": fNHits:" << fPanel.fNHits << std::endl;
    }
//-----------------------------------------------------------------------------
// 3) reconstruct the line
//    construct the segment - need to fake a list of ComboHits
//-----------------------------------------------------------------------------
    if (fPanel.fNHits > 2) {
      for (int i=0; i<fPanel.fNHits; ++i) {
        Straw_t* s          = HitStraw(i);
        double   drift_time = D2T(s->fDoca);
        
        mu2e::ComboHit ch;
        ch._eend                       = mu2e::StrawEnd::cal;
        ch._etime[mu2e::StrawEnd::cal] = fT0+drift_time;
        ch._ptime                      = 0;

        fPanel.fListOfComboHits.emplace_back(ch);

        SegmentHit sgh(&ch);
        fSegment.fListOfHits.emplace_back(sgh);
      }
    }
//-----------------------------------------------------------------------------
// 4) at this point, should be ready to fit the segment
//-----------------------------------------------------------------------------
    TrkSegment::Par_t par;

    SegmentFit sfitter(&fSegment);
    sfitter.DefineDriftDirections();
    int niter(6);
    int converged = sfitter.Fit(niter,0,nullptr,&par);
    std::cout << std::format("converged:{}\n",converged);
  }

  fill_histograms(&fHist);

  return 0;
}

//-----------------------------------------------------------------------------
int simres::display_event() {
//-----------------------------------------------------------------------------
// plot the results after the last iteration, but before updating the points
//-----------------------------------------------------------------------------
  if (fDebugMode & 0x1) std::cout << std::format("-- {}:{} START\n",__func__,__LINE__);

  int    nhits = fPanel.fListOfHitStraws.size();

  TCanvas* c = new TCanvas("c","c",1500,600);
  TH2F* h2 = new TH2F("h2","",500,-25,25,200,-10,10);
  h2->SetStats(kFALSE);
  h2->Draw();
  //  c->Range(-150,-600,150,60);
  c->cd();
//-----------------------------------------------------------------------------
// draw the line
//-----------------------------------------------------------------------------
  double xmin(-150.), xmax(150.);
  double ymin = fY0+(xmin-fX0)*fDyDx;
  double ymax = fY0+(xmax-fX0)*fDyDx;
  TLine* l = new TLine(xmin,ymin,xmax,ymax);
  l->SetLineColor(kRed+2);
  l->Draw();

  for (int i=0; i<96; i++) {
                                        // draw all straws, all - in mm
    Straw_t* s = Straw(i);
    TEllipse* e = new TEllipse(s->fX,s->fY,fStrawRadius,0,0,360);
    e->SetLineColor(kBlack);
    e->SetFillStyle(0);
    e->Draw();
  }
                                        // draw hits
  for (int i=0; i<nhits; i++) {
    Straw_t* s = HitStraw(i);
                                        // hit straw

    TEllipse* e = new TEllipse(s->fX,s->fY,fStrawRadius,0,0,360);
    e->SetLineColor(kRed);
    e->SetFillStyle(0);
    e->Draw();
                                        // hit
    double r = s->fDoca;

    TEllipse* e2 = new TEllipse(s->fX,s->fY,r,0,0,360);
    if (s->fMask == 0) {
      e2->SetFillColor(kBlue-10);
      e2->SetFillStyle(3003);
    }
    else {
      e2->SetFillColor(kRed-10);
      e2->SetFillStyle(3001);
    }
    e2->Draw();
  }

  gPad->Modified();
  gPad->Update();

  if (fDebugMode & 0x1) std::cout << std::format("-- simres::{}:{} END\n",__func__,__LINE__);
  return 0;

}
