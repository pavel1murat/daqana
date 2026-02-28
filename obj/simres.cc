///////////////////////////////////////////////////////////////////////////////
// 0x1 : intermittent
// 0x2 : hits
// 0x4 : geometry initialization
///////////////////////////////////////////////////////////////////////////////

#include "daqana/obj/simres.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "daqana/obj/TrkSegment.hh"
#include <regex>

#include "TROOT.h"
#include "TProfile.h"
#include "TFile.h"

//-----------------------------------------------------------------------------
// the simulated timing resolution of 2 ns translates into 125 um
//-----------------------------------------------------------------------------
simres::simres(int DebugMode) {
  fDebugMode    = DebugMode;
  fStrawRadius  = 2.5;                   // mm
  fVDrift       = fStrawRadius/40.;      // max 'legal' drift time of 40 ns ... 62.5 um/ns
  fSigmaT0      = 2.;                    // [ns]

  fEvent.number = -1;

  fMaxGoodDy    =  1.;
  fBeamMode     =  0;

  for (int i=0; i<96; i++) {
    fPanel.fStrawOffsetX[i] = 0;
    fPanel.fStrawOffsetY[i] = 0;
  }

  SegmentHit::SetVDrift (fVDrift);
  SegmentHit::SetTOffset(0);

  fFolder  = gROOT->GetRootFolder()->AddFolder("simres","simres");
}

//-----------------------------------------------------------------------------
void simres::Panel_t::Clear() {
  fListOfHitStraws.clear();

  int nch = fListOfComboHits.size();
  for (int i=0; i<nch; ++i) {
    delete fListOfComboHits[i];
  }
  fListOfComboHits.clear();
}

//-----------------------------------------------------------------------------
void simres::print_(const std::string& Message, const std::source_location& location) {

  struct xx_t {
    std::vector<std::string> split(const std::string& str, const std::string& delimiter) {
      std::vector<std::string> result;
      std::regex re(delimiter);
      std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
      std::sregex_token_iterator end;
      while (it != end) {
        result.push_back(*it++);
      }
      return result;
    }
  } xx;

  std::string s = std::format("event:{:6d}",fEvent.number);

  std::vector<std::string> ss = xx.split(location.file_name(),"/");

  std::cout << s << " " << ss.back() << ":" << location.line()
    //            << location.function_name()
            << " : " << Message;
}

//-----------------------------------------------------------------------------
// R in mm, t in ns
//-----------------------------------------------------------------------------
double simres::D2T(double R) {
  double t = R/fVDrift;
  return t;
}

//-----------------------------------------------------------------------------
// R in mm, t in ns - need to use _RECONSTRUCTION T2D
//-----------------------------------------------------------------------------
double simres::T2D(double T) {
  double r = T*fVDrift;
  return r;
}

//-----------------------------------------------------------------------------
int simres::init_geometry() {

  double pitch = 6.25;              // step between straws in the layer
  for (int is=0; is<96; ++is) {
    mu2e::StrawId sid(0,0,is);
    Straw_t* s   = new Straw_t(sid);
    int      lay = sid.layer();
    int      num = is/2;
//-----------------------------------------------------------------------------
// TODO the code needs to be modified to consistently use fPanel
//-----------------------------------------------------------------------------
    s->fY      = 3.*(2*lay-1)+fPanel.fStrawOffsetY[is];
    s->fX      = (num-24)*pitch + pitch/2*lay + fPanel.fStrawOffsetX[is]; // layers are staggered

    fPanel.fListOfStraws.emplace_back(s);

    if (fDebugMode & 0x4) {
      print_(std::format("is:{:2d} layer:{} num:{:2d} y:{:10.3f} x:{:10.3f}\n",is, lay,num, s->fX,s->fY));
    }
  }
  return 0;
}




//-----------------------------------------------------------------------------
int simres::book_hit_histograms(HitHist_t* Hist, int ISet, const char* Folder) {
  
  HBook1F(Hist->r         ,Form("r_%02i"         ,ISet),"rdrift, mm", 100,  0, 5, Folder);
  HBook1F(Hist->doca      ,Form("doca_%02i"      ,ISet),"doca, mm"  , 200,  -5, 5, Folder);
  HBook1F(Hist->dr        ,Form("dr_%02i"        ,ISet),"dr, mm"    , 500,  -2.5, 2.5, Folder);
  HBook1F(Hist->drho      ,Form("drho_%02i"      ,ISet),"drho, mm"  , 500,  -2.5, 2.5, Folder);
  HBook2F(Hist->dr_vs_r   ,Form("dr_vs_r_%02i"   ,ISet),"dr_vs_r"   , 100,  0,  5, 250,-2.5,2.5,Folder);
  HBook2F(Hist->drho_vs_r ,Form("drho_vs_r_%02i" ,ISet) ,"drho_vs_r", 100,  0,  5, 250,-2.5,2.5,Folder);
  HBook2F(Hist->dr_vs_is  ,Form("dr_vs_is_%02i"  ,ISet)  ,"dr_vs_is", 100,  0,  100, 250,-2.5,2.5,Folder);
  HBook2F(Hist->drho_vs_is,Form("drho_vs_is_%02i",ISet),"drho_vs_is", 100,  0,  100, 250,-2.5,2.5,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int simres::book_segment_histograms(SegmentHist_t* Hist, int ISet, const char* Folder) {
  
  HBook1F(Hist->dy , Form("dy_%i" ,ISet), "Y0(reco)-Y0(sim)"    ,100, -5, 5,Folder);
  HBook1F(Hist->ddy, Form("ddy_%i",ISet), "DyDx(reco)-DyDx(sim)",200, -0.1, 0.1,Folder);
  HBook1F(Hist->dt0, Form("dt0_%i",ISet), "T0(reco)-T0(sim)"    ,200,-10,10,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int simres::book_event_histograms(EventHist_t* Hist, int ISet, const char* Folder) {
  
  HBook1F(Hist->nhits,Form("nhits_%02i",ISet),"nhits",20,0,20,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int simres::book_histograms() {
  char folder_name[128];
  TFolder* fol;

  TH1::AddDirectory(0);

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fEvt[i] = new EventHist_t;
      book_event_histograms(fHist.fEvt[i],i,Form("%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book segment histograms
//-----------------------------------------------------------------------------
  int book_segment_histset[kNSegmentHistSets];
  for (int i=0; i<kNSegmentHistSets; i++) book_segment_histset[i] = 0;

  book_segment_histset[ 0] = 1;		// all
  book_segment_histset[ 4] = 1;		// segments with 4 hits
  book_segment_histset[ 5] = 1;		// all
  book_segment_histset[ 6] = 1;		// all
  book_segment_histset[ 7] = 1;		// all
  book_segment_histset[ 8] = 1;		// all
  book_segment_histset[ 9] = 1;		// all
  book_segment_histset[10] = 1;		// segments with 10 hits

  for (int i=0; i<kNSegmentHistSets; i++) {
    if (book_segment_histset[i] != 0) {
      sprintf(folder_name,"seg_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fSeg[i] = new SegmentHist_t;
      book_segment_histograms(fHist.fSeg[i],i,Form("%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book hit histograms
//-----------------------------------------------------------------------------
  int book_hit_histset[kNHitHistSets];
  for (int i=0; i<kNHitHistSets; i++) book_hit_histset[i] = 0;

  book_hit_histset[ 0] = 1;		// all
  book_hit_histset[ 4] = 1;		// all  4-hit segments
  book_hit_histset[ 5] = 1;		// all  5-hit segments
  book_hit_histset[ 6] = 1;		// all  6-hit segments
  book_hit_histset[ 7] = 1;		// all  7-hit segments
  book_hit_histset[ 8] = 1;		// all  8-hit segments
  book_hit_histset[ 9] = 1;		// all  9-hit segments
  book_hit_histset[10] = 1;		// all 10-hit segments

  book_hit_histset[40] = 1;		// hit#0 in a segment
  book_hit_histset[41] = 1;		// hit#1 in a segment
  book_hit_histset[42] = 1;		// hit#2 in a segment
  book_hit_histset[43] = 1;		// hit#3 in a segment

  book_hit_histset[50] = 1;		// hit#0 in a segment
  book_hit_histset[51] = 1;		// hit#1 in a segment
  book_hit_histset[52] = 1;		// hit#2 in a segment
  book_hit_histset[53] = 1;		// hit#3 in a segment
  book_hit_histset[54] = 1;		// hit#4 in a segment

  book_hit_histset[60] = 1;		// hit#0 in a segment
  book_hit_histset[61] = 1;		// hit#1 in a segment
  book_hit_histset[62] = 1;		// hit#2 in a segment
  book_hit_histset[63] = 1;		// hit#3 in a segment
  book_hit_histset[64] = 1;		// hit#4 in a segment
  book_hit_histset[65] = 1;		// hit#5 in a segment

  book_hit_histset[70] = 1;		// hit#0 in a segment
  book_hit_histset[71] = 1;		// hit#1 in a segment
  book_hit_histset[72] = 1;		// hit#2 in a segment
  book_hit_histset[73] = 1;		// hit#3 in a segment
  book_hit_histset[74] = 1;		// hit#4 in a segment
  book_hit_histset[75] = 1;		// hit#5 in a segment
  book_hit_histset[76] = 1;		// hit#6 in a segment

  book_hit_histset[80] = 1;		// hit#0 in a segment
  book_hit_histset[81] = 1;		// hit#1 in a segment
  book_hit_histset[82] = 1;		// hit#2 in a segment
  book_hit_histset[83] = 1;		// hit#3 in a segment
  book_hit_histset[84] = 1;		// hit#4 in a segment
  book_hit_histset[85] = 1;		// hit#5 in a segment
  book_hit_histset[86] = 1;		// hit#6 in a segment
  book_hit_histset[87] = 1;		// hit#7 in a segment

  book_hit_histset[90] = 1;		// hit#0 in a segment
  book_hit_histset[91] = 1;		// hit#1 in a segment
  book_hit_histset[92] = 1;		// hit#2 in a segment
  book_hit_histset[93] = 1;		// hit#3 in a segment
  book_hit_histset[94] = 1;		// hit#4 in a segment
  book_hit_histset[95] = 1;		// hit#5 in a segment
  book_hit_histset[96] = 1;		// hit#6 in a segment
  book_hit_histset[97] = 1;		// hit#7 in a segment
  book_hit_histset[98] = 1;		// hit#8 in a segment

  book_hit_histset[100] = 1;		// hit#0 in a segment
  book_hit_histset[101] = 1;		// hit#1 in a segment
  book_hit_histset[102] = 1;		// hit#2 in a segment
  book_hit_histset[103] = 1;		// hit#3 in a segment
  book_hit_histset[104] = 1;		// hit#4 in a segment
  book_hit_histset[105] = 1;		// hit#5 in a segment
  book_hit_histset[106] = 1;		// hit#6 in a segment
  book_hit_histset[107] = 1;		// hit#7 in a segment
  book_hit_histset[108] = 1;		// hit#8 in a segment
  book_hit_histset[109] = 1;		// hit#9 in a segment


  for (int i=0; i<kNHitHistSets; i++) {
    if (book_hit_histset[i] != 0) {
      sprintf(folder_name,"hit_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fHit[i] = new HitHist_t;
      book_hit_histograms(fHist.fHit[i],i,Form("%s",folder_name));
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int simres::fill_event_histograms(EventHist_t* Hist) {
  Hist->nhits->Fill(fSegment.nHits());
  return 0;
}

//-----------------------------------------------------------------------------
// ## TODO
//-----------------------------------------------------------------------------
int simres::fill_hit_histograms(HitHist_t* Hist, HitData_t* Hd) {
  Hist->r->Fill(Hd->r);
  Hist->doca->Fill(Hd->doca);

  Hist->dr->Fill(Hd->dr);
  Hist->drho->Fill(Hd->drho);

  Hist->dr_vs_r->Fill(Hd->r,Hd->dr);
  Hist->drho_vs_r->Fill(Hd->r,Hd->drho);

  Hist->dr_vs_is->Fill(Hd->is,Hd->dr);
  Hist->drho_vs_is->Fill(Hd->is,Hd->drho);

  return 0;
}

//-----------------------------------------------------------------------------
int simres::fill_segment_histograms(SegmentHist_t* Hist) {
  Hist->dy->Fill(fEvent.dy);
  Hist->ddy->Fill(fEvent.ddy);
  Hist->dt0->Fill(fEvent.dt0);
  return 0;
}

//-----------------------------------------------------------------------------
int simres::fill_histograms(Hist_t* Hist) {

//-----------------------------------------------------------------------------
// fill event-level histograms
//-----------------------------------------------------------------------------
  fill_event_histograms(Hist->fEvt[0]);

//-----------------------------------------------------------------------------
// fill segment histograms
//-----------------------------------------------------------------------------
  int nhits = fSegment.nHits();

  fEvent.dy    = fSegment.y0_rec()-fY0;
  fEvent.ddy   = fSegment.DyDx()-fDyDx;
  fEvent.dt0   = fSegment.t0_rec()-fT0;
  fEvent.nhits = nhits;


  for (int nh=4; nh<11; nh++) {
    if (nhits == nh) {
      fill_segment_histograms(Hist->fSeg[nh]);
      break;
    }
  }
//-----------------------------------------------------------------------------
// fill hit histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; ++i) {
    SegmentHit* sgh = fSegment.Hit(i);
    const mu2e::ComboHit* ch = sgh->ComboHit();
    fHitData.index = i;
    fHitData.is    = ch->strawId().straw();
    fHitData.r     = fSegment.R(sgh);
    fHitData.doca  = fSegment.Doca(sgh);
    fHitData.drho  = fabs(fHitData.doca)-fHitData.r;
    fHitData.dr    = fSegment.Dr(sgh);

    fill_hit_histograms(Hist->fHit[   0],&fHitData);
    if ((nhits >= 4) and (nhits <= 10)) {
//-----------------------------------------------------------------------------
// split by hit index within the segment
// for nhits= 4: fHit[ 4] and fHit[ 40- 43]
// for nhits=10: fHit[10] and fHit[100-109]
//-----------------------------------------------------------------------------
      fill_hit_histograms(Hist->fHit[ nhits    ],&fHitData);
      fill_hit_histograms(Hist->fHit[10*nhits+i],&fHitData);
    }
  }

  if (fabs(fEvent.dy) > fMaxGoodDy) {
    print_(std::format("dy:{:10.4f} dt0:{:10.4f}\n",fEvent.dy,fEvent.dt0));
  }
  return 0;
}

//-----------------------------------------------------------------------------
// need a mode in which would just use fixed line parameters (for debugging)
//-----------------------------------------------------------------------------
int simres::simulate_line_parameters() {
  int rc(0);

  fT0      = 0.;                        // trn3.Gaus(0,10); // turn off the T0 smearing to begin with
  fDyDx    = 0.4;                       // trn3.Gaus(0.4,0.1);
  fX0      = 0.;
  double x = 0.; // (2*fTrn3.Rndm()-1)*10.0;    // [-10,10] mm should cover pretty much the full range
  fY0      = -x/fDyDx;                  // y(x=0)

  if (fBeamMode == 1) {
                                        // these parameters lead to odd reco
    fY0   = -22.94303316;
    fDyDx = 0.4;
    fT0   = 0.;
  }

  if (fDebugMode & 0x8) {
    print_(std::format("fY0:{:14.8f} fDyDx:{:14.8f} fT0:{:14.8f}\n",fY0,fDyDx,fT0));
  }

  fNx  = 1./sqrt(1+fDyDx*fDyDx);
  fNy  = fDyDx*fNx;
  fNux = -fNy;
  fNuy = fNx;

  if (fDebugMode & 0x1) {
    print_(std::format("nx:{:10.6f} ny:{:10.6f} nux:{:10.6f} nuy:{:10.6f}\n",fNx,fNy,fNux,fNuy));
  }

  return rc;
}

//-----------------------------------------------------------------------------
int simres::simulate_hits() {
  int rc(0);
  for (int is=0; is<96; is++) {
    if (fDebugMode & 0x1) {
      std::cout << __func__ << ":" << __LINE__ << ": ------------------------- is:" << is << std::endl;
    }
    Straw_t* s = Straw(is);
 
    double dr = (fY0-s->fY)*fNuy+(fX0-s->fX)*fNux; // distance (signed) from the track to the wire
 
    if (fDebugMode & 0x1) {
      print_(std::format("s->fX:{:12.6f} s->fY:{:12.6f} dr:{:12.6f}\n",s->fX,s->fY,dr));
    }
 
    if (fabs(dr) < fStrawRadius) {
      s->fDoca = fabs(dr);
      fPanel.fListOfHitStraws.push_back(s);
      if (fDebugMode & 0x2) {
        print_(std::format("-- hit in straw:{:2d} x0:{:10.3f}  y0:{:10.3f} doca:{:10.3f}\n",
                           is,s->fX,s->fY,s->fDoca));
      }
    }
    if (fDebugMode & 0x1) {
      std::cout <<  __func__ << ":" << __LINE__ << ": fNHits:" << fPanel.NHits() << std::endl;
    }
  }
  return rc;
}


//-----------------------------------------------------------------------------
int simres::process_event(int IEvent) {
  int rc(0);

  fEvent.number = IEvent;
  fPanel.Clear();
  fSegment.Clear();
//----------------------------------------------------------------------------
// 1) simulate the line parameters
//-----------------------------------------------------------------------------
  rc = simulate_line_parameters();
  if (rc != 0) return rc;
//-----------------------------------------------------------------------------
// 2) simulate hits
//-----------------------------------------------------------------------------
  rc = simulate_hits();
  if (rc != 0) return rc;
//-----------------------------------------------------------------------------
// 3) reconstruct the line
//    construct the segment - need to fake a list of ComboHits
//-----------------------------------------------------------------------------
  int nhits = fPanel.NHits();
  if (nhits > 2) {
    for (int i=0; i<nhits; ++i) {
      Straw_t* s          = HitStraw(i);
      double   drift_time = D2T(s->fDoca);
      
      mu2e::ComboHit* ch = new mu2e::ComboHit;
      ch->_sid                        = s->id;
      ch->_eend                       = mu2e::StrawEnd::cal;
//-----------------------------------------------------------------------------
// simulate drift time resolution - 2 ns
//-----------------------------------------------------------------------------
      // double sig_t0                  = fTrn3.Gaus(0.,fSigmaT0);
      // double tmeas                   = fT0+drift_time+sig_t0;
      double sig_t0 = 0;
      if (i == 1) sig_t0             = 10.; // fSigmaT0; // fTrn3.Gaus(0.,fSigmaT0);
      double tmeas                   = fT0+drift_time+sig_t0;
      if (tmeas < 0) tmeas = 0;
      ch->_etime[mu2e::StrawEnd::cal] = tmeas;
      ch->_ptime                      = 0;
      ch->_pos                        = CLHEP::Hep3Vector(0,s->fX,s->fY);
      
      fPanel.fListOfComboHits.emplace_back(ch);

      SegmentHit sgh(ch);
      sgh.fSigmaR = 0.125;            // 2ns*62.5 um/ns
      fSegment.fListOfHits.emplace_back(sgh);
    }
  }
//-----------------------------------------------------------------------------
// 4) make a segment a segment, don't rotate
//-----------------------------------------------------------------------------
  int do_transform(0);
  fSegment.InitHits(nullptr,-1,-1,do_transform);
//-----------------------------------------------------------------------------
// 5) at this point, should be ready to fit the segment
//-----------------------------------------------------------------------------
  TrkSegment::Par_t par;

  SegmentFit sfitter(&fSegment);
  sfitter.DefineDriftDirections();
  
  if (fDebugMode & 0x8) {
    fSegment.print();
  }
  
  int niter(6);
  int converged = sfitter.Fit(niter,0,nullptr,&par);
  if (fDebugMode & 0x1) {
    std::cout << std::format("-- event:{} converged:{}\n",IEvent,converged);
  }
  return rc;
}

//-----------------------------------------------------------------------------
int simres::run(int NEvents, int DebugMode) {
  int rc(0);

  init_geometry();

  book_histograms();
                                        // simulation per se

  for (int ievent=0; ievent<NEvents; ++ievent) {
//-----------------------------------------------------------------------------
// process events consists of two steps: simulate, then reconstruct
//-----------------------------------------------------------------------------
    rc = process_event(ievent);

    fill_histograms(&fHist);
  }
  return rc;
}

//-----------------------------------------------------------------------------
int simres::display_event() {
//-----------------------------------------------------------------------------
// plot the results after the last iteration, but before updating the points
//-----------------------------------------------------------------------------
  if (fDebugMode & 0x1) print_("START\n");

  int    nhits = fPanel.fListOfHitStraws.size();

  TCanvas* c = new TCanvas("c","c",1500,600);
  TH2F* h2 = new TH2F("h2","",500,-125,125,200,-50,50);
  h2->SetStats(kFALSE);
  h2->Draw();
  //  c->Range(-150,-600,150,60);
  c->cd();
//-----------------------------------------------------------------------------
// draw simulated trajectory
//-----------------------------------------------------------------------------
  double xmin(-150.), xmax(150.);
  double ymin = fY0+(xmin-fX0)*fDyDx;
  double ymax = fY0+(xmax-fX0)*fDyDx;
  TLine* l = new TLine(xmin,ymin,xmax,ymax);
  l->SetLineColor(kRed+2);
  l->Draw();
//-----------------------------------------------------------------------------
// draw reconstructed segment
//-----------------------------------------------------------------------------
  ymin = fSegment.y0_rec()+(xmin-fX0)*fSegment.DyDx();
  ymax = fSegment.y0_rec()+(xmax-fX0)*fSegment.DyDx();
  TLine* l_rec = new TLine(xmin,ymin,xmax,ymax);
  l_rec->SetLineColor(kBlue+2);
  l_rec->Draw();

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
                                        // simulated hit
    mu2e::ComboHit* ch = ComboHit(i);
    double          t  = ch->time();
    double          r  = SegmentHit::T2D(t);

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

  if (fDebugMode & 0x1) print_("-- END\n");
  return 0;

}


//_____________________________________________________________________________
void     simres::AddHistogram(TObject* hist, const char* FolderName) {
  TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
  fol->Add(hist); 
}

//_____________________________________________________________________________
void simres::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void simres::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const float* LowEdge,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,LowEdge);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void simres::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // this is introduced specifically for weighting the DIO Mu2e spectrum
  // book 1D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void simres::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const double* LowEdge,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,LowEdge);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void simres::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 Int_t Ny, Double_t YMin, Double_t YMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void simres::HProf(TProfile*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       Double_t YMin, Double_t YMax,
		       const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
  AddHistogram(Hist,FolderName);
}




//_____________________________________________________________________________
void simres::DeleteHistograms(TFolder* Folder) {
  // internal method...

  if (((long int) Folder) == -1) Folder = fFolder;

  TObject  *o1;

  TIter    it1(Folder->GetListOfFolders());

  while ((o1 = it1.Next())) {
    if (o1->InheritsFrom("TFolder")) {
      DeleteHistograms((TFolder*) o1);
    }
    else if (o1->InheritsFrom("TH1")) {
      Folder->Remove(o1);
      delete o1;
    }
  }
}

//_____________________________________________________________________________
int  simres::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

  TDirectory*  dir;
  TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
  Dir->cd();
  //  dir = new TDirectory(Folder->GetName(),Folder->GetName(),"");
  dir = Dir->mkdir(Folder->GetName(),Folder->GetTitle());
  dir->cd();

//   printf(" ------------------- Dir: %s, new dir: %s\n",
// 	 Dir->GetName(),dir->GetName());


  TIter  it(Folder->GetListOfFolders());
  while ((o = it.Next())) {
//     printf(" o->GetName, o->ClassName : %-20s %-20s\n",
// 	   o->GetName(),
// 	   o->ClassName());

    if (strcmp(o->ClassName(),"TFolder") == 0) {
      simres::SaveFolder((TFolder*) o, dir);
      //      dir->cd();
    }
    else if (! o->InheritsFrom("TStnModule")) {
      //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
      o->Write();
      //      gDirectory->GetListOfKeys()->Print();
    }
  }

  Dir->cd();
  return 0;
}

//_____________________________________________________________________________
void simres::SaveHist(const char* Filename) {
  // save histograms booked by all the modules into a file with the given name
  // Mode = 1: save folders

  TFile* f = new TFile(Filename,"recreate");

  simres::SaveFolder(fFolder,f);

  f->Close();
  delete f;
}


