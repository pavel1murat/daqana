//
#include <iostream>
#include <format>
#include "TH2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"
// #include "TGeoMatrix.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/obj/TrkSegment.hh"

#include "daqana/obj/SegmentFit.hh"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"

#include "daqana/test/test3.hh"
#include "daqana/test/station0_heptransform.hh"

int readDataFile(const char* Fn, std::vector<std::string>& Names, std::vector<const mu2e::ComboHit*>& Data, int Plane, int Panel);
void printData(const std::vector<std::string>& Names, std::vector<const mu2e::ComboHit*> Data);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
int test3::init_segment_geometry(TrkSegment* Seg) {
  int rc(0);
  
  int unique_panel = Seg->Plane()*6 + Seg->Panel();
  std::cout << "unique_panel:" << unique_panel << std::endl;
  
  HepTransformData_t* dt = &ht_data[unique_panel];
  CLHEP::Hep3Vector trn(dt->trn[0],dt->trn[1],dt->trn[2]);
  std::cout << "read translation:" << trn << std::endl;
  
  CLHEP::HepRotation r(dt->rot);
  mu2e::HepTransform ht(trn,r);
  std::cout << "read rotation:" << trn << std::endl;
  
  Seg->fTrkPanel = new mu2e::Panel();
  Seg->fTrkPanel->setPanelToDS(ht);
  
  std::cout << "-- transformation matrix:" << Seg->fTrkPanel->dsToPanel();
  return rc;
}

//-----------------------------------------------------------------------------
int test3::test_fit_line(const char* Fn, int Plane, int Panel, int NIter) {
  int          rc(0);
  double const epsilon(1.e-5);

  std::cout << std::format("---------------------------------- {}:{}: START\n",__func__,__LINE__);

  TrkSegment::fgDebugMode = 1;
  SegmentFit::fgDebugMode = 1;
//-----------------------------------------------------------------------------
// read combohits from a text file
//-----------------------------------------------------------------------------
  std::vector<std::string>     vnames;
  std::vector<const mu2e::ComboHit*> chhits;

  std::cout << "-- before readDataFile" << std::endl;
  readDataFile(Fn,vnames,chhits,Plane,Panel);
  printData(vnames,chhits);
  
  // std::cout << "-- after readDataFile" << std::endl;
  //                                       // sort by straw
  // std::sort(chhits.begin(),chhits.end(), [](const mu2e::ComboHit* a, const mu2e::ComboHit* b) {
  //   mu2e::StrawId sa(a->strawId()), sb(b->strawId());
  //   return sa.getStraw() < sb.getStraw();
  // });
  
  int nhits = chhits.size();
  std::cout << "-- hits sorted, nhits:" << nhits << std::endl;
  if (nhits < 2) return -1;
//-----------------------------------------------------------------------------
// step 1: initialize the segment and the hit drift signs
//-----------------------------------------------------------------------------
  if (mTs) delete mTs;
  mTs = new TrkSegment(Plane,Panel);
  init_segment_geometry(mTs);
//-----------------------------------------------------------------------------
// InitHits also defines two seed hits - can draw a tangent line at this point
//-----------------------------------------------------------------------------
  mTs->InitHits(&chhits);
//-----------------------------------------------------------------------------
  if (sfitter != nullptr) delete sfitter;
  sfitter = new SegmentFit(mTs);

  mTs->print("-- segment after creating sfitter");

                                        // this needs to be done just once
  sfitter->DefineDriftDirections();

  mTs->print("-- segment after SegmentFit::DefineDriftDirections");
//-----------------------------------------------------------------------------  
// ready to iterate, before that print segment and pop up a canvas
//-----------------------------------------------------------------------------
  mCanvas = new TCanvas("c","c",1000,1000);

  // set range of the 2D plot

  double xm = (mTs->points[nhits-1].x+mTs->points[0].x)/2;
  double dx = (mTs->points[nhits-1].x-mTs->points[0].x);
  if (dx < 10) dx = 10;

  mH2 = new TH2F(Form("plane_%02d_panel_%d",Plane,Panel),"",
                       1000,xm-dx/2-5,xm+dx/2+5,1000,-dx/2-5,dx/2+5);
  mH2->SetStats(0);
  mH2->SetTitle(Form("data:%s plane:%i panel:%i",Fn,Plane,Panel));
  mH2->Draw();
//-----------------------------------------------------------------------------
// iterations, at this point all drift signs are defined...
// pin=nullptr: use fSegment->fPar to start
//-----------------------------------------------------------------------------
  Par_t par;
  int converged = sfitter->Fit(NIter,0, nullptr, &par);
  mTs->print(Form("-- after SegmentFit::Fit converged:%i",converged));

  mH2->Draw();
  sfitter->DisplaySegment();

  return rc;
}
