//
#include <iostream>
#include <format>
#include "TH2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TGeoMatrix.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/obj/obj/TrkSegment.hh"

#include "daqana/obj/obj/SegmentFit.hh"
#include "daqana/scripts/panel_geometry.C"
#include "daqana/scripts/read_combohits_2.C"

//-----------------------------------------------------------------------------
// global variables
//-----------------------------------------------------------------------------
TCanvas*    mCanvas;
TH2F*       mH2;
TrkSegment* mTs (nullptr);
SegmentFit* sfitter;

//-----------------------------------------------------------------------------
int init_segment_geometry(TrkSegment* Seg) {
  int rc(0);

  int loc              = 12*Seg->Plane()+2*Seg->Panel();
  PanelGeometry_t* pg0 = &PanelGeometry[loc  ]; // layer 0
  PanelGeometry_t* pg1 = &PanelGeometry[loc+1]; // layer 1

  double phix(0), phiy(90), phiz(0), thetax(90), thetay(90), thetaz(0);
  if (Seg->Panel() < 6) { 
    phix   = atan2(pg0->yc,pg0->xc)*180./M_PI;
    phiy   = phix+90;
    phiz   =  0;
    thetax = 90;
    thetay = 90;
    thetaz =  0;
  }

  TGeoRotation* r  = new TGeoRotation  ("a",thetax,phix,thetay,phiy,thetaz,phiz);
  Seg->fCombiTrans = new TGeoCombiTrans(Form("plane_%02i_%i",Seg->Plane(),Seg->Panel()),0,0,0,r);
  std::cout << "-- transformation matrix\n";
  Seg->fCombiTrans->Print();
  return rc;
}

//-----------------------------------------------------------------------------
int test_fit_line(const char* Fn, int Plane, int Panel, int NIter=0) {
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
