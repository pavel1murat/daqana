///////////////////////////////////////////////////////////////////////////////
// SimRes: reconstruct the segment of the simulated hits (see SimHits) with SegmentFit,
//         put the segment into the event, and book histograms of the residuals
//         the segment hits refer to the ComboHits of the event by their index
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <format>
#include <iostream>
#include <memory>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/alg/SegmentFit.hh"
#include "daqana/art_obj/TrkSegment.hh"
#include "SimPanel.hh"

namespace mu2e {

class SimRes : public art::EDProducer {
public:
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int>          debugMode      {Name("debugMode"      ), Comment("debug bits, see SimPanel.hh")};
    fhicl::Atom<double>       strawRadius    {Name("strawRadius"    ), Comment("straw radius, mm")};
    fhicl::Atom<double>       vDrift         {Name("vDrift"         ), Comment("drift velocity, mm/ns")};
    fhicl::Atom<int>          beamMode       {Name("beamMode"       ), Comment("1: line parameters of the odd-reco study, as in SimHits")};
    fhicl::Atom<double>       maxGoodDy      {Name("maxGoodDy"      ), Comment("print the events with |dy| above this value, mm")};
    fhicl::Atom<int>          nFitIterations {Name("nFitIterations" ), Comment("number of iterations of the segment fit")};
    fhicl::Atom<bool>         displayEvent   {Name("displayEvent"   ), Comment("draw each event")};
    fhicl::Atom<art::InputTag> chCollTag    {Name("chCollTag"     ), Comment("the ComboHits of the simulated hits")};
    fhicl::Sequence<double> strawOffsetX{Name("strawOffsetX"), Comment("misalignment of the straws in X, 96 values, mm")};
    fhicl::Sequence<double> strawOffsetY{Name("strawOffsetY"), Comment("misalignment of the straws in Y, 96 values, mm")};
  };
  using Parameters = art::EDProducer::Table<Config>;

  explicit SimRes(Parameters const& config);

  void beginJob() override;
  void produce(art::Event& event) override;

private:
  enum {
    kNHitHistSets     = 1000,           // used up to 109
    kNSegmentHistSets =   10,
    kNEventHistSets   =   10,
  };

  struct HitHist_t {
    TH1F* r;                            // drift radius
    TH1F* doca;                         // [signed] track-to-wire distance
    TH1F* dr;
    TH1F* drho;
    TH2F* dr_vs_r;
    TH2F* drho_vs_r;
    TH2F* dr_vs_is;
    TH2F* drho_vs_is;
  };

  struct EventHist_t {
    TH1F* nhits;
  };

  struct SegmentHist_t {
    TH1F* dy;
    TH1F* ddy;
    TH1F* dt0;
  };

  struct Hist_t {
    EventHist_t*   fEvt[kNEventHistSets  ] = {};
    SegmentHist_t* fSeg[kNSegmentHistSets] = {};
    HitHist_t*     fHit[kNHitHistSets    ] = {};
  };

  struct HitData_t {
    int    index;                       // hit index within the segment
    int    is;                          // straw number
    double r;                           // drift radius
    double dr;                          // residual
    double drho;                        // rdrift - track_doca
    double doca;
  };

  struct EventData_t {
    int    nhits;
    double dy;
    double ddy;
    double dt0;
  };

  int                       fDebugMode;
  double                    fStrawRadius;
  double                    fVDrift;
  int                       fBeamMode;
  double                    fMaxGoodDy;
  int                       fNFitIterations;
  bool                      fDisplayEvent;
  std::vector<double>       fStrawOffsetX;
  std::vector<double>       fStrawOffsetY;
  art::InputTag             fComboHitTag;
  int                       fEventNumber;

  simpanel::Line_t          fLine;            // the simulated line, the truth for the residuals
  TrkSegment                fSegment;         // reconstructed segment, its hits refer to the ComboHits
  Hist_t                    fHist;
  HitData_t                 fHitData;
  EventData_t               fEvent;

  void  print_(const std::string& Message);

  void  book_event_histograms  (EventHist_t*   Hist, art::TFileDirectory const& Dir, int ISet);
  void  book_hit_histograms    (HitHist_t*     Hist, art::TFileDirectory const& Dir, int ISet);
  void  book_segment_histograms(SegmentHist_t* Hist, art::TFileDirectory const& Dir, int ISet);
  void  book_histograms        (art::TFileService* Tfs);

  void  fill_event_histograms  (EventHist_t*   Hist);
  void  fill_hit_histograms    (HitHist_t*     Hist, HitData_t* Hd);
  void  fill_segment_histograms(SegmentHist_t* Hist);
  void  fill_histograms        (Hist_t*        Hist);

  void  display_event(mu2e::ComboHitCollection const& Hits);
};

//-----------------------------------------------------------------------------
SimRes::SimRes(Parameters const& config) : art::EDProducer(config) {
  fDebugMode      = config().debugMode();
  fStrawRadius    = config().strawRadius();
  fVDrift         = config().vDrift();
  fBeamMode       = config().beamMode();
  fMaxGoodDy      = config().maxGoodDy();
  fNFitIterations = config().nFitIterations();
  fDisplayEvent   = config().displayEvent();
  fComboHitTag    = config().chCollTag();
  fStrawOffsetX   = config().strawOffsetX();
  fStrawOffsetY   = config().strawOffsetY();
  fEventNumber    = -1;

  SegmentHit::SetVDrift (fVDrift);
  SegmentHit::SetTOffset(0);

  produces<TrkSegmentCollection>();
}

//-----------------------------------------------------------------------------
void SimRes::print_(const std::string& Message) {
  std::cout << std::format("event:{:6d} : {}", fEventNumber, Message);
}

//-----------------------------------------------------------------------------
void SimRes::book_event_histograms(EventHist_t* Hist, art::TFileDirectory const& Dir, int ISet) {
  Hist->nhits = Dir.make<TH1F>(std::format("nhits_{:02}",ISet).c_str(),"nhits",20,0,20);
}

//-----------------------------------------------------------------------------
void SimRes::book_hit_histograms(HitHist_t* Hist, art::TFileDirectory const& Dir, int ISet) {
  Hist->r         = Dir.make<TH1F>(std::format("r_{:02}"          ,ISet).c_str(),"rdrift, mm", 100,  0,  5);
  Hist->doca      = Dir.make<TH1F>(std::format("doca_{:02}"       ,ISet).c_str(),"doca, mm"  , 200, -5,  5);
  Hist->dr        = Dir.make<TH1F>(std::format("dr_{:02}"         ,ISet).c_str(),"dr, mm"    , 500, -2.5, 2.5);
  Hist->drho      = Dir.make<TH1F>(std::format("drho_{:02}"       ,ISet).c_str(),"drho, mm"  , 500, -2.5, 2.5);
  Hist->dr_vs_r   = Dir.make<TH2F>(std::format("dr_vs_r_{:02}"    ,ISet).c_str(),"dr_vs_r"   , 100, 0, 5, 250, -2.5, 2.5);
  Hist->drho_vs_r = Dir.make<TH2F>(std::format("drho_vs_r_{:02}"  ,ISet).c_str(),"drho_vs_r" , 100, 0, 5, 250, -2.5, 2.5);
  Hist->dr_vs_is  = Dir.make<TH2F>(std::format("dr_vs_is_{:02}"   ,ISet).c_str(),"dr_vs_is"  , 100, 0, 100, 250, -2.5, 2.5);
  Hist->drho_vs_is= Dir.make<TH2F>(std::format("drho_vs_is_{:02}" ,ISet).c_str(),"drho_vs_is", 100, 0, 100, 250, -2.5, 2.5);
}

//-----------------------------------------------------------------------------
void SimRes::book_segment_histograms(SegmentHist_t* Hist, art::TFileDirectory const& Dir, int ISet) {
  Hist->dy  = Dir.make<TH1F>(std::format("dy_{}" ,ISet).c_str(), "Y0(reco)-Y0(sim)"    , 100, -5, 5);
  Hist->ddy = Dir.make<TH1F>(std::format("ddy_{}",ISet).c_str(), "DyDx(reco)-DyDx(sim)", 200, -0.1, 0.1);
  Hist->dt0 = Dir.make<TH1F>(std::format("dt0_{}",ISet).c_str(), "T0(reco)-T0(sim)"    , 200, -10, 10);
}

//-----------------------------------------------------------------------------
// the sets booked: event set 0; segment sets 0 and 4-10; hit sets 0, 4-10, and the
// sets split by the hit index within the segment, see fill_histograms
//-----------------------------------------------------------------------------
void SimRes::book_histograms(art::TFileService* Tfs) {

  TH1::AddDirectory(kFALSE);

  art::TFileDirectory evt_dir = Tfs->mkdir("evt_0");
  fHist.fEvt[0] = new EventHist_t;
  book_event_histograms(fHist.fEvt[0],evt_dir,0);

  for (int i=0; i<kNSegmentHistSets; i++) {
    if ((i != 0) and ((i < 4) or (i > 10))) continue;
    art::TFileDirectory dir = Tfs->mkdir(std::format("seg_{}",i));
    fHist.fSeg[i] = new SegmentHist_t;
    book_segment_histograms(fHist.fSeg[i],dir,i);
  }

  for (int i=0; i<kNHitHistSets; i++) {
    bool booked = (i == 0) or ((i >= 4) and (i <= 10)) or
      ((i >= 40) and (i <= 43)) or ((i >= 50) and (i <= 54)) or ((i >= 60) and (i <= 65)) or
      ((i >= 70) and (i <= 76)) or ((i >= 80) and (i <= 87)) or ((i >= 90) and (i <= 98)) or
      ((i >= 100) and (i <= 109));
    if (not booked) continue;
    art::TFileDirectory dir = Tfs->mkdir(std::format("hit_{}",i));
    fHist.fHit[i] = new HitHist_t;
    book_hit_histograms(fHist.fHit[i],dir,i);
  }
}

//-----------------------------------------------------------------------------
void SimRes::fill_event_histograms(EventHist_t* Hist) {
  Hist->nhits->Fill(fSegment.NHitsWithComboHit());
}

//-----------------------------------------------------------------------------
void SimRes::fill_hit_histograms(HitHist_t* Hist, HitData_t* Hd) {
  Hist->r->Fill(Hd->r);
  Hist->doca->Fill(Hd->doca);

  Hist->dr->Fill(Hd->dr);
  Hist->drho->Fill(Hd->drho);

  Hist->dr_vs_r->Fill(Hd->r,Hd->dr);
  Hist->drho_vs_r->Fill(Hd->r,Hd->drho);

  Hist->dr_vs_is->Fill(Hd->is,Hd->dr);
  Hist->drho_vs_is->Fill(Hd->is,Hd->drho);
}

//-----------------------------------------------------------------------------
void SimRes::fill_segment_histograms(SegmentHist_t* Hist) {
  Hist->dy->Fill(fEvent.dy);
  Hist->ddy->Fill(fEvent.ddy);
  Hist->dt0->Fill(fEvent.dt0);
}

//-----------------------------------------------------------------------------
void SimRes::fill_histograms(Hist_t* Hist) {

  fill_event_histograms(Hist->fEvt[0]);

  int n_hits = fSegment.NHitsWithComboHit();

  fEvent.dy    = fSegment.y0_rec()-fLine.fY0;
  fEvent.ddy   = fSegment.DyDx()-fLine.fDyDx;
  fEvent.dt0   = fSegment.t0_rec()-fLine.fT0;
  fEvent.nhits = n_hits;

  for (int nh=4; nh<11; nh++) {
    if (n_hits == nh) {
      fill_segment_histograms(Hist->fSeg[nh]);
      break;
    }
  }

  int nstraws = fSegment.NStraws();
  for (int i=0; i<nstraws; ++i) {
    SegmentHit* sgh = fSegment.Hit(i);
    const mu2e::ComboHit* ch = sgh->ComboHit();
    if (ch == nullptr) continue;
    fHitData.index = i;
    fHitData.is    = ch->strawId().straw();
    fHitData.r     = fSegment.R(sgh);
    fHitData.doca  = fSegment.Doca(sgh);
    fHitData.drho  = fabs(fHitData.doca)-fHitData.r;
    fHitData.dr    = fSegment.Dr(sgh);

    fill_hit_histograms(Hist->fHit[   0],&fHitData);
    if ((n_hits >= 4) and (n_hits <= 10)) {
//-----------------------------------------------------------------------------
// split by hit index within the segment
// for nhits= 4: fHit[ 4] and fHit[ 40- 43]
// for nhits=10: fHit[10] and fHit[100-109]
//-----------------------------------------------------------------------------
      fill_hit_histograms(Hist->fHit[ n_hits    ],&fHitData);
      fill_hit_histograms(Hist->fHit[10*n_hits+i],&fHitData);
    }
  }

  if (fabs(fEvent.dy) > fMaxGoodDy) {
    print_(std::format("dy:{:10.4f} dt0:{:10.4f}\n",fEvent.dy,fEvent.dt0));
  }
}

//-----------------------------------------------------------------------------
void SimRes::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  book_histograms(&*tfs);
}

//-----------------------------------------------------------------------------
// one segment per event: fit the ComboHits of the event
//-----------------------------------------------------------------------------
void SimRes::produce(art::Event& event) {
  fEventNumber = event.event();
  fSegment.Clear();
  fLine = simpanel::line(fBeamMode,fDebugMode);

  auto const chH = event.getValidHandle<mu2e::ComboHitCollection>(fComboHitTag);
  int nhits = chH->size();

  auto segments = std::make_unique<TrkSegmentCollection>();
  if (nhits > 2) {
    for (int i=0; i<nhits; ++i) {
      // model plane=0 panel=0
      SegmentHit sgh(i,StrawId(0,0,i));

      // assume ~2ns resolution in drift time :2ns*62.5 um/ns
      sgh.fSigmaR = 0.125;

      fSegment.fListOfHits.emplace_back(sgh);
    }
//-----------------------------------------------------------------------------
// make the segment, don't rotate
//-----------------------------------------------------------------------------
    fSegment.setComboHits(chH.product());
    fSegment.InitHits(0);

    TrkSegment::Par_t par;
    SegmentFit sfitter(&fSegment);
    sfitter.DefineDriftDirections();

    if (fDebugMode & 0x8) {
      fSegment.print();
    }

    int converged = sfitter.Fit(fNFitIterations,0,nullptr,&par);
    if (fDebugMode & 0x1) {
      print_(std::format("converged:{}\n",converged));
    }

    fill_histograms(&fHist);
    TrkSegment copy = fSegment;
    copy.setComboHits(nullptr);        // the pointer is valid only while the event is processed
    segments->push_back(copy);
  }
  event.put(std::move(segments));

  if (fDisplayEvent) display_event(*chH);
}

//-----------------------------------------------------------------------------
void SimRes::display_event(mu2e::ComboHitCollection const& Hits) {

  auto const straws    = simpanel::straws(fDebugMode,fStrawOffsetX,fStrawOffsetY);
  auto const hitStraws = simpanel::hitStraws(straws,fLine,fStrawRadius,fDebugMode);
  int nhits = hitStraws.size();

  TCanvas* c = new TCanvas("c","c",1500,600);
  TH2F* h2 = new TH2F("h2","",500,-125,125,200,-50,50);
  h2->SetStats(kFALSE);
  h2->Draw();
  c->cd();
//-----------------------------------------------------------------------------
// draw simulated trajectory
//-----------------------------------------------------------------------------
  double xmin(-150.), xmax(150.);
  double ymin = fLine.fY0+(xmin-fLine.fX0)*fLine.fDyDx;
  double ymax = fLine.fY0+(xmax-fLine.fX0)*fLine.fDyDx;
  TLine* l = new TLine(xmin,ymin,xmax,ymax);
  l->SetLineColor(kRed+2);
  l->Draw();
//-----------------------------------------------------------------------------
// draw reconstructed segment
//-----------------------------------------------------------------------------
  ymin = fSegment.y0_rec()+(xmin-fLine.fX0)*fSegment.DyDx();
  ymax = fSegment.y0_rec()+(xmax-fLine.fX0)*fSegment.DyDx();
  TLine* l_rec = new TLine(xmin,ymin,xmax,ymax);
  l_rec->SetLineColor(kBlue+2);
  l_rec->Draw();

  for (auto const& s : straws) {
    TEllipse* e = new TEllipse(s.fX,s.fY,fStrawRadius,0,0,360);
    e->SetLineColor(kBlack);
    e->SetFillStyle(0);
    e->Draw();
  }
//-----------------------------------------------------------------------------
// draw hits: the straw, and the simulated hit with the radius from its time
//-----------------------------------------------------------------------------
  for (int i=0; i<nhits; i++) {
    simpanel::Straw_t const& s = hitStraws[i];

    TEllipse* e = new TEllipse(s.fX,s.fY,fStrawRadius,0,0,360);
    e->SetLineColor(kRed);
    e->SetFillStyle(0);
    e->Draw();

    double t = Hits[i].time();
    double r = SegmentHit::T2D(t);

    TEllipse* e2 = new TEllipse(s.fX,s.fY,r,0,0,360);
    e2->SetFillColor(kBlue-10);
    e2->SetFillStyle(3003);
    e2->Draw();
  }

  gPad->Modified();
  gPad->Update();
}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimRes)
