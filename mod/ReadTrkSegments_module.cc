///////////////////////////////////////////////////////////////////////////////
// ReadTrkSegments: test module, reads a collection of TrkSegments from the input file
//                  and prints the segments with their hits; the hits are resolved
//                  through their index in the ComboHitCollection of the file.
//
// histograms: booked in beginRun, one directory per run, run_<N>, with a set of
//             histograms per event (evt_<i>) and per segment (seg_<i>) for each set i.
//             The sets are switched on with the booking flags in book_histograms;
//             they are filled in analyze.
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "cetlib_except/exception.h"
#include "TH1.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "daqana/art_obj/TrkSegment.hh"

namespace mu2e {

class ReadTrkSegments : public art::EDAnalyzer {
public:
  enum {
    kNEffHistSets = 100,
    kNEvtHistSets = 100,
    kNSegHistSets = 100,
  };

//-----------------------------------------------------------------------------
// FCL parameters
//-----------------------------------------------------------------------------
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> tsCollTag{Name("tsCollTag"), Comment("the collection of TrkSegments to read")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  art::InputTag fTsCollTag;
//-----------------------------------------------------------------------------
// parameters of an event and of a segment, to be used by the histograms
//-----------------------------------------------------------------------------
  struct EffPar_t {
    int  ip216;                         // offline panel number in the range [0,215]
  };

  struct EvtPar_t {
    int  nseg;
  };

  struct SegPar_t {
    float x;
  };

//-----------------------------------------------------------------------------
// histogram records
//-----------------------------------------------------------------------------
  struct EffHist_t {
    TH1F* h_num  [216];                          // number of segments per event
    TH1F* h_den  [216];                          // number of segments per event
    TH1F* h_r_num[216];
    TH1F* h_r_den[216];
  };

  struct EvtHist_t {
    TH1F* h_nseg;                          // number of segments per event
  };

  struct SegHist_t {
    TH1F* h_nhits;                         // number of hits with a ComboHit per segment
    TH1F* h_chi2d;                         // chi2/dof per segment
  };

  struct Hist_t {
    EffHist_t* eff[kNEffHistSets];
    EvtHist_t* evt[kNEvtHistSets];
    SegHist_t* seg[kNSegHistSets];
  };
//-----------------------------------------------------------------------------
// oter parameters
//-----------------------------------------------------------------------------
  int                         fRunNumber;
  Hist_t*                     fHist;

  const TrkSegmentCollection* fTsColl;  // the segments of the current event
  int                         fNTs;     // number of track segments
  const ComboHitCollection*   fChColl;  // the ComboHits of the current event

  std::vector<EvtPar_t>       fEvtPar;  // event parameters, for the event histograms
  std::vector<SegPar_t>       fSegPar;  // segment parameters, one per segment of the current event

  std::vector<TrkSegment>     clone_;
//-----------------------------------------------------------------------------
// other parameters
//-----------------------------------------------------------------------------
  ProditionsHandle<Tracker>            fAlignedTracker;
  ProditionsHandle<TrackerStatus>      fTrackerStatusH;
  std::shared_ptr<const TrackerStatus> fTrackerStatus;
  const Tracker*                       fTracker;


  int                                  debugBit_[100]{};

  explicit ReadTrkSegments(Parameters const& config);

  void beginRun(art::Run const& ArtRun) override;
  void analyze(art::Event const& event) override;

  art::InputTag comboHitTagOf(art::ValidHandle<TrkSegmentCollection> const& SegH, art::Event const& Event) const;

  void book_eff_histograms(EffHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_evt_histograms(EvtHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_seg_histograms(SegHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_histograms(int RunNumber);

  void fill_eff_histograms(EffHist_t* Hist, SegmentHit* Hit, EffPar_t* EffP);
  void fill_evt_histograms(EvtHist_t* Hist, EvtPar_t* EvtP);
  void fill_seg_histograms(SegHist_t* Hist, TrkSegment const* Seg, SegPar_t* SegP);
  void fill_histograms();
};

//-----------------------------------------------------------------------------
ReadTrkSegments::ReadTrkSegments(Parameters const& config) : art::EDAnalyzer(config), fRunNumber(-1), fHist(nullptr),
  fTsColl(nullptr), fNTs(0), fChColl(nullptr) {
  fTsCollTag = config().tsCollTag();
  consumes<TrkSegmentCollection>(fTsCollTag);
                                        // the ComboHit collection is found from the provenance, see comboHitTagOf
  consumesMany<mu2e::ComboHitCollection>();
}

//-----------------------------------------------------------------------------
// the ComboHit collection the segments were made from: the parents of the segment collection
// are the products read by its producer, the ComboHit collection is one of them
//-----------------------------------------------------------------------------
art::InputTag ReadTrkSegments::comboHitTagOf(art::ValidHandle<TrkSegmentCollection> const& SegH, art::Event const& Event) const {
  art::TypeID const chType(typeid(mu2e::ComboHitCollection));
  for (art::ProductID const& pid : SegH.provenance()->parentage().parents()) {
    auto const desc = Event.getProductDescription(pid);
    if (desc and (desc->friendlyClassName() == chType.friendlyClassName())) {
      return art::InputTag(desc->moduleLabel(),desc->productInstanceName(),desc->processName());
    }
  }
  throw cet::exception("ReadTrkSegments") << "no ComboHitCollection among the parents of " << fTsCollTag << "\n";
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::book_eff_histograms(EffHist_t* Hist, int RunNumber, art::TFileDirectory* Dir) {

  std::string prefix = std::format("run number %6d",RunNumber);

  for (int i=0; i<216; i++) {
    Hist->h_num[i] = Dir->make<TH1F>(Form("num_%03d",i),Form("%s: panel %03d (num)",prefix.data(),i),100,0,100);
    Hist->h_den[i] = Dir->make<TH1F>(Form("den_%03d",i),Form("%s: panel %03d (den)",prefix.data(),i),100,0,100);

    Hist->h_r_num[i] = Dir->make<TH1F>(Form("r_num_%03d",i),Form("%s: panel %03d (num)",prefix.data(),i),50,0,2.5);
    Hist->h_r_den[i] = Dir->make<TH1F>(Form("r_den_%03d",i),Form("%s: panel %03d (den)",prefix.data(),i),50,0,2.5);
  }
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::book_evt_histograms(EvtHist_t* Hist, int RunNumber, art::TFileDirectory* Dir) {
  Hist->h_nseg = Dir->make<TH1F>("nseg","number of segments per event",100,0,100);
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::book_seg_histograms(SegHist_t* Hist, int RunNumber, art::TFileDirectory* Dir) {
  Hist->h_nhits = Dir->make<TH1F>("nhits","number of hits with a ComboHit per segment",60,0,60);
  Hist->h_chi2d = Dir->make<TH1F>("chi2d","chi2/dof per segment",100,0,50);
}

//-----------------------------------------------------------------------------
// the sets switched on: event set 0 and segment set 0; the other sets are booked by
// setting their flag to 1, and then filled in fill_histograms
//-----------------------------------------------------------------------------
void ReadTrkSegments::book_histograms(int RunNumber) {

  art::ServiceHandle<art::TFileService> tfs;
  TH1::AddDirectory(kFALSE);

  art::TFileDirectory run_dir = tfs->mkdir(std::format("run_{}",RunNumber));

  fHist = new Hist_t();
  for (int i=0; i<kNEvtHistSets; i++) fHist->evt[i] = nullptr;
  for (int i=0; i<kNSegHistSets; i++) fHist->seg[i] = nullptr;

  
//-----------------------------------------------------------------------------
// book efficiency histograms
//-----------------------------------------------------------------------------
  int book_eff_histset[kNEffHistSets];
  for (int i=0; i<kNEffHistSets; i++) book_eff_histset[i] = 0;

  book_eff_histset[0] = 1;              // all events
  book_eff_histset[1] = 0;              // reserved

  for (int i=0; i<kNEffHistSets; i++) {
    if (book_eff_histset[i] != 0) {
      art::TFileDirectory dir = run_dir.mkdir(std::format("eff_{}",i));
      fHist->eff[i] = new EffHist_t;
      book_eff_histograms(fHist->eff[i],RunNumber,&dir);
    }
  }

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_evt_histset[kNEvtHistSets];
  for (int i=0; i<kNEvtHistSets; i++) book_evt_histset[i] = 0;

  book_evt_histset[0] = 1;              // all events
  book_evt_histset[1] = 0;              // reserved

  for (int i=0; i<kNEvtHistSets; i++) {
    if (book_evt_histset[i] != 0) {
      art::TFileDirectory dir = run_dir.mkdir(std::format("evt_{}",i));
      fHist->evt[i] = new EvtHist_t;
      book_evt_histograms(fHist->evt[i],RunNumber,&dir);
    }
  }

//-----------------------------------------------------------------------------
// book segment histograms
//-----------------------------------------------------------------------------
  int book_seg_histset[kNSegHistSets];
  for (int i=0; i<kNSegHistSets; i++) book_seg_histset[i] = 0;

  book_seg_histset[0] = 1;              // all segments
  book_seg_histset[1] = 1;              // nhits > 10

  for (int i=0; i<kNSegHistSets; i++) {
    if (book_seg_histset[i] != 0) {
      art::TFileDirectory dir = run_dir.mkdir(std::format("seg_{}",i));
      fHist->seg[i] = new SegHist_t;
      book_seg_histograms(fHist->seg[i],RunNumber,&dir);
    }
  }
}

//-----------------------------------------------------------------------------
// art calls beginRun for each new run in the input
//-----------------------------------------------------------------------------
void ReadTrkSegments::beginRun(art::Run const& ArtRun) {
  int rn = ArtRun.run();
  if (rn == fRunNumber) return;
  fRunNumber = rn;
  book_histograms(rn);

  art::EventID eid(rn,1,1);
  fTracker = fAlignedTracker.getPtr(eid).get();

}

//-----------------------------------------------------------------------------
  void ReadTrkSegments::fill_eff_histograms(EffHist_t* Hist, SegmentHit* Hit, EffPar_t* EffP) {
    int panel = EffP->ip216;
    int straw = Hit->Straw();
    Hist->h_den[panel]->Fill(straw);
    if (Hit->ComboHit() != nullptr) {
      // there is a [some] hit in this straw
      Hist->h_num[panel]->Fill(straw);
    }

    float doca = Hit->Doca();
    Hist->h_r_den[panel]->Fill(doca);
    if (Hit->ComboHit() != nullptr) {
      Hist->h_r_num[panel]->Fill(doca);
    }
  }

//-----------------------------------------------------------------------------
void ReadTrkSegments::fill_evt_histograms(EvtHist_t* Hist, EvtPar_t* EvtP) {
  Hist->h_nseg->Fill(EvtP->nseg);
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::fill_seg_histograms(SegHist_t* Hist, TrkSegment const* Seg, SegPar_t* SegP) {
  Hist->h_nhits->Fill(Seg->NHitsWithComboHit());
  Hist->h_chi2d->Fill(Seg->Chi2Dof());
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::fill_histograms() {

  fEvtPar[0].nseg = fNTs;
  fill_evt_histograms(fHist->evt[0],&fEvtPar[0]);

  for (int is=0; is<fNTs; ++is) {
    TrkSegment const* seg = &fTsColl->at(is);
    SegPar_t* segp = &fSegPar[is];
    fill_seg_histograms(fHist->seg[0],seg,segp);
    int nh = seg->NHitsWithComboHit();
    if (nh >= 10) {
      fill_seg_histograms(fHist->seg[1],seg,segp);
    }
  }
//-----------------------------------------------------------------------------
// now - fill efficiency histograms
//-----------------------------------------------------------------------------
  EffPar_t effp;
  for (int is=0; is<fNTs; ++is) {
    TrkSegment const* seg = &clone_[is];
    if (seg->NHitsWithComboHit() < 10) continue;
    if (seg->Chi2Dof()           > 10) continue;

    int plane = seg->Plane();
    int panel = seg->Panel();
    effp.ip216 = 6*plane+panel;

    int nstraws = seg->NStraws();
    for (int is=0; is<nstraws; ++is) {
      SegmentHit* hit = seg->Hit(is);
      // a straw with a long-term or short-term status flag is filled with light gray
      bool good_straw = fTrackerStatus->strawStatus(StrawId(plane,panel,hit->Straw())).empty();
      if (not good_straw) continue;
      double doca = hit->Doca();
      if (doca < 2.5) {
        //expect to have a hit
        fill_eff_histograms(fHist->eff[0],hit,&effp);
      }
    }
  }
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::analyze(art::Event const& event) {

  fTrackerStatus  = fTrackerStatusH.getPtr(event.id());

  auto const segH = event.getValidHandle<TrkSegmentCollection>(fTsCollTag);
                                        // read the ComboHits explicitly, the segment hits refer to them
  art::InputTag const chCollTag = comboHitTagOf(segH,event);
  auto const chH  = event.getValidHandle<mu2e::ComboHitCollection>(chCollTag);
  std::cout << std::format("ReadTrkSegments: ComboHits {} nComboHits {}\n",chCollTag.encode(),chH->size());

  std::cout << std::format("ReadTrkSegments: event {} run {} subrun {} nsegments {}\n",
                           event.event(),event.run(),event.subRun(),segH->size());

  clone_.clear();
  for (size_t i=0; i<segH->size(); ++i) {
    TrkSegment seg = segH->at(i);
    seg.setComboHits(chH.product());
    clone_.push_back(seg);
  }

  fTsColl = segH.product();
  fNTs    = fTsColl->size();
  fChColl = chH.product();

  fEvtPar.resize(1);
  fSegPar.resize(fNTs);

  fill_histograms();
}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadTrkSegments)
