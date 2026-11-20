///////////////////////////////////////////////////////////////////////////////
// ReadTrkSegments: test module, reads a collection of TrkSegments from the input file
//                  and prints the segments with their hits; the hits are resolved
//                  through their index in the ComboHitCollection of the file.
//
// histograms: booked in beginRun, one directory per run, run_<N>, with a set of
//             histograms per event (evt_<i>) and per segment (seg_<i>) for each set i.
//             The sets are switched on with the booking flags in book_histograms;
//             they are filled in analyze.
// debugBit_[3] : print segments
// debugBit_[4] : debug, per fill_pnl_histograms entry
// debugBit_[5] : debug, print outliers indr vs doca
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
#include "TH2.h"

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
    kNPnlHistSets = 216,
    kNSegHistSets = 100,
  };

//-----------------------------------------------------------------------------
// FCL parameters
//-----------------------------------------------------------------------------
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag>   tsCollTag{Name("tsCollTag"), Comment("the collection of TrkSegments to read")};
    fhicl::Atom<int>             minPanel {Name("minPanel" ), Comment("min panel number (0-215)")};
    fhicl::Atom<int>             maxPanel {Name("maxPanel" ), Comment("max panel number (0-215)")};
    fhicl::Atom<bool>            debugMode{Name("debugMode"), Comment("ON/OFF switch for all debug printouts")};
    fhicl::Sequence<std::string> debugBits{Name("debugBits"), Comment("Debug switches in bitN:value form")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  art::InputTag fTsCollTag;
  int           fMinPanel;
  int           fMaxPanel;
  bool          fDebugMode;
  int           debugBit_[100]{};
//-----------------------------------------------------------------------------
// parameters of an event and of a segment, to be used by the histograms
//-----------------------------------------------------------------------------
  struct EffPar_t {
    int  ip216;                         // offline panel number in the range [0,215]
  };

  struct EvtPar_t {
    int  nseg;
  };

  struct HitPar_t {
    int   hitp;
    float seg_dydx;                       // segment dydx, need for alignment
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
    TH1F* h_nseg;                         // number of segments per event
  };

  struct PnlHist_t {
    TH1F* h_dr [96][2];                   // left and right residuals
    TH1F* h_dr2[96][2];                   // up and down residuals
    TH2F* h_dr_vs_doca;                   // drift vel
  };

  struct SegHist_t {
    TH1F* h_nhits_all;                     // number of hits with a ComboHit per segment
    TH1F* h_chi2d_all;                     // chi2/dof per segment
    TH1F* h_nhits[216];                    // number of hits with a ComboHit per segment
    TH1F* h_chi2d[216];                    // chi2/dof per segment
  };

  struct Hist_t {
    EffHist_t* eff[kNEffHistSets];
    EvtHist_t* evt[kNEvtHistSets];
    PnlHist_t* pnl[216];
    SegHist_t* seg[kNSegHistSets];
  };
//-----------------------------------------------------------------------------
// other parameters
//-----------------------------------------------------------------------------
  int                         fRunNumber;
  Hist_t*                     fHist;
  bool                        fHistogramsBooked;

  int                         fNTs;     // number of track segments
  const ComboHitCollection*   fChColl;  // the ComboHits of the current event

  std::vector<EvtPar_t>       fEvtPar;  // event parameters, for the event histograms
  std::vector<SegPar_t>       fSegPar;  // segment parameters, one per segment of the current event

  std::vector<TrkSegment>     fListOfSegments;
  const art::Event*           fEvent;
//-----------------------------------------------------------------------------
// other parameters
//-----------------------------------------------------------------------------
  ProditionsHandle<Tracker>            fAlignedTracker;
  ProditionsHandle<TrackerStatus>      fTrackerStatusH;
  std::shared_ptr<const TrackerStatus> fTrackerStatus;
  const Tracker*                       fTracker;


  //  int                                  debugBit_[100]{};

  explicit ReadTrkSegments(Parameters const& config);

  void beginRun(art::Run const& ArtRun) override;
  void analyze(art::Event const& event) override;

  art::InputTag comboHitTagOf(art::ValidHandle<TrkSegmentCollection> const& SegH, art::Event const& Event) const;

  void book_eff_histograms(EffHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_evt_histograms(EvtHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_pnl_histograms(PnlHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_seg_histograms(SegHist_t* Hist, int RunNumber, art::TFileDirectory* Dir);
  void book_histograms(int RunNumber);

  void fill_eff_histograms(EffHist_t* Hist, SegmentHit* Hit, EffPar_t* EffP);
  void fill_evt_histograms(EvtHist_t* Hist, EvtPar_t* EvtP);
  void fill_pnl_histograms(PnlHist_t* Hist, const SegmentHit* Hit, HitPar_t* HitP);
  void fill_seg_histograms(SegHist_t* Hist, const TrkSegment* Seg, SegPar_t* SegP);
  void fill_histograms();

  void print(int DebugBit);
};

//-----------------------------------------------------------------------------
ReadTrkSegments::ReadTrkSegments(Parameters const& config) :
  art::EDAnalyzer(config),
  fRunNumber(-1),
  fHist(nullptr),
  fNTs(0),
  fChColl(nullptr) {
  
  fTsCollTag = config().tsCollTag();
  fMinPanel  = config().minPanel();
  fMaxPanel  = config().maxPanel();
  fDebugMode = config().debugMode();

  fHistogramsBooked = false;
  
  for (auto const& setting : config().debugBits()) {
    int index = -1;
    int value = 0;
    int consumed = 0;
    if (std::sscanf(setting.c_str(), "bit_%02d:%d%n", &index, &value, &consumed) != 2 ||
        setting[consumed] != '\0' || index < 0 || index >= 100) {
      std::cout << "ReadTrkFragments" << "Ignoring invalid debug bit setting '" << setting
                                       << "' (expected bitN:value with 0 <= N < 100)";
      continue;
    }
    debugBit_[index] = value;
  }
  
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
void ReadTrkSegments::book_pnl_histograms(PnlHist_t* Hist, int RunNumber, art::TFileDirectory* Dir) {
  std::string name;
  std::string title;
  for (int is=0; is<96; is++) {
    for (int i=0; i<2; i++) {
      name  = std::format("dr_{:02d}_{}",is,i);
      title = std::format("dr straw:{:02d} dright_sign:{}",is,i);
      Hist->h_dr[is][i] = Dir->make<TH1F>(name.data(),title.data(),200,-5,5);
      
      name  = std::format("dr2_{:02d}_{}",is,i);
      if   (i == 0) title = std::format("dr2 straw:{:02d} seg_dydx > 0",is);
      else          title = std::format("dr2 straw:{:02d} seg_dydx < 0",is);
      Hist->h_dr2[is][i] = Dir->make<TH1F>(name.data(),title.data(),200,-5,5);
    }
  }
  
  name   = std::format("dr_vs_doca");
  title  = std::format("dr vs doca");
  Hist->h_dr_vs_doca = Dir->make<TH2F>(name.data(),title.data(),300,0,3,200,-5,5);
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::book_seg_histograms(SegHist_t* Hist, int RunNumber, art::TFileDirectory* Dir) {
  std::string name, title;
  
  Hist->h_nhits_all = Dir->make<TH1F>("nhits_all","number of hits with a ComboHit per segment, all panels",60,0,60);
  Hist->h_chi2d_all = Dir->make<TH1F>("chi2d_all","chi2/dof per segment, all panels",100,0,50);

  for (int i=0; i<216; i++) {
    name  = std::format("nhits_{:03d}",i);
    title = std::format("nhits panel:{:03d}",i);
    Hist->h_nhits[i] = Dir->make<TH1F>(name.data(),title.data(),60,0,60);
    name  = std::format("chi2d_{:03d}",i);
    title = std::format("segment chi2d panel:{:03d}",i);
    Hist->h_chi2d[i] = Dir->make<TH1F>(name.data(),title.data(),100,0,50);
  }
}

//-----------------------------------------------------------------------------
// the sets switched on: event set 0 and segment set 0; the other sets are booked by
// setting their flag to 1, and then filled in fill_histograms
//-----------------------------------------------------------------------------
void ReadTrkSegments::book_histograms(int RunNumber) {

  art::ServiceHandle<art::TFileService> tfs;
  TH1::AddDirectory(kFALSE);

  //  art::TFileDirectory run_dir = tfs->mkdir(std::format("run_{}",RunNumber));

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
      //      art::TFileDirectory dir = run_dir.mkdir(std::format("eff_{}",i));
      art::TFileDirectory dir = tfs->mkdir(std::format("eff_{}",i));
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
      art::TFileDirectory dir = tfs->mkdir(std::format("evt_{}",i));
      fHist->evt[i] = new EvtHist_t;
      book_evt_histograms(fHist->evt[i],RunNumber,&dir);
    }
  }

//-----------------------------------------------------------------------------
// book panel histograms
//-----------------------------------------------------------------------------
  int book_pnl_histset[kNPnlHistSets];
  for (int i=0; i<kNPnlHistSets; i++) book_pnl_histset[i] = 0;

  for (int i=fMinPanel; i<fMaxPanel+1; i++) {
    book_pnl_histset[i] = 1;              // all events
  }

  for (int i=0; i<kNPnlHistSets; i++) {
    if (book_pnl_histset[i] != 0) {
      art::TFileDirectory dir = tfs->mkdir(std::format("pnl_{}",i));
      fHist->pnl[i] = new PnlHist_t;
      book_pnl_histograms(fHist->pnl[i],RunNumber,&dir);
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
      art::TFileDirectory dir = tfs->mkdir(std::format("seg_{}",i));
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
  
  // book histograms only once
  if (not fHistogramsBooked) {
    book_histograms(rn);
    fHistogramsBooked = true;
  }

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
  void ReadTrkSegments::fill_pnl_histograms(PnlHist_t* Hist, const SegmentHit* Hit, HitPar_t* HitP) {
    int is = Hit->Straw();
    float doca = fabs(Hit->Doca());
    float dr = Hit->R()-doca;
    if (fDebugMode and debugBit_[4]) {
      std::cout << std::format("Hit plane:panel:straw:{:02d}:{:d}:{:02d} R:{:10.3f} doca:{:10.3f} dr:{:10.3f}\n",
                               Hit->Sid().plane(),Hit->Sid().panel(),Hit->Sid().straw(),
                               Hit->R(),Hit->Doca(),dr);
    }

    Hist->h_dr_vs_doca->Fill(doca,dr);

    if (doca > 0.5) {
      // plot residuals only for hits with well defined radius and drift sign
      if (Hit->Drs() > 0) {
        Hist->h_dr[is][0]->Fill(dr);
      }
      else if (Hit->Drs() < 0) {
        Hist->h_dr[is][1]->Fill(dr);
      }

      if (HitP->seg_dydx > 0) {
        Hist->h_dr2[is][0]->Fill(dr);
      }
      else {
        Hist->h_dr2[is][1]->Fill(dr);
      }
    }
    
    if (fDebugMode and debugBit_[5]) {
      if ((doca > 2.5) and (dr < 0.6)) {
        std::cout << std::format("ReadTrkSegments bit_05 : {:6d}:{:6d}:{:8d} doca:{:10.3f} dr:{:10.3f}\n",
                                 fEvent->run(),fEvent->subRun(),fEvent->event(),doca,dr);
      }
    }
  }

//-----------------------------------------------------------------------------
  void ReadTrkSegments::fill_seg_histograms(SegHist_t* Hist, const TrkSegment* Seg, SegPar_t* SegP) {
    Hist->h_nhits_all->Fill(Seg->NHitsWithComboHit());
    Hist->h_chi2d_all->Fill(Seg->Chi2Dof());
    
    int p216 = Seg->Plane()*6+Seg->Panel();
    
    Hist->h_nhits[p216]->Fill(Seg->NHitsWithComboHit());
    Hist->h_chi2d[p216]->Fill(Seg->Chi2Dof());
  }

//-----------------------------------------------------------------------------
void ReadTrkSegments::fill_histograms() {

  fEvtPar[0].nseg = fNTs;
  fill_evt_histograms(fHist->evt[0],&fEvtPar[0]);

  for (int is=0; is<fNTs; ++is) {
    TrkSegment const* seg = &fListOfSegments.at(is);
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
  HitPar_t hitp;
  for (int is=0; is<fNTs; ++is) {
    TrkSegment const* seg = &fListOfSegments[is];
    if (seg->NHitsWithComboHit() < 10) continue;
    if (seg->Chi2Dof()           > 10) continue;
    // skip vertical segments
    if (fabs(seg->Ny()) < 0.01) continue;

    int plane = seg->Plane();
    int panel = seg->Panel();
    int p216  = 6*plane+panel;

    effp.ip216 = p216;
//-----------------------------------------------------------------------------
// the first and the last hits are always present, therefore they should not be used
// for the efficiency determination
//-----------------------------------------------------------------------------
    hitp.seg_dydx = seg->DyDx(); 
    int nstraws = seg->NStraws();
    for (int is=1; is<nstraws-1; ++is) {
      SegmentHit* hit = seg->Hit(is);
      // a straw with a long-term or short-term status flag is filled with light gray
      bool good_straw = fTrackerStatus->strawStatus(StrawId(plane,panel,hit->Straw())).empty();
      
      if (not good_straw) continue;
      double doca = hit->Doca();
      if (doca < 2.5) {
        //expect to have a hit
        fill_eff_histograms(fHist->eff[0],hit,&effp);
      }

      // these histograms have to be filled for straws with hits
      if (hit->ComboHit() == nullptr) continue;
   
      if ((p216 >= fMinPanel) and (p216 <= fMaxPanel)) {
        fill_pnl_histograms(fHist->pnl[p216],hit,&hitp);
      }
    }
  }
}

//-----------------------------------------------------------------------------
void ReadTrkSegments::analyze(art::Event const& event) {

  fTrackerStatus  = fTrackerStatusH.getPtr(event.id());

  auto const segH = event.getValidHandle<TrkSegmentCollection>(fTsCollTag);

  fEvent = &event;
  
                                        // read the ComboHits explicitly, the segment hits refer to them
  art::InputTag const chCollTag = comboHitTagOf(segH,event);
  auto          const chH       = event.getValidHandle<mu2e::ComboHitCollection>(chCollTag);
  
  // std::cout << std::format("ReadTrkSegments: ComboHits {} nComboHits {}\n",chCollTag.encode(),chH->size());

  // std::cout << std::format("ReadTrkSegments: event {} run {} subrun {} nsegments {}\n",
  //                          event.event(),event.run(),event.subRun(),segH->size());

  fListOfSegments.clear();
  for (size_t i=0; i<segH->size(); ++i) {
    TrkSegment seg = segH->at(i);
    seg.setComboHits(chH.product());
    fListOfSegments.push_back(seg);
  }

  //  fTsColl = segH.product();

  if (fDebugMode and (debugBit_[3] == 1)) {
    print(3);
  }
  
  fNTs    = segH->size();
  fChColl = chH.product();

  fEvtPar.resize(1);
  fSegPar.resize(fNTs);

  fill_histograms();
}

//-----------------------------------------------------------------------------
  void ReadTrkSegments::print(int DebugBit) {

    if (not fDebugMode) return;
    
    std::cout << std::format("---- ReadTrkFragments:: bit:{}\n",DebugBit);

    if (DebugBit == 3) {
      int nseg = fListOfSegments.size();
      for (int i = 0; i < nseg; ++i) {
        auto* seg = &fListOfSegments.at(i);
        seg->print("on read");
      }
    }
    else {
      std::cout << std::format("bit {}: printout undefined\n",DebugBit);
    }

  }
}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadTrkSegments)
