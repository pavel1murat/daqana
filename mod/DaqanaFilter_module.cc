// Select events with at least the configured number of reconstructed tracks
// and CRV coincidence clusters.  Also write multiplicity histograms for all
// events and for events passing the filter.

// assume maxDt is large enough to be used for any selection w/o introducing a timing window per selection type

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"

#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Provenance/ProductID.h"

#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "TH1D.h"

#include "TRACE/trace.h"
#define TRACE_NAME "DaqanaFilter"


namespace mu2e {

  class DaqanaFilter : public art::EDFilter {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> ksCollTag    {Name("ksCollTag"  ), Comment("Tag of the KalSeedCollection")};
      fhicl::Atom<art::InputTag> calcCollTag  {Name("calcCollTag"), Comment("Tag of the CaloClusterClusterCollection")};
      fhicl::Atom<art::InputTag> crvcCollTag  {Name("crvcCollTag"), Comment("Tag of the CrvCoincidenceClusterCollection")};
      fhicl::Atom<art::InputTag> tcCollTag    {Name("tcCollTag"  ), Comment("Tag of the TimeClusterCollection")};
      fhicl::Atom<art::InputTag> chCollTag    {Name("chCollTag"  ), Comment("Tag of the ComboHitCollection")};

      fhicl::Atom<std::string>   selectionType{Name("selectionType"), Comment("selection type:tc_crvc/tc_calc/calc_only/tc_only")};

      // calorimeter cluster (CALC) selection
      fhicl::Atom<float>         minCalcEnergy{Name("minCalcEnergy"), Comment("Minimum CALC energy")};
      fhicl::Atom<int>           minCalcSize  {Name("minCalcSize"  ), Comment("Minimum CALC size (ncrystals)")};
      fhicl::Atom<int>           minNGoodCalc {Name("minNGoodCalc" ), Comment("Minimum number of calo clusters")};

      // muon stub (CRVC) selection - no real cuts 
      fhicl::Atom<unsigned>      minNGoodCrvc {Name("minNGoodCrvc" ), Comment("Minimum number of good CRVC")};

      // tracker time cluster (TC) selection
      fhicl::Atom<unsigned>      minNTcHits   {Name("minNTcHits"   ), Comment("Minimum number of time cluster hits")};
      fhicl::Atom<int>           minNTcPlanes {Name("minNTcPlanes" ), Comment("if > 0, minimum number of planes")};
      fhicl::Atom<int>           maxNTcPlanes {Name("maxNTcPlanes" ), Comment("if > 0, maximum number of planes")};
      fhicl::Atom<unsigned>      minNGoodTc   {Name("minNGoodTc"   ), Comment("Minimum number of time clusters")};
      
      // track selection (not used yet)
      fhicl::Atom<unsigned>      minNTrkHits  {Name("minNTrkHits"  ), Comment("Minimum number of track hits")};
      fhicl::Atom<unsigned>      minNTrk      {Name("minNTrk"      ), Comment("Minimum number of reconstructed tracks")};

      // timing window - common
      fhicl::Atom<float>         maxDt        {Name("maxDt"        ), Comment("Max allowed |DT(O1-O2)|")};
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit DaqanaFilter(Parameters const& p);

    // std::string findComboHitCollTag(art::Event const& event, art::Handle<TimeClusterCollection> const& TcHandle);
    bool        filter             (art::Event& event) override;
    
    bool filter_calc_only(art::Event& event);
    bool filter_tc_calc  (art::Event& event);
    bool filter_tc_crvc  (art::Event& event);
    bool filter_tc_only  (art::Event& event);
    bool filter_two_disks(art::Event& event);

    bool tc_is_good(const TimeCluster& Tc);
    
  private:
    art::InputTag calcCollTag_;
    art::InputTag chCollTag_;
    art::InputTag crvcCollTag_;
    art::InputTag ksCollTag_;
    art::InputTag tcCollTag_;

    const ComboHitCollection*                ch_coll_  ;
    const TimeClusterCollection*             tc_coll_  ;
    const CaloClusterCollection*             calc_coll_;
    const CrvCoincidenceClusterCollection*   crvc_coll_;

    std::string   selectionType_;
    
    float         minCalcEnergy_;
    int           minCalcSize_;
    int           minNGoodCalc_;

    unsigned      minNGoodCrvc_;

    unsigned      minNTrkHits_;
    unsigned      minNTrk_    ;

    unsigned      minNTcHits_ ;
    int           minNTcPlanes_;
    int           maxNTcPlanes_;
    unsigned      minNGoodTc_ ;

    
    float         maxDt_;
   
    TH1D* nTracksAll_         = nullptr;
    TH1D* nTracksPassed_      = nullptr;
    TH1D* nCrvClustersAll_    = nullptr;
    TH1D* nCrvClustersPassed_ = nullptr;
  };

//-----------------------------------------------------------------------------
  DaqanaFilter::DaqanaFilter(Parameters const& p)
    : art::EDFilter{p}
    , calcCollTag_{p().calcCollTag()}
    , chCollTag_{p().chCollTag()}
    , crvcCollTag_{p().crvcCollTag()}
    , ksCollTag_{p().ksCollTag()}
    , tcCollTag_{p().tcCollTag()}
      
    , selectionType_{p().selectionType()}

    , minCalcEnergy_{p().minCalcEnergy()}
    , minCalcSize_  {p().minCalcSize  ()}
    , minNGoodCalc_ {p().minNGoodCalc ()}
      
    , minNGoodCrvc_{p().minNGoodCrvc()}
      
    , minNTrkHits_{p().minNTrkHits()}
    , minNTrk_   {p().minNTrk()}
      
    , minNTcHits_  {p().minNTcHits()  }
    , minNTcPlanes_{p().minNTcPlanes()}
    , maxNTcPlanes_{p().maxNTcPlanes()}
    , minNGoodTc_  {p().minNGoodTc()  }

    , maxDt_{p().maxDt()}

  {
    // consumes<KalSeedCollection>(ksCollTag_);
    consumes<CaloClusterCollection>          (calcCollTag_);
    consumes<CrvCoincidenceClusterCollection>(crvcCollTag_);
    consumes<TimeClusterCollection>          (tcCollTag_  );

    art::ServiceHandle<art::TFileService> tfs;
    auto all    = tfs->mkdir("all");
    auto passed = tfs->mkdir("passed");

    // The overflow bins retain events whose multiplicity exceeds the plotted
    // range.  The range can be enlarged if the expected occupancy requires it.
    nTracksAll_         = all.make<TH1D>("nTracks", "Number of reconstructed tracks;N tracks;Events", 100, -0.5, 99.5);
    nCrvClustersAll_    = all.make<TH1D>("nCrvClusters", "Number of CRV coincidence clusters;N clusters;Events",100, -0.5, 99.5);
    nTracksPassed_      = passed.make<TH1D>("nTracks", "Number of reconstructed tracks;N tracks;Events", 100, -0.5, 99.5);
    nCrvClustersPassed_ = passed.make<TH1D>("nCrvClusters", "Number of CRV coincidence clusters;N clusters;Events",100, -0.5, 99.5);
  }

//-----------------------------------------------------------------------------
// generated by Claude
// tc_handle: handle of the TimeClusterCollection produced by the time-cluster module
// returns the InputTag string "label:instance:process" of the ComboHitCollection used to build it,
// or an empty string if no ComboHitCollection is among the parents
//-----------------------------------------------------------------------------
  // std::string DaqanaFilter::findComboHitCollTag(art::Event const& event, art::Handle<TimeClusterCollection> const& tc_handle) {
  //   std::string tag;

  //   art::Provenance const* prov = tc_handle.provenance();
  //   for (art::ProductID const& pid : prov->parents()) {
  //     cet::exempt_ptr<art::BranchDescription const> desc = event.getProductDescription(pid);
  //     if (desc.get() == nullptr) continue;

  //     if (desc->producedClassName() == "mu2e::ComboHitCollection") {
  //       tag = desc->moduleLabel() + ":" + desc->productInstanceName() + ":" + desc->processName();
  //       break;
  //     }
  //   }
    
  //   return tag;
  // }
//-----------------------------------------------------------------------------
// require at least one cluster above minCalcEnergy_ in each disk
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter_calc_only(art::Event& event) {
    bool passed  = false;

    calc_coll_ = event.getValidHandle<CaloClusterCollection>(calcCollTag_);
    int nCalc = calc_coll_->size();

    int n_good_cl(0);
    for (int i=0; i<nCalc; i++) {
      const CaloCluster* c = &calc_coll_->at(i);
      if (c->energyDep() < minCalcEnergy_) continue;
      if (c->size     () < minCalcSize_  ) continue;
      n_good_cl++;
    }
    passed = (n_good_cl > minNGoodCalc_);
    return passed;
  }

//-----------------------------------------------------------------------------
// require at least one cluster above minCalcEnergy_ in each disk
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter_two_disks(art::Event& event) {
    bool passed  = false;

    calc_coll_ = event.getValidHandle<CaloClusterCollection>(calcCollTag_);
    int nCalc = calc_coll_->size();

    int nmatches(0);
    for (int i1=0; i1<nCalc; i1++) {
      const CaloCluster* c1 = &calc_coll_->at(i1);
      if (c1->energyDep() < minCalcEnergy_) continue;
      if (c1->size     () < minCalcSize_  ) continue;

      for (int i2=i1+1; i2<nCalc; i2++) {
        const CaloCluster* c2 = &calc_coll_->at(i2);
        if (c2->energyDep() < minCalcEnergy_) continue;
        if (c2->size     () < minCalcSize_  ) continue;
            
        // clusters in different disks
        if (c1->diskID() == c2->diskID()) continue;
            
        float dt = c1->time()-c2->time();
        if (fabs(dt) < maxDt_) {
          nmatches++;
        }
      }
    }
    passed = (nmatches > 0);
    return passed;
  }

//-----------------------------------------------------------------------------
  bool DaqanaFilter::tc_is_good(const TimeCluster& Tc) {
    bool good{false};

    int nh_plane[36]{};                 // all elements initialized to zero
    int n_tc_planes{0};

    unsigned tc_nhits = Tc.nhits();
    if (tc_nhits < minNTcHits_)    return good;
    // count number of planes

    for (size_t ih=0; ih<tc_nhits; ih++) {
      StrawHitIndex hit_index   = Tc.hits().at(ih);
      const mu2e::ComboHit* hit = &ch_coll_->at(hit_index);
      int plane = hit->strawId().plane();
      nh_plane[plane]++;
      if (nh_plane[plane] == 1) {
        n_tc_planes++;
      }
    }
    
    if ((minNTcPlanes_ > 0) and (n_tc_planes < minNTcPlanes_)) return good;
    if ((maxNTcPlanes_ > 0) and (n_tc_planes > maxNTcPlanes_)) return good;

    good = true;
    return good;
  }

  
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter_tc_calc(art::Event& event) {
    bool passed  = false;
    
    calc_coll_ = event.getValidHandle<CaloClusterCollection>(calcCollTag_);
    ch_coll_   = event.getValidHandle<ComboHitCollection>(chCollTag_);
    tc_coll_   = event.getValidHandle<TimeClusterCollection>(tcCollTag_);


    size_t n_good = 0;
    for (auto const& tc : *tc_coll_) {
      // these are combo hits
      if (not tc_is_good(tc)) continue;

      int n_good_calc = 0;
      for (auto const& calc: *calc_coll_) {
        if (calc.energyDep() > minCalcEnergy_) continue;
        if (calc.size()      > minCalcSize_  ) continue;

        float dt = tc.t0().t0()-calc.time();
        if (fabs(dt) < maxDt_) {
          // "good" combination found
          n_good_calc++;
        }
      }

      if (n_good_calc >= minNGoodCalc_) {
        n_good++;
      }
    }

    passed = (n_good >= minNGoodTc_); 
    
    return passed;
  }
    
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter_tc_crvc(art::Event& event) {
    bool passed  = false;
    
    ch_coll_   = event.getValidHandle<ComboHitCollection>(chCollTag_);
    crvc_coll_ = event.getValidHandle<CrvCoincidenceClusterCollection>(crvcCollTag_);
    tc_coll_   = event.getValidHandle<TimeClusterCollection>(tcCollTag_);

    size_t n_good   = 0;
    for (auto const& tc : *tc_coll_) {
      // these are the combo hits
      if (not tc_is_good(tc)) continue;

      size_t n_good_crvc = 0;

      for (auto const& crvc: *crvc_coll_) {
        float tc_crvc_dt = tc.t0().t0()-crvc.GetStartTime();
        if (fabs(tc_crvc_dt) < maxDt_) {
          ++n_good_crvc;
        }
      }
      // number of CRV clusters within the window
      if (n_good_crvc >= minNGoodCrvc_) {
        n_good++;
      }
    }

    passed = (n_good >= minNGoodTc_);
    
    return passed;
  }
    
//-----------------------------------------------------------------------------
// no need to get ChCollTag from FCL
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter_tc_only(art::Event& event) {
    bool passed  = false;
    
    ch_coll_   = event.getValidHandle<ComboHitCollection>(chCollTag_);
    tc_coll_   = event.getValidHandle<TimeClusterCollection>(tcCollTag_);

    size_t n_good_tc   = 0;
    for (auto const& tc : *tc_coll_) {
      if (not tc_is_good(tc)) continue;
      n_good_tc++;
    }
    
    passed = (n_good_tc >= minNGoodTc_);
    
    return passed;
  }
  
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter(art::Event& event) {
    //    auto ks_coll   = event.getValidHandle<KalSeedCollection>(ksCollTag_);

    //    size_t nTrk = 0;

    bool passed(false);
    
    if (selectionType_ == "calc_only") {
      passed = filter_calc_only(event);
    }
    else if (selectionType_ == "tc_calc") {
      passed = filter_tc_calc(event);
    }
    else if (selectionType_ == "tc_crvc") {
      passed = filter_tc_crvc(event);
    }
    else if (selectionType_ == "tc_only") {
      passed = filter_tc_only(event);
    }
    else if (selectionType_ == "two_disks") {
      passed = filter_two_disks(event);
    }
    else {
      TLOG(TLVL_ERROR) << std::format("unknown selection type:{}",selectionType_);
    }

    // // histogramming
    // nTracksAll_->Fill(static_cast<double>(nTrk));
    // nCrvClustersAll_->Fill(static_cast<double>(nCrvc));

    // if (passed) {
    //   nTracksPassed_->Fill(static_cast<double>(nTrk));
    //   nCrvClustersPassed_->Fill(static_cast<double>(nCrvc));
    // }

    return passed;
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DaqanaFilter)
