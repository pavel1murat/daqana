// Select events with at least the configured number of reconstructed tracks
// and CRV coincidence clusters.  Also write multiplicity histograms for all
// events and for events passing the filter.

// assume maxDt is large enough to be used for any selection w/o introducing a timing window per selection type

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
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

      fhicl::Atom<std::string>   selectionType{Name("selectionType"), Comment("selection type:tc_crvc/tc_calc/calc_only/tc_only")};

      // calorimeter cluster (CALC) selection
      fhicl::Atom<float>         minCalcEnergy{Name("minCalcEnergy"), Comment("Minimum CALC energy")};
      fhicl::Atom<unsigned>      minNGoodCalc {Name("minNGoodCalc" ), Comment("Minimum number of calo clusters")};

      // muon stub (CRVC) selection - no real cuts 
      fhicl::Atom<unsigned>      minNGoodCrvc {Name("minNGoodCrvc" ), Comment("Minimum number of good CRVC")};

      // tracker time cluster (TC) selection
      fhicl::Atom<unsigned>      minNTcHits   {Name("minNTcHits"   ), Comment("Minimum number of time cluster hits")};
      fhicl::Atom<unsigned>      minNGoodTc   {Name("minNGoodTc"   ), Comment("Minimum number of time clusters")};
      
      // track selection (not used yet)
      fhicl::Atom<unsigned>      minNTrkHits  {Name("minNTrkHits"  ), Comment("Minimum number of track hits")};
      fhicl::Atom<unsigned>      minNTrk      {Name("minNTrk"      ), Comment("Minimum number of reconstructed tracks")};

      // timing window - common
      fhicl::Atom<float>         maxDt        {Name("maxDt"        ), Comment("Max allowed |DT(O1-O2)|")};
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit DaqanaFilter(Parameters const& p);

    bool filter          (art::Event& event) override;
    
    bool select_calc_only(art::Event& event);
    bool select_tc_calc  (art::Event& event);
    bool select_tc_crvc  (art::Event& event);
    bool select_tc_only  (art::Event& event);
    
  private:
    art::InputTag calcCollTag_;
    art::InputTag crvcCollTag_;
    art::InputTag ksCollTag_;
    art::InputTag tcCollTag_;

    const TimeClusterCollection*            tc_coll_  ;
    const CaloClusterCollection*            calc_coll_;
    const CrvCoincidenceClusterCollection*  crvc_coll_;

    std::string   selectionType_;
    
    float         minCalcEnergy_;
    unsigned      minNGoodCalc_;

    unsigned      minNGoodCrvc_;

    unsigned      minNTrkHits_;
    unsigned      minNTrk_    ;

    unsigned      minNTcHits_ ;
    unsigned      minNGoodTc_ ;

    
    float         maxDt_;
   
    TH1D* nTracksAll_         = nullptr;
    TH1D* nTracksPassed_      = nullptr;
    TH1D* nCrvClustersAll_    = nullptr;
    TH1D* nCrvClustersPassed_ = nullptr;
  };

  DaqanaFilter::DaqanaFilter(Parameters const& p)
    : art::EDFilter{p}
    , calcCollTag_{p().calcCollTag()}
    , crvcCollTag_{p().crvcCollTag()}
    , ksCollTag_{p().ksCollTag()}
    , tcCollTag_{p().tcCollTag()}
      
    , selectionType_{p().selectionType()}

    , minCalcEnergy_{p().minCalcEnergy()}
    , minNGoodCalc_ {p().minNGoodCalc() }
      
    , minNGoodCrvc_{p().minNGoodCrvc()}
      
    , minNTrkHits_{p().minNTrkHits()}
    , minNTrk_   {p().minNTrk()}
      
    , minNTcHits_{p().minNTcHits()}
    , minNGoodTc_{p().minNGoodTc()}

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
// require at least one cluster above minCalcEnergy_ in each disk
//-----------------------------------------------------------------------------
  bool DaqanaFilter::select_calc_only(art::Event& event) {
    bool passed  = false;

    int nCalc = calc_coll_->size();

    int nmatches(0);
    for (int i1=0; i1<nCalc; i1++) {
      const CaloCluster* c1 = &calc_coll_->at(i1);
      if (c1->energyDep() < minCalcEnergy_) continue;

      for (int i2=i1+1; i2<nCalc; i2++) {
        const CaloCluster* c2 = &calc_coll_->at(i2);
        if (c2->energyDep() < minCalcEnergy_) continue;
            
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
  bool DaqanaFilter::select_tc_calc(art::Event& event) {
    bool passed  = false;
    
    size_t n_good_tc   = 0;
    size_t n_good_calc = 0;
    for (auto const& tc : *tc_coll_) {
      if (tc.nhits() < minNTcHits_)     continue;
      n_good_tc++;
        
      for (auto const& calc: *calc_coll_) {
        if (calc.energyDep() < minCalcEnergy_) continue;
          
        float dt = tc.t0().t0()-calc.time();
        if (fabs(dt) < maxDt_) {
          ++n_good_calc;
        }
      }
    }
    passed = (n_good_tc >= minNGoodTc_) and (n_good_calc >= minNGoodCalc_);
    
    return passed;
  }
  
//-----------------------------------------------------------------------------
  bool DaqanaFilter::select_tc_crvc(art::Event& event) {
    bool passed  = false;
    
    size_t n_good_tc   = 0;
    size_t n_good_crvc = 0;
    for (auto const& tc : *tc_coll_) {
      // these are combo hits
      if (tc.nhits() < minNTcHits_) continue;
      n_good_tc++;
      for (auto const& crvc: *crvc_coll_) {
        float tc_crvc_dt = tc.t0().t0()-crvc.GetStartTime();
        if (fabs(tc_crvc_dt) < maxDt_) {
          ++n_good_crvc;
        }
      }
    }

    passed = (n_good_tc>= minNGoodTc_) and (n_good_crvc >= minNGoodCrvc_);
    
    return passed;
  }
    
//-----------------------------------------------------------------------------
  bool DaqanaFilter::select_tc_only(art::Event& event) {
    bool passed  = false;
    
    size_t n_good_tc   = 0;
    for (auto const& tc : *tc_coll_) {
      if (tc.nStrawHits() < minNTcHits_)     continue;
      n_good_tc++;
    }
    
    passed = (n_good_tc >= minNGoodTc_);
    
    return passed;
  }
  
//-----------------------------------------------------------------------------
  bool DaqanaFilter::filter(art::Event& event) {
    //    auto ks_coll   = event.getValidHandle<KalSeedCollection>(ksCollTag_);
    crvc_coll_ = event.getValidHandle<CrvCoincidenceClusterCollection>(crvcCollTag_);
    calc_coll_ = event.getValidHandle<CaloClusterCollection>(calcCollTag_);
    tc_coll_   = event.getValidHandle<TimeClusterCollection>(tcCollTag_);

    //    size_t nTrk = 0;

    bool passed(false);
    
    if (selectionType_ == "calc_only") {
      passed = select_calc_only(event);
    }
    else if (selectionType_ == "tc_calc") {
      passed = select_tc_calc(event);
    }
    else if (selectionType_ == "tc_crvc") {
      passed = select_tc_calc(event);
    }
    else if (selectionType_ == "calc_only") {
      passed = select_calc_only(event);
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
