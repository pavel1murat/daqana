// Select events with at least the configured number of reconstructed tracks
// and CRV coincidence clusters.  Also write multiplicity histograms for all
// events and for events passing the filter.

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

namespace mu2e {

  class TrkCrvFilter : public art::EDFilter {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> ksCollTag  {Name("ksCollTag"  ), Comment("Tag of the KalSeedCollection")};
      fhicl::Atom<art::InputTag> crvcCollTag{Name("crvcCollTag"), Comment("Tag of the CrvCoincidenceClusterCollection")};
      fhicl::Atom<art::InputTag> tcCollTag  {Name("tcCollTag"  ), Comment("Tag of the TimeClusterCollection")};
      fhicl::Atom<unsigned>      minNCrvc   {Name("minNCrvc"   ), Comment("Minimum number of CRV coincidence clusters")};
      fhicl::Atom<unsigned>      minNTrkHits{Name("minNTrkHits"), Comment("Minimum number of track hits")};
      fhicl::Atom<unsigned>      minNTrk    {Name("minNTrk"    ), Comment("Minimum number of reconstructed tracks")};
      fhicl::Atom<unsigned>      minNTcHits {Name("minNTcHits" ), Comment("Minimum number of time cluster hits")};
      fhicl::Atom<float>         maxTcCrvcDt{Name("maxTcCrvcDt"), Comment("Max delta T between TC and CRVC")};
      fhicl::Atom<unsigned>      minNGoodTc {Name("minNGoodTc" ), Comment("Minimum number of time clusters")};
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit TrkCrvFilter(Parameters const& p);

    bool filter(art::Event& event) override;

  private:
    art::InputTag crvcCollTag_;
    art::InputTag ksCollTag_;
    art::InputTag tcCollTag_;
    
    unsigned      minNCrvc_;

    unsigned      minNTrkHits_;
    unsigned      minNTrk_    ;
    
    unsigned      minNTcHits_ ;
    float         maxTcCrvcDt_;
    unsigned      minNGoodTc_ ;

    TH1D* nTracksAll_         = nullptr;
    TH1D* nTracksPassed_      = nullptr;
    TH1D* nCrvClustersAll_    = nullptr;
    TH1D* nCrvClustersPassed_ = nullptr;
  };

  TrkCrvFilter::TrkCrvFilter(Parameters const& p)
    : art::EDFilter{p}
    , crvcCollTag_{p().crvcCollTag()}
    , ksCollTag_{p().ksCollTag()}
    , tcCollTag_{p().tcCollTag()}
    , minNCrvc_{p().minNCrvc()}
    , minNTrkHits_{p().minNTrkHits()}
    , minNTrk_   {p().minNTrk()}
    , minNTcHits_{p().minNTcHits()}
    , maxTcCrvcDt_{p().maxTcCrvcDt()}
    , minNGoodTc_{p().minNGoodTc()}
  {
    // consumes<KalSeedCollection>(ksCollTag_);
    consumes<CrvCoincidenceClusterCollection>(crvcCollTag_);
    consumes<TimeClusterCollection>(tcCollTag_);

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

  bool TrkCrvFilter::filter(art::Event& event) {
    //    auto ks_coll   = event.getValidHandle<KalSeedCollection>(ksCollTag_);
    auto crvc_coll = event.getValidHandle<CrvCoincidenceClusterCollection>(crvcCollTag_);
    auto tc_coll   = event.getValidHandle<TimeClusterCollection>(tcCollTag_);

    size_t nTrk = 0;
    // for (auto const& track : *ks_coll) {
    //   if (track.nHits() >= minNTrkHits_) ++nTrk;
    // }

    auto const nCrvc = crvc_coll->size();

    size_t n_good_tc = 0;
    for (auto const& tc : *tc_coll) {
      if (tc.nStrawHits() >= minNTcHits_) {
        for (auto const& crvc: *crvc_coll) {
          float tc_crvc_dt = tc.t0().t0()-crvc.GetStartTime();
          if (fabs(tc_crvc_dt) < maxTcCrvcDt_) {
            ++n_good_tc;
          }
        }
      }
    }

    bool const pass = (nTrk       >= minNTrk_   )
      &&              (nCrvc      >= minNCrvc_  )
      &&              (n_good_tc  >= minNGoodTc_);
    
    // histogramming
    nTracksAll_->Fill(static_cast<double>(nTrk));
    nCrvClustersAll_->Fill(static_cast<double>(nCrvc));

    if (pass) {
      nTracksPassed_->Fill(static_cast<double>(nTrk));
      nCrvClustersPassed_->Fill(static_cast<double>(nCrvc));
    }

    return pass;
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrkCrvFilter)
