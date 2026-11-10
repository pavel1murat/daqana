// Reconstruct per-panel tracker segments from ComboHits assigned to time clusters.
// The segments are event-local in this first version and are not persisted.

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "daqana/obj/SegmentFit.hh"
#include "daqana/obj/TrkSegment.hh"

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace mu2e {

class PSegFinder : public art::EDProducer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> comboHitTag{fhicl::Name("comboHitTag"), fhicl::Comment("ComboHit collection")};
    fhicl::Atom<art::InputTag> timeClusterTag{fhicl::Name("timeClusterTag"), fhicl::Comment("TimeCluster collection")};
    fhicl::Atom<double> driftVelocity{fhicl::Name("driftVelocity"), fhicl::Comment("Straw drift velocity in mm/ns")};
    fhicl::Atom<double> timeOffset{fhicl::Name("timeOffset"), fhicl::Comment("Drift-time offset in ns")};
    fhicl::Atom<int> minimumHits{fhicl::Name("minimumHits"), fhicl::Comment("Minimum hits required to fit a segment")};
    fhicl::Atom<int> maximumStrawGap{fhicl::Name("maximumStrawGap"), fhicl::Comment("Largest gap between hits in one subsegment")};
    fhicl::Atom<bool> verbose{fhicl::Name("verbose"), fhicl::Comment("Print the number of reconstructed panel segments")};
    fhicl::Sequence<std::string> debugBits{fhicl::Name("debugBits"), fhicl::Comment("Debug switches in bitN:value form")};
  };
  using Parameters = art::EDProducer::Table<Config>;


private:

  art::InputTag             comboHitTag_;
  art::InputTag             timeClusterTag_;
  int                       minimumHits_;
  int                       maximumStrawGap_;
  bool                      verbose_;
  ProditionsHandle<Tracker> alignedTracker_;
  int                       debugBit_[100]{};
public:

//-----------------------------------------------------------------------------
  explicit PSegFinder(Parameters const& config)
    : art::EDProducer(config)
    , comboHitTag_(config().comboHitTag())
    , timeClusterTag_(config().timeClusterTag())
    , minimumHits_(config().minimumHits())
    , maximumStrawGap_(config().maximumStrawGap())
    , verbose_(config().verbose())
  {
    for (auto const& setting : config().debugBits()) {
      int index = -1;
      int value = 0;
      int consumed = 0;
      if (std::sscanf(setting.c_str(), "bit%d:%d%n", &index, &value, &consumed) != 2 ||
          setting[consumed] != '\0' || index < 0 || index >= 100) {
        mf::LogWarning("PSegFinder") << "Ignoring invalid debug bit setting '" << setting
                                      << "' (expected bitN:value with 0 <= N < 100)";
        continue;
      }
      debugBit_[index] = value;
    }

    consumes<ComboHitCollection>(comboHitTag_);
    consumes<TimeClusterCollection>(timeClusterTag_);
    SegmentHit::SetVDrift(config().driftVelocity());
    SegmentHit::SetTOffset(config().timeOffset());
  }

//-----------------------------------------------------------------------------
  void produce(art::Event& event) override {
    auto const& comboHits = *event.getValidHandle<ComboHitCollection>(comboHitTag_);
    auto const& timeClusters = *event.getValidHandle<TimeClusterCollection>(timeClusterTag_);
    auto tracker = alignedTracker_.getPtr(event.id());

    // Keep one working segment per panel, as in the starting algorithm.
    std::vector<TrkSegment> panelSegments;
    panelSegments.reserve(36 * 6);
    for (int plane = 0; plane < 36; ++plane) {
      for (int panel = 0; panel < 6; ++panel) {
        panelSegments.emplace_back(plane, panel);
        auto const& geomPanel = tracker->getPlane(plane).getPanel(panel);
        panelSegments.back().fTrkPanel = const_cast<Panel*>(&geomPanel);
      }
    }

    // print the input time clusters with all their ComboHits
    if (debugBit_[9] > 0) {
      for (size_t icl = 0; icl < timeClusters.size(); ++icl) {
        auto const& timeCluster = timeClusters[icl];
        mf::LogVerbatim("PSegFinder") << "Time cluster " << icl
                                      << " nhits=" << timeCluster.nhits()
                                      << " t0=" << timeCluster.t0().t0()
                                      << " t0err=" << timeCluster.t0().t0Err()
                                      << " pos=(" << timeCluster.position().X() << ","
                                      << timeCluster.position().Y() << ","
                                      << timeCluster.position().Z() << ")";
        for (auto hitIndex : timeCluster.hits()) {
          auto const& comboHit = comboHits.at(hitIndex);
          mf::LogVerbatim("PSegFinder") << "  hit " << hitIndex
                                        << " plane=" << comboHit.strawId().plane()
                                        << " panel=" << comboHit.strawId().panel()
                                        << " layer=" << comboHit.strawId().layer()
                                        << " straw=" << comboHit.strawId().straw()
                                        << " time=" << comboHit.time()
                                        << " driftTime=" << comboHit.driftTime()
                                        << " edep=" << comboHit.energyDep();
        }
      }
    }

    std::vector<TrkSegment*> segments;
    for (auto const& timeCluster : timeClusters) {
      for (auto hitIndex : timeCluster.hits()) {
        auto const& comboHit = comboHits.at(hitIndex);
        auto const strawId = comboHit.strawId();
        int const plane = strawId.plane();
        int const panel = strawId.panel();
        if (plane < 0 || plane >= 36 || panel < 0 || panel >= 6) continue;

        auto& segment = panelSegments[plane * 6 + panel];
        if (segment.nHits() == 0) segments.push_back(&segment);
        segment.fListOfHits.emplace_back(&comboHit);
      }
    }

    // Sort by straw and retain only the largest contiguous subsegment in each panel.
    struct SubSegment { int first; int last; int size() const { return last - first + 1; } };
    for (auto* segment : segments) {
      auto& hits = segment->fListOfHits;
      std::sort(hits.begin(), hits.end(), [](SegmentHit const& a, SegmentHit const& b) {
        return a.ComboHit()->strawId().straw() < b.ComboHit()->strawId().straw();
      });
      if (static_cast<int>(hits.size()) < minimumHits_) {
        segment->fMask |= TrkSegment::kNHitsBit;
        continue;
      }

      std::vector<SubSegment> subsegments;
      int first = 0;
      for (int i = 1; i < static_cast<int>(hits.size()); ++i) {
        int const gap = (hits[i].ComboHit()->strawId().straw()
                         - hits[i - 1].ComboHit()->strawId().straw() - 2) / 2;
        if (gap > maximumStrawGap_) {
          subsegments.push_back({first, i - 1});
          first = i;
        }
      }
      subsegments.push_back({first, static_cast<int>(hits.size()) - 1});

      auto const best = std::max_element(subsegments.begin(), subsegments.end(),
        [](SubSegment const& a, SubSegment const& b) { return a.size() < b.size(); });
      for (auto const& subsegment : subsegments) {
        if (&subsegment == &*best) continue;
        for (int i = subsegment.first; i <= subsegment.last; ++i) {
          hits[i].fMask |= TrkSegment::kSubsegmentBit;
        }
      }
    }

    if (debugBit_[10] > 0) {
      for (auto* segment : segments) {
        mf::LogVerbatim("PSegFinder") << "Panel candidate plane=" << segment->Plane()
                                      << " panel=" << segment->Panel()
                                      << " hits=" << segment->nHits();
        for (auto const& hit : segment->fListOfHits) {
          mf::LogVerbatim("PSegFinder") << "  straw=" << hit.ComboHit()->strawId().straw()
                                        << " mask=0x" << std::hex << hit.fMask << std::dec;
        }
      }
    }
//-----------------------------------------------------------------------------
// here the fit starts
//-----------------------------------------------------------------------------
    TrkSegment::Par_t fitParameters{};
    std::vector<TrkSegment*> fittedSegments;
    for (auto* segment : segments) {
      if (static_cast<int>(segment->nHits()) < minimumHits_) continue;
      if (segment->InitHits(nullptr) != 0) continue;

      SegmentFit fitter(segment);
      fitter.DefineDriftDirections();
      if (debugBit_[10] > 0) segment->print("after drift-sign assignment");
      fitter.Fit(4, 0, nullptr, &fitParameters);
      fittedSegments.push_back(segment);
    }

    if (debugBit_[10] > 0) {
      for (auto* segment : fittedSegments) segment->print("final fit result");
    }

    if (verbose_) {
      mf::LogInfo("PSegFinder")
        << "Reconstructed " << segments.size() << " panel candidates in event " << event.id();
    }
    // Intentionally no event.put(): segment persistence is deferred to a later version.
  }
};

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PSegFinder)
