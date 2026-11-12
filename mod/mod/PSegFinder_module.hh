#ifndef __daqana_mod_PSegFinder_module_hh__
#define __daqana_mod_PSegFinder_module_hh__

// Reconstruct per-panel tracker segments from ComboHits assigned to time clusters.
// The segments are event-local in this first version and are not persisted.

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/DelegatedParameter.h"
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

#include "TFolder.h"

#include "daqana/mod/PSegFinderHolder.hh"

#include "Stntuple/mod/TAnaRint.hh"
#include "Stntuple/print/TAnaDump.hh"

namespace mu2e {

class PSegFinder : public art::EDProducer {
public:
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag>   comboHitTag    {Name("comboHitTag"    ), Comment("ComboHit collection")};
    fhicl::Atom<art::InputTag>   tcCollTag      {Name("tcCollTag"      ), Comment("TimeCluster collection tag")};
    fhicl::Atom<double>          driftVelocity  {Name("driftVelocity"  ), Comment("Straw drift velocity in mm/ns")};
    fhicl::Atom<double>          timeOffset     {Name("timeOffset"     ), Comment("Drift-time offset in ns")};
    fhicl::Atom<int>             minimumHits    {Name("minimumHits"    ), Comment("Minimum hits required to fit a segment")};
    fhicl::Atom<int>             maximumStrawGap{Name("maximumStrawGap"), Comment("Largest gap between hits in one subsegment")};
    fhicl::Atom<int>             interactiveMode{Name("interactiveMode"), Comment("0: none, 1: stop after each event; 2: stop at endJob")};
    fhicl::Atom<bool>            verbose        {Name("verbose"        ), Comment("Print the number of reconstructed panel segments")};
    fhicl::Sequence<std::string> debugBits      {Name("debugBits"      ), Comment("Debug switches in bitN:value form")};
    fhicl::DelegatedParameter    anaDump        {Name("TAnaDump"       ), Comment("TAnaDump parameters"  ) };
  };
  using Parameters = art::EDProducer::Table<Config>;

private:
//-----------------------------------------------------------------------------
// FCL parameters
//-----------------------------------------------------------------------------
  art::InputTag             comboHitTag_;
  art::InputTag             tcCollTag_;
  int                       minimumHits_;
  int                       maximumStrawGap_;
  int                       interactiveMode_;
  bool                      verbose_;
  ProditionsHandle<Tracker> alignedTracker_;
  int                       debugBit_[100]{};
//-----------------------------------------------------------------------1------
// the rest
//-----------------------------------------------------------------------------
  int                       runNumber_{-1};
  TAnaRint*                 fAnaRint;
  TAnaDump*                 anaDump_;

  TFolder*                  folder_;

  std::vector<TCanvas*>     canvases_;          // canvases of the segment plots, closed at the end of the event

  std::vector<TrkSegment>   panelSegments_;     // static allocation, to be done once
  std::vector<TrkSegment*>  segments_;          // segments with hits
  const Tracker*            tracker_;

  PSegFinderHolder*         moduleHolder_;

  const TimeClusterCollection* tcColl_;
  const ComboHitCollection*    chColl_;

public:
  explicit      PSegFinder(Parameters const& config);

  virtual void  beginRun(art::Run&   ArtRun  ) override;
  virtual void  produce (art::Event& ArtEvent) override;

  void           print(int DebugBit);

  TrkSegment*    segment(int I) { return segments_[I]; }

  std::vector<TrkSegment*>* segments() { return &segments_; }

  void           plot_segment(int I);
  void           closeCanvases();

  const Tracker* tracker() { return tracker_; }

};

} // namespace mu2e
#endif
