// This module produces "background frames" that mimic Mu2e
// microbunch-events from single particle simualtion inputs.  The
// number of particles to mix is determined based on the input
// ProtonBunchIntensity object that models beam intensity
// fluctuations.  There is a random Poisson process that is on top of
// the beam intensity fluctuations, which represents the probability
// of a secondary from a given proton creating a hit in a collection
// to be mixed.  This Poisson is sampled by the module.
//
// Andrei Gaponenko, 2018

#include <random>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art_root_io/RootIOPolicy.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"

#include "artdaq-core/Data/Fragment.hh"

#include "daqana/obj/FragmentMixer.hh"
//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  // Our "detail" class for art/Framework/Modules/MixFilter.h
  class MixFragmentsDetail {
    FragmentMixer spm_;
    const int debugLevel_;
    const unsigned maxEventsToSkip_;


    bool writeEventIDs_;
    art::EventIDSequence idseq_;

  public:

    struct FragmentMixerConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<FragmentMixer::Config> products { Name("products"),
          Comment("A table specifying products to be mixed.  For each supported data type\n"
                  "there is a mixingMap sequence that defines mapping of inputs to outputs.\n"
                  "Each entry in the top-level mixingMap sequence is a sequence of two strings:\n"
                  "    [ \"InputTag\", \"outputInstanceName\" ]\n"
                  "The output instance name colon \":\" is special: it means take instance name from the input tag.\n"
                  "For example, with this config:\n"
                  "   mixingMap: [ [ \"detectorFilter:tracker\", \"tracker\" ], [ \"detectorFilter:virtualdetector\", \":\" ] ]\n"
                  "the outputs will be named \"tracker\" and \"virtualdetector\"\n"
                  )
          };

      fhicl::Atom<int> debugLevel { Name("debugLevel"),
          Comment("control the level of debug output"),
          0u
          };
      fhicl::Atom<unsigned> maxEventsToSkip { Name("MaxEventsToSkip"),
          Comment("Maximum number of events to skip at the beginning of the first secondary input file in sequential readMode.\n"),
          0u };

      fhicl::Atom<bool> writeEventIDs { Name("writeEventIDs"),
          Comment("Write out IDs of events on the secondary input stream."),
          false
          };

    };

    // The ".mu2e" in FHICL parameters like
    // physics.filters.somemixer.mu2e.meanEventsPerProton clearly
    // separates experiment specific settings from those provided by
    // the art framework (like "somemixer.wrapFiles").
    struct Config {
      fhicl::Table<FragmentMixerConfig> mu2e { fhicl::Name("mu2e") };
    };

    using Parameters = art::MixFilterTable<Config>;
    explicit MixFragmentsDetail(const Parameters& pars, art::MixHelper& helper);


    size_t nSecondaries();
    size_t eventsToSkip();

    void processEventIDs(const art::EventIDSequence& seq);

    void beginSubRun(const art::SubRun& sr);
    void startEvent(const art::Event& evt);
    void finalizeEvent(art::Event& e);
    void endSubRun(art::SubRun& sr);

  };

  //================================================================
  MixFragmentsDetail::MixFragmentsDetail(const Parameters& pars, art::MixHelper& helper)
    : spm_{ pars().mu2e().products(), helper }
    , debugLevel_{ pars().mu2e().debugLevel() }
    , maxEventsToSkip_{ pars().mu2e().maxEventsToSkip() }
    , writeEventIDs_{ pars().mu2e().writeEventIDs() }
  {
    if(writeEventIDs_) {
      helper.produces<art::EventIDSequence>();
    }
  }

  //================================================================
  void MixFragmentsDetail::beginSubRun(const art::SubRun& sr) {
    spm_.beginSubRun(sr);
  }

  //================================================================
  void MixFragmentsDetail::endSubRun(art::SubRun& sr) {
    spm_.endSubRun(sr);
  }

  //================================================================
  void MixFragmentsDetail::startEvent(const art::Event& event) {
  // call down to product mixer
    spm_.startEvent(event);

  }

  //================================================================
  size_t MixFragmentsDetail::eventsToSkip() {
    return 0;
  }

  //================================================================
  size_t MixFragmentsDetail::nSecondaries() {
    return 1;
  }

  //================================================================
  void MixFragmentsDetail::processEventIDs(art::EventIDSequence const& seq) {

    spm_.processEventIDs(seq);

    if(writeEventIDs_) {
      idseq_ = seq;
    }

  }

  //================================================================
  void MixFragmentsDetail::finalizeEvent(art::Event& e) {
    std::cout << std::format(" MixFragmentsDetail::finalizeEvent: finalizing: writeEventIDs_:{}\n",writeEventIDs_);
    if(writeEventIDs_) {
      auto o = std::make_unique<art::EventIDSequence>();
      o->swap(idseq_);
      e.put(std::move(o));
    }
  }

  //================================================================

  //================================================================
  // This is the module class.
  typedef art::MixFilter<MixFragmentsDetail,art::RootIOPolicy> MixFragments;
}

DEFINE_ART_MODULE(mu2e::MixFragments)
