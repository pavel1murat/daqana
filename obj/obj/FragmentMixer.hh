// A class that encapsulates event mixing code that is common for
// different Mu2e use cases.  It provides and registers callbacks to
// mix various data products, and is supposed to be used by MixFilter
// "detail" classes.  More documentation can be found in
//
//    art/Framework/Modules/MixFilter.h
//    art/Framework/IO/ProductMix/MixHelper.h
//    art/Framework/Core/PtrRemapper.h
//    art/Persistency/Common/CollectionUtilities.h
//
//
// cloned from : Andrei Gaponenko, 2018

#ifndef daqana_mod_FragmentMixer_hh
#define daqana_mod_FragmentMixer_hh

#include <string>
#include <vector>
#include <optional>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"

#include "artdaq-core/Data/Fragment.hh"

typedef std::vector<artdaq::Fragment> FragmentCollection;

//================================================================
namespace mu2e {

  class FragmentMixer {
  public:

    // Configuration for mixing one type of data products.
    struct CollectionMixerConfig {
      struct Entry {
        art::InputTag inTag;
        std::string outInstance;

        // Some outInstance inputs should not be treated literally,
        // this function is responsible for interpreting them.
        std::string resolvedInstanceName() const;

        Entry(const art::InputTag& i, const std::string& o): inTag(i), outInstance(o) {}
      };

      fhicl::Sequence<fhicl::TupleAs<Entry(art::InputTag,std::string)> >
      mixingMap { fhicl::Name("mixingMap"),
          fhicl::Comment("A sequence of InputTag to outputInstanceName"
                         " mappings for collections to be mixed."),
          std::vector<Entry>()
          };
    };

    // Configuration for the Mu2eProductMixing helper
    struct Config {
      fhicl::Table<CollectionMixerConfig> fragmentMixer { fhicl::Name("fragmentMixer") };
    };

    FragmentMixer(const Config& conf, art::MixHelper& helper);

    void startEvent(art::Event const& e);
    void processEventIDs(const art::EventIDSequence& seq);
    void beginSubRun(const art::SubRun& sr);
    void endSubRun(art::SubRun& sr);

  private:

    bool mixFragments(std::vector<FragmentCollection const*> const& in,
                      FragmentCollection& out,
                      art::PtrRemapper const& remap);

    bool mixEventIDs(std::vector<art::EventIDSequence const*> const &in,
                     art::EventIDSequence& out,
                     art::PtrRemapper const& remap);
    

    unsigned int generatedEvents_ = 0;
    unsigned int resampledEvents_ = 0;
  };

}

#endif
