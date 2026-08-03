// from: Andrei Gaponenko, 2018

#include "daqana/obj/FragmentMixer.hh"

#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include "cetlib_except/exception.h"

#include "art/Persistency/Common/CollectionUtilities.h"

//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  std::string FragmentMixer::CollectionMixerConfig::Entry::resolvedInstanceName() const {
    return (outInstance == ":") ? inTag.instance() :  outInstance;
  }

  //----------------------------------------------------------------
  namespace {

    // newEntryIndex is an address in the output (flattened)
    // collection and the offsets are those that were recorded by the
    // flattenCollections() call
    template<typename IndexType, typename OFFSETS>
    typename OFFSETS::size_type
    getInputEventIndex(IndexType newEntryIndex, const OFFSETS& offsets) {
      auto ub = std::upper_bound(offsets.begin(), offsets.end(), newEntryIndex);
      if(ub == offsets.begin()) {
        throw cet::exception("RANGE")<<"getInputEventIndex(): newEntryIndex="
                                     <<newEntryIndex
                                     <<" is below the first offset="
                                     <<*offsets.begin()
                                     <<std::endl;
      }
      return std::distance(offsets.begin(), --ub);
    }
  }

  //----------------------------------------------------------------
  FragmentMixer::FragmentMixer(const Config& conf, art::MixHelper& helper) {
    for(const auto& e: conf.fragmentMixer().mixingMap()) {
      helper.declareMixOp(e.inTag, e.resolvedInstanceName(), &FragmentMixer::mixFragments, *this);
    }
  }


  //================================================================
  void FragmentMixer::startEvent(art::Event const& e) {
    resampledEvents_++;
  }

  //----------------------------------------------------------------
  void FragmentMixer::processEventIDs(const art::EventIDSequence& seq)  {
  }

  //----------------------------------------------------------------
  void FragmentMixer::beginSubRun(const art::SubRun& s) {

  }

  //----------------------------------------------------------------
  void FragmentMixer::endSubRun(art::SubRun& sr) {
  }

 

  //-----------------------------------------------------------------------------  
  bool FragmentMixer::mixFragments(std::vector<FragmentCollection const*> const& in,
                                   FragmentCollection& out,
                                   art::PtrRemapper const& remap)
  {
    std::cout << std::format(" mixFragments: flattening : size:{}\n",in.size());
    int ncolls = in.size();
    for (int i=0; i<ncolls; i++) {
      const std::vector<artdaq::Fragment>* coll = in[i];
      int nf = -1;
      if (coll != nullptr) nf = coll->size();
      std::cout << std::format("FragmentMixer::mixFragments: coll->size():{}\n",nf);
    }
    art::flattenCollections(in, out);
    return true;
  }


  //----------------------------------------------------------------
  bool FragmentMixer::mixEventIDs(std::vector<art::EventIDSequence const*> const &in,
                                     art::EventIDSequence& out,
                                     art::PtrRemapper const& remap)
  {
    std::cout << std::format(" mixEventIDs: flattening\n");
    
    art::flattenCollections(in, out);
    return true;
  }

}
//================================================================
