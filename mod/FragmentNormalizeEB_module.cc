///////////////////////////////////////////////////////////////////////////////
// FragmentNormalizeEB_module.cc, mostly writen by codex-5.3
///////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "artdaq-core/Data/Fragment.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <memory>
#include <string>
#include <vector>
#include <cstdio>

namespace mu2e {

class FragmentNormalizeEB : public art::EDProducer {
public:
  struct Config {
    fhicl::Atom<std::string>     moduleLabel{fhicl::Name("moduleLabel"), "daq"};
    fhicl::Atom<std::string>     ebLabel    {fhicl::Name("ebLabel"), "eb"};
    fhicl::Sequence<std::string> instances{fhicl::Name("instances")};
    fhicl::Atom<bool>            allowMissing{fhicl::Name("allowMissing"), true};
    fhicl::Atom<bool>            verbose{fhicl::Name("verbose"), false};
  };
  using Parameters = art::EDProducer::Table<Config>;

  explicit FragmentNormalizeEB(Parameters const& cfg)
    : art::EDProducer(cfg)
    , moduleLabel_ (cfg().moduleLabel())
    , ebLabel_     (cfg().ebLabel())
    , instances_   (cfg().instances())
    , allowMissing_(cfg().allowMissing())
    , verbose_     (cfg().verbose())
  {
    for (int i = 1; i <= 17; ++i) {
      std::string b = std::format("{}{:02d}", ebLabel_.data(),i);
      processes_.emplace_back(b);
    }

    tags_.resize(instances_.size());
    for (size_t iInst = 0; iInst < instances_.size(); ++iInst) {
      for (auto const& p : processes_) {
        std::cout << std::format("process: {}\n",p);
        tags_[iInst].emplace_back(moduleLabel_, instances_[iInst], p);
      }
      produces<std::vector<artdaq::Fragment>>(instances_[iInst]);
    }
  }

  void produce(art::Event& evt) override {
    for (size_t iInst = 0; iInst < instances_.size(); ++iInst) {
      bool found = false;
      size_t foundProc = 0;
      art::Handle<std::vector<artdaq::Fragment>> h;

      for (size_t ip = 0; ip < processes_.size(); ++ip) {
        evt.getByLabel(tags_[iInst][ip], h);
        if (h.isValid()) {
          found = true;
          foundProc = ip;
          break;
        }
      }

      auto out = std::make_unique<std::vector<artdaq::Fragment>>();
      if (found) {
        out->insert(out->end(), h->begin(), h->end());
      } else if (!allowMissing_) {
        throw cet::exception("FragmentNormalizeEB")
          << "Missing " << moduleLabel_ << ":" << instances_[iInst]
          << " for eb01..eb17 in event " << evt.id();
      }

      evt.put(std::move(out), instances_[iInst]);

      if (verbose_) {
        mf::LogInfo("FragmentNormalizeEB")
          << "Event " << evt.id()
          << " instance " << instances_[iInst]
          << (found ? " from " + processes_[foundProc] : " missing (wrote empty)");
      }
    }
  }

private:
  std::string moduleLabel_;
  std::string ebLabel_;
  std::vector<std::string> instances_;
  bool allowMissing_;
  bool verbose_;

  std::vector<std::string> processes_;                    // eb01..eb17
  std::vector<std::vector<art::InputTag>> tags_;          // [instance][process]
};

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FragmentNormalizeEB)
