///////////////////////////////////////////////////////////////////////////////
// SimHits: simulate one straight track per event through a single tracker panel and
//          put the straw hits of the track into the event as a ComboHitCollection.
//          the reconstruction of the segment is done by the SimRes module
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include <memory>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "SimPanel.hh"

namespace mu2e {

class SimHits : public art::EDProducer {
public:
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int>    debugMode  {Name("debugMode"  ), Comment("debug bits, see SimPanel.hh")};
    fhicl::Atom<double> strawRadius{Name("strawRadius"), Comment("straw radius, mm")};
    fhicl::Atom<double> vDrift     {Name("vDrift"     ), Comment("drift velocity, mm/ns")};
    fhicl::Atom<int>    beamMode   {Name("beamMode"   ), Comment("1: line parameters of the odd-reco study")};
    fhicl::Sequence<double> strawOffsetX{Name("strawOffsetX"), Comment("misalignment of the straws in X, 96 values, mm")};
    fhicl::Sequence<double> strawOffsetY{Name("strawOffsetY"), Comment("misalignment of the straws in Y, 96 values, mm")};
  };
  using Parameters = art::EDProducer::Table<Config>;

  explicit SimHits(Parameters const& config);

  void produce(art::Event& event) override;

private:
  int    fDebugMode;
  double fStrawRadius;
  double fVDrift;
  int    fBeamMode;
  std::vector<double> fStrawOffsetX;
  std::vector<double> fStrawOffsetY;
};

//-----------------------------------------------------------------------------
SimHits::SimHits(Parameters const& config) : art::EDProducer(config) {
  fDebugMode   = config().debugMode();
  fStrawRadius = config().strawRadius();
  fVDrift      = config().vDrift();
  fBeamMode    = config().beamMode();
  fStrawOffsetX = config().strawOffsetX();
  fStrawOffsetY = config().strawOffsetY();

  produces<mu2e::ComboHitCollection>();
}

//-----------------------------------------------------------------------------
void SimHits::produce(art::Event& event) {
  auto const straws = simpanel::straws(fDebugMode,fStrawOffsetX,fStrawOffsetY);
  auto const line   = simpanel::line(fBeamMode,fDebugMode);
  auto const hits   = simpanel::hitStraws(straws,line,fStrawRadius,fDebugMode);

  auto chc = std::make_unique<mu2e::ComboHitCollection>();
  for (size_t i=0; i<hits.size(); ++i) {
    simpanel::Straw_t const& s = hits[i];
                                        // the second hit has an extra 10 ns, as before
    double sig_t0 = (i == 1) ? 10. : 0.;
    double tmeas  = line.fT0+s.fDoca/fVDrift+sig_t0;
    if (tmeas < 0) tmeas = 0;

    mu2e::ComboHit ch;
    ch._sid                        = s.id;
    ch._eend                       = mu2e::StrawEnd::cal;
    ch._etime[mu2e::StrawEnd::cal] = tmeas;
    ch._ptime                      = 0;
    ch._pos                        = CLHEP::Hep3Vector(0,s.fX,s.fY);
    chc->push_back(ch);
  }

  event.put(std::move(chc));
}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimHits)
