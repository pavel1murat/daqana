// ======================================================================
//
// MakeDigiNtuple:  assume that digis have been produced, make a hit ntuple
//                  for the cross-subsystem timing studies
// tracker       :  assume that the waveforms are made
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"

#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include "otsdaq-mu2e-tracker/Nt/DaqEvent.hh"
#include "otsdaq-mu2e-tracker/Nt/DaqStrawDigi.hh"

#include <iostream>
#include <string>
#include <memory>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
// #define TRACEMF_USE_VERBATIM 1

// #include "TRACE/tracemf.h"
// #define TRACE_NAME "MakeDigiNtuple"


namespace mu2e {
  class MakeDigiNtuple;
}

// ======================================================================
class mu2e::MakeDigiNtuple : public art::EDAnalyzer {

public:

  DaqEvent*         _event;
  const art::Event* _art_event;

  int              _nstrawdigis;
  int              _ncalodigis;
  int              _ncrvdigis;
  int              _nstmdigis;

  TFile*           _file;
  TTree*           _tree;
  TBranch*         _branch;

  struct Config {
    fhicl::Atom<art::InputTag>  sdCollTag    {fhicl::Name("sdCollTag"    ), fhicl::Comment("straw digi coll tag"       ),"undefined"};
    fhicl::Atom<int>            debugMode    {fhicl::Name("debugMode"    ), fhicl::Comment("debug mode"                ),0};
    fhicl::Atom<int>            diagLevel    {fhicl::Name("diagLevel"    ), fhicl::Comment("diagnostic level"          ),0};
    fhicl::Atom<std::string>    outputDir    {fhicl::Name("outputDir"    ), fhicl::Comment("output directory"          ),"./"};
    fhicl::Atom<int>            saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms"),0};
    fhicl::Atom<int>            ewLength     {fhicl::Name("ewLength"     ), fhicl::Comment("event window length, in units of 25 ns"),1000};
  };

  // --- C'tor/d'tor:
  explicit MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config);
  virtual ~MakeDigiNtuple() {}

  int      getData(const art::Event& ArtEvent);
  
  void     print_(const std::string&  Message, int DiagLevel = -1,
                  const std::source_location& location = std::source_location::current());

//-----------------------------------------------------------------------------
// overloaded virtual functions of EDAnalyzer
//-----------------------------------------------------------------------------
  virtual void beginRun(const art::Run&   r);
  virtual void analyze (const art::Event& e);
  virtual void endRun  (const art::Run&   r);
  virtual void beginJob();
  virtual void endJob  ();

  int              debugMode_;
  int              diagLevel_;
  art::InputTag    sdCollTag_;          // straw digi collection tag
  std::string      outputDir_;
  int              saveWaveforms_;
  
  int              n_adc_samples_;
  int              ewLength_;           // it is up to the user to make sure it is set correctly
  double           tdc_bin_;            // TDC bin, in us
  double           tdc_bin_ns_;         // TDC bin, ns

  const mu2e::StrawDigiCollection*             _sdc;
  const mu2e::StrawDigiADCWaveformCollection*  _sdawfc;

}; // MakeDigiNtuple

// ======================================================================

mu2e::MakeDigiNtuple::MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    diagLevel_    (config().diagLevel    ()),
    sdCollTag_    (config().sdCollTag    ()),
    outputDir_    (config().outputDir    ()),
    saveWaveforms_(config().saveWaveforms()),
    ewLength_     (config().ewLength     ())
{
  if (saveWaveforms_ == 0) n_adc_samples_ = 0;

  tdc_bin_             = (5/256.*1e-3);       // TDC bin width (Richie), in us
  tdc_bin_ns_          = tdc_bin_*1e3;        // convert to ns

}


//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::print_(const std::string& Message, int DiagLevel,
                                 const std::source_location& location) {
  
  if (DiagLevel > diagLevel_) return;
  std::cout << std::format(" event:{}:{}:{}",_art_event->run(),_art_event->subRun(),_art_event->event())
            << " " << location.file_name() << ":" << location.line()
    //            << location.function_name()
            << ": " << Message << std::endl;
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::beginRun(const art::Run& ArtRun) {

  art::ServiceHandle<art::TFileService> tfs;
  TH1::AddDirectory(kFALSE);

  
  //  _file = new TFile(Form("%s/make_digi_ntuple_%06i.root",outputDir_.data(),ArtRun.run()),"RECREATE");
  TTree::SetMaxTreeSize(8000000000LL);

  //  _tree = new TTree("digis","digis");
  _tree = tfs->make<TTree>("digis","digis");

  _event = new DaqEvent();
  
  _branch = _tree->Branch("evt","DaqEvent",_event,32000,99);
  _branch->SetAutoDelete(kFALSE);

  if (_branch) { 
    // _event->strawdigis = new TClonesArray("DaqStrawDigi",100);
    // _event->strawdigis->BypassStreamer(kFALSE);          // the whole point is to split everything
    // _event->calodigis  = new TClonesArray("DaqStrawDigi",100);
    // _event->crvdigis   = new TClonesArray("DaqStrawDigi",100);
    // _event->stmdigis   = new TClonesArray("DaqStrawDigi",100);
  }
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::endRun(const art::Run& ArtRun) {
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::beginJob() {
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::endJob() {
  // _file->Write();
  // _file->Close();
  // delete _file;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::getData(const art::Event& ArtEvent) {
  int rc(0);

//-----------------------------------------------------------------------------
// tracker
//-----------------------------------------------------------------------------
  art::Handle<mu2e::StrawDigiCollection>            sdch;
  art::Handle<mu2e::StrawDigiADCWaveformCollection> sdawfch;

  _nstrawdigis = 0;
  // _ncalodigis  = 0;
  // _ncrvdigis   = 0;
  // _nstmdigis   = 0;

  _sdc         = nullptr;
  _sdawfc      = nullptr;
  
  bool ok = ArtEvent.getByLabel(sdCollTag_,sdch);
  if (ok) { 
    _sdc         = sdch.product();
    _nstrawdigis = _sdc->size();
  }
  else {
    print_(std::format("ERROR: StrawDigiCollection:{:s} is not available. Bail out",sdCollTag_.encode().data()));
    return -1;
  }

  ok =  ArtEvent.getByLabel(sdCollTag_,sdawfch);
  if (ok) { 
    _sdawfc = sdawfch.product();
  }
  else {
    print_(std::format("WARNING: StrawDigiADCWaveformCollection:{:s} is not available. Bail out",
                       sdCollTag_.encode().data()),1);
    return -1;
  }
//-----------------------------------------------------------------------------
// calorimeter
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// STM
//-----------------------------------------------------------------------------

  return rc;
}

// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::analyze(const art::Event& ArtEvent) {
  //  int const packet_size(16); // in bytes

  _art_event = &ArtEvent;

  if (debugMode_ > 0) {
    print_(std::format("-- START event:{}:{}:{}",ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event()),1);
  }

  int rc = getData(ArtEvent);
  if (rc < 0) return;
//-----------------------------------------------------------------------------
// clear event
//-----------------------------------------------------------------------------
  _event->nsd = _nstrawdigis;
  // _event->ncalodigis  = _ncalodigis;
  // _event->ncrvdigis   = _ncrvdigis;
  // _event->nstmdigis   = _nstmdigis;
  //  _event->sdgis->Clear();
  // _event->calodigis->Clear();
  // _event->crvdigis->Clear();
  // _event->stmdigis->Clear();
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
  _event->run = ArtEvent.run();
  _event->srn = ArtEvent.subRun();
  _event->evn = ArtEvent.event();
  _event->sd->Clear();
  
  DaqStrawDigi nt_sd;
  //   DaqEvent::fgSd.clear();
  if (debugMode_ > 0) {
    print_(std::format("_nstrawdigis:{}",_nstrawdigis),1);
  }
  for (int i=0; i<_nstrawdigis; i++) {
    const mu2e::StrawDigi*            sd    = &_sdc->at(i);
    const mu2e::StrawDigiADCWaveform* sdawf = &_sdawfc->at(i);
    int ns = sdawf->samples().size();

    DaqStrawDigi* nt_sd = new ((*_event->sd)[i]) DaqStrawDigi(ns);
    nt_sd->sid          = sd->strawId().asUint16();
    nt_sd->tdc0         = sd->TDC(mu2e::StrawEnd::cal);
    nt_sd->tdc1         = sd->TDC(mu2e::StrawEnd::hv );
    nt_sd->tot0         = sd->TOT(mu2e::StrawEnd::cal);
    nt_sd->tot1         = sd->TOT(mu2e::StrawEnd::hv );
    nt_sd->pmp          = sd->PMP();
    nt_sd->flag         = *((uint8_t*) &sd->digiFlag());
//-----------------------------------------------------------------------------
// store the waveform
//-----------------------------------------------------------------------------
    // for (int is=0; is<ns; is++) {
    //   nt_sd.adc[is] = sdawf->samples()[is];
    // }

    if (debugMode_  > 0) {
      if (diagLevel_ == 11) {
//-----------------------------------------------------------------------------
// for all hits, print hit times assuming contiguous timing
//-----------------------------------------------------------------------------
        double t0_offset  = _event->evn*ewLength_*25;                  // in ns
        double t0         = t0_offset + nt_sd->tdc0*tdc_bin_ns_;
        double t1         = t0_offset + nt_sd->tdc1*tdc_bin_ns_;
        printf("%8i %5i %8i %8i %12.4lf %12.4lf %6i %6i %6i 0x%04x\n",
               _event->evn,
               (int) nt_sd->sid,
               nt_sd->tdc0, nt_sd->tdc1,
               t0, t1,
               nt_sd->tot0, nt_sd->tot1,
               nt_sd->pmp, nt_sd->flag);
      }
    }
  }

  //  _event->sd = DaqEvent::fgSd.data();

  if (debugMode_ > 0) print_(std::format("_event->strawdigis->GetEntries():{}",_event->nsd),1);

  _tree->Fill();

  if (debugMode_ > 0) print_("-- END",1);
}



// ======================================================================

DEFINE_ART_MODULE(mu2e::MakeDigiNtuple)

// ======================================================================
