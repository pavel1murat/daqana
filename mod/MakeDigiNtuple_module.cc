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
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Offline/DAQ/inc/TrkPanelMap_t.hh"

#include "daqana/obj/DaqEvent.hh"
#include "daqana/obj/DaqStrawDigi.hh"

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

  struct Config {
    fhicl::Atom<art::InputTag>   sdCollTag    {fhicl::Name("sdCollTag"    ), fhicl::Comment("straw digi coll tag"       ),""};
    fhicl::Atom<art::InputTag>   shCollTag    {fhicl::Name("shCollTag"    ), fhicl::Comment("straw hit  coll tag"       ),""};
    fhicl::Atom<int>             debugMode    {fhicl::Name("debugMode"    ), fhicl::Comment("debug mode"                ),0};
    fhicl::Sequence<std::string> debugBits    {fhicl::Name("debugBits"    ), fhicl::Comment("debug bits"                )};
    fhicl::Atom<std::string>     outputDir    {fhicl::Name("outputDir"    ), fhicl::Comment("output directory"          )};
    fhicl::Atom<int>             saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms"),0};
    fhicl::Atom<int>             makeSD       {fhicl::Name("makeSD"       ), fhicl::Comment("make straw digi branch"     ),1};
    fhicl::Atom<int>             makeSH       {fhicl::Name("makeSH"       ), fhicl::Comment("make straw hit branch"      ),1};
    fhicl::Atom<int>             ewLength     {fhicl::Name("ewLength"     ), fhicl::Comment("event window length, in units of 25 ns"),1000};
  };

  // --- C'tor/d'tor:
  explicit MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config);
  virtual ~MakeDigiNtuple() {}

  int      getData(const art::Event& ArtEvent);
  
  void     print_(const std::string&  Message, const std::source_location& location = std::source_location::current());

  int      fillSD();
  int      fillSH();
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDAnalyzer
//-----------------------------------------------------------------------------
  virtual void beginRun(const art::Run&   r);
  virtual void analyze (const art::Event& e);
  virtual void endRun  (const art::Run&   r);
  virtual void beginJob();
  virtual void endJob  ();

  int                      _debugMode;
  std::vector<std::string> _debugBits;
  int                      _debugBit[100];
  art::InputTag            _sdCollTag;          // straw digi collection tag
  art::InputTag            _shCollTag;          // straw hit  collection tag
  std::string              _outputDir;
  int                      _saveWaveforms;
  int                      _makeSD;
  int                      _makeSH;
  int                      _ewLength;           // it is up to the user to make sure it is set correctly
  
  int                      _n_adc_samples;
  double                   _tdc_bin;            // TDC bin, in us
  double                   _tdc_bin_ns;         // TDC bin, ns

  DaqEvent*                _event;
  const art::Event*        _art_event;

  int                     _nstrawdigis;
  int                     _nstrawhits;
  int                     _ncalodigis;
  int                     _ncrvdigis;
  int                     _nstmdigis;
                          
  TFile*                  _file;
  TTree*                  _tree;
  TBranch*                _branch;
  const TrkPanelMap_t*    _panel_map[36][6];   // indexing - offline: [plane][panel]

  const mu2e::StrawDigiCollection*             _sdc;
  const mu2e::StrawDigiADCWaveformCollection*  _sdawfc;
  const mu2e::StrawHitCollection*              _shc;
  
}; // MakeDigiNtuple

// ======================================================================

mu2e::MakeDigiNtuple::MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _debugMode    (config().debugMode    ()),
    _debugBits    (config().debugBits    ()),
    _sdCollTag    (config().sdCollTag    ()),
    _shCollTag    (config().shCollTag    ()),
    _outputDir    (config().outputDir    ()),
    _saveWaveforms(config().saveWaveforms()),
    _makeSD       (config().makeSD       ()),
    _makeSH       (config().makeSH       ()),
    _ewLength     (config().ewLength     ()),
    _art_event    (nullptr)
{
  if (_saveWaveforms == 0) _n_adc_samples = 0;

  _tdc_bin             = (5/256.*1e-3);       // TDC bin width (Richie), in us
  _tdc_bin_ns          = _tdc_bin*1e3;        // convert to ns
//-----------------------------------------------------------------------------
// parse debug bits
//-----------------------------------------------------------------------------
  const char* key;
                                        // a flag is an integer!
  int nbits = _debugBits.size();
  for (int i=0; i<nbits; i++) {
    int index(0), value(0);
    key               = _debugBits[i].data();
    sscanf(key,"bit%i:%i",&index,&value);
    _debugBit[index]  = value;
    
    print_(std::format("... StrawDigisFromArtdaqFragments: bit={:4d} is set to {}\n",index,_debugBit[index]));
  }
//-----------------------------------------------------------------------------
// initialize the panel map
//-----------------------------------------------------------------------------
  for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
    int plane = tpm->plane;
    int panel = tpm->panel;
    _panel_map[plane][panel] = tpm;
  }
}


//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::print_(const std::string& Message, const std::source_location& location) {
  if (_art_event) {
    std::cout << std::format(" event:{}:{}:{}",
                             _art_event->run(),_art_event->subRun(),_art_event->event());
  }
  std::cout << " " << location.file_name() << ":" << location.line()
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
  art::Handle<mu2e::StrawHitCollection>             shch;

  _nstrawdigis = 0;
  _nstrawhits  = 0;
  // _ncalodigis  = 0;
  // _ncrvdigis   = 0;
  // _nstmdigis   = 0;

  _sdc         = nullptr;
  _sdawfc      = nullptr;
  _shc         = nullptr;
  
  bool ok = ArtEvent.getByLabel(_shCollTag,shch);
  if (ok) { 
    _shc         = shch.product();
    _nstrawhits = _shc->size();
  }
  else {
    print_(std::format("ERROR: StrawHitCollection:{:s} is not available. Bail out",_shCollTag.encode().data()));
    return -1;
  }

  ok = ArtEvent.getByLabel(_sdCollTag,sdch);
  if (ok) { 
    _sdc         = sdch.product();
    _nstrawdigis = _sdc->size();
  }
  else {
    print_(std::format("ERROR: StrawDigiCollection:{:s} is not available. Bail out",_sdCollTag.encode().data()));
    return -1;
  }

  ok =  ArtEvent.getByLabel(_sdCollTag,sdawfch);
  if (ok) { 
    _sdawfc = sdawfch.product();
  }
  else {
    print_(std::format("WARNING: StrawDigiADCWaveformCollection:{:s} is not available. Bail out",
                       _sdCollTag.encode().data()));
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

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSD() {
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

    if (_debugMode  > 0) {
      if (_debugBit[1] != 0) {
//-----------------------------------------------------------------------------
// for all hits, print hit times assuming contiguous timing
//-----------------------------------------------------------------------------
        double t0_offset  = _event->evn*_ewLength*25;                  // in ns
        double t0         = t0_offset + nt_sd->tdc0*_tdc_bin_ns;
        double t1         = t0_offset + nt_sd->tdc1*_tdc_bin_ns;
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
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSH() {
  for (int i=0; i<_nstrawhits; i++) {
    const mu2e::StrawHit* sh = &_shc->at(i);
    int pln = sh->strawId().plane();
    int pnl = sh->strawId().panel();
    const TrkPanelMap_t* tpm = _panel_map[pln][pnl];
    
   
    DaqStrawHit* nt_sh = new ((*_event->sh)[i]) DaqStrawHit();
    nt_sh->sid         = sh->strawId().asUint16();
    nt_sh->zface       = tpm->zface;
    nt_sh->mnid        = tpm->mnid;
    nt_sh->time        = sh->time(mu2e::StrawEnd::cal);
    nt_sh->dt          = sh->dt();
    nt_sh->tot0        = sh->TOT(mu2e::StrawEnd::cal);
    nt_sh->tot1        = sh->TOT(mu2e::StrawEnd::hv );
    nt_sh->edep        = sh->energyDep();

    if (_debugMode  > 0) {
      if (_debugBit[1] != 0) {
        printf("%8i %5i %12.4f %12.4f %6.1f %6.1f %7.4f\n",
               _event->evn,
               (int) nt_sh->sid,
               nt_sh->time, nt_sh->dt,
               nt_sh->tot0, nt_sh->tot1,
               nt_sh->edep);
      }
    }
  }
  return 0;
}

// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::analyze(const art::Event& ArtEvent) {
  //  int const packet_size(16); // in bytes

  _art_event = &ArtEvent;

  if (_debugMode > 0) {
    print_(std::format("-- START event:{}:{}:{}",ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event()));
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
  
  //   DaqEvent::fgSd.clear();
  if (_debugMode > 0) {
    print_(std::format("_nstrawdigis:{}",_nstrawdigis));
  }

  if (_makeSD) fillSD();
  if (_makeSH) fillSH();

  //  _event->sd = DaqEvent::fgSd.data();

  if (_debugMode > 0) {
    print_(std::format("_event->strawdigis->GetEntries():{}",_event->nsd));
  }

  _tree->Fill();

  if (_debugMode > 0) print_("-- END");
}



// ======================================================================

DEFINE_ART_MODULE(mu2e::MakeDigiNtuple)

// ======================================================================
