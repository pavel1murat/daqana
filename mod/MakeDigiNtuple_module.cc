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
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "daqana/mod/TrkPanelMap_t.hh"

#include "daqana/obj/DaqEvent.hh"
#include "daqana/mod/WfParam_t.hh"

#include <regex>
#include <ranges>
#include <string>
#include <vector>
#include <algorithm>

#include <iostream>
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

using namespace fhicl;

// ======================================================================
class mu2e::MakeDigiNtuple : public art::EDAnalyzer {

public:

  struct Config {
    
    Atom<art::InputTag>   sdCollTag     {Name("sdCollTag"     ), Comment("straw digi coll tag"       ),""};
    Atom<art::InputTag>   shCollTag     {Name("shCollTag"     ), Comment("straw hit  coll tag"       ),""};
    Atom<art::InputTag>   tcCollTag     {Name("tcCollTag"     ), Comment("time cluster coll tag"     ),""};
    Atom<int>             debugMode     {Name("debugMode"     ), Comment("debug mode"                ),0};
    Sequence<std::string> debugBits     {Name("debugBits"     ), Comment("debug bits"                )};
    Atom<std::string>     outputDir     {Name("outputDir"     ), Comment("output directory"          )};
    Atom<int>             saveWaveforms {Name("saveWaveforms" ), Comment("save StrawDigiADCWaveforms"),0};
    Atom<int>             makeSD        {Name("makeSD"        ), Comment("make straw digi branch"     ),1};
    Atom<int>             makeSH        {Name("makeSH"        ), Comment("make straw hit branch"      ),1};
    Atom<int>             makeTC        {Name("makeTC"        ), Comment("make time cluster branch"   ),1};
    Atom<int>             ewLength      {Name("ewLength"      ), Comment("event window length, in units of 25 ns"),1000};
    Atom<int>             nSamplesBL    {Name("nSamplesBL"    ), Comment("n(samples) to determine the BL"),6};
    Atom<float>           minPulseHeight{Name("minPulseHeight"), Comment("min height of the first non-BL sample"),5};
  };

  // --- C'tor/d'tor:
  explicit MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config);
  virtual ~MakeDigiNtuple() {}

  int      getData(const art::Event& ArtEvent);
  
  void     print_(const std::string&  Message, const std::source_location& location = std::source_location::current());

  int      process_adc_waveform(float* Wf, WfParam_t* Wp);

  int      fillSD();
  int      fillSH();
  int      fillTC();
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
  art::InputTag            _shCollTag;          // straw hit collection tag
  art::InputTag            _tcCollTag;          // time cluster collection tag
  std::string              _outputDir;
  int                      _saveWaveforms;
  int                      _makeSD;
  int                      _makeSH;
  int                      _makeTC;
  int                      _ewLength;           // it is up to the user to make sure it is set correctly
  int                      _nSamplesBL;
  int                      _minPulseHeight;
    
  
  int                      _n_adc_samples;
  double                   _tdc_bin;            // TDC bin, in us
  double                   _tdc_bin_ns;         // TDC bin, ns

  DaqEvent*                _event;
  const art::Event*        _art_event;

  int                     _nstrawdigis;
  int                     _nstrawhits;
  int                     _ntimeclusters;
  
  int                     _ncalodigis;
  int                     _ncrvdigis;
  int                     _nstmdigis;
                          
  TFile*                  _file;
  TTree*                  _tree;
  TBranch*                _branch;

  int                     _hist_booked;
  const TrkPanelMap_t*    _panel_map[200][6];   // indexing - offline: [plane][panel]

  const mu2e::StrawDigiCollection*             _sdc;
  const mu2e::StrawDigiADCWaveformCollection*  _sdawfc;
  const mu2e::StrawHitCollection*              _shc;
  const mu2e::TimeClusterCollection*           _tcc;
  
}; // MakeDigiNtuple

// ======================================================================

mu2e::MakeDigiNtuple::MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _debugMode    (config().debugMode    ()),
    _debugBits    (config().debugBits    ()),
    _sdCollTag    (config().sdCollTag    ()),
    _shCollTag    (config().shCollTag    ()),
    _tcCollTag    (config().tcCollTag    ()),
    _outputDir    (config().outputDir    ()),
    _saveWaveforms(config().saveWaveforms()),
    _makeSD       (config().makeSD       ()),
    _makeSH       (config().makeSH       ()),
    _makeTC       (config().makeTC       ()),
    _ewLength     (config().ewLength     ()),
    _nSamplesBL   (config().nSamplesBL   ()),
    _minPulseHeight(config().minPulseHeight   ()),
    _art_event    (nullptr)
{
  _n_adc_samples = -1;

  _tdc_bin             = (5/256.*1e-3);       // TDC bin width (Richie), in us
  _tdc_bin_ns          = _tdc_bin*1e3;        // convert to ns
  _hist_booked         = 0;
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

std::vector<std::string> splitString(const std::string& str, const std::string& delimiter) {
    std::vector<std::string> result;
    std::regex re(delimiter);
    std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
    std::sregex_token_iterator end;
    while (it != end) {
        result.push_back(*it++);
    }
    return result;
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::print_(const std::string& Message, const std::source_location& location) {
  if (_art_event) {
    std::cout << std::format(" event:{}:{}:{}",
                             _art_event->run(),_art_event->subRun(),_art_event->event());
  }


  std::vector<std::string> ss = splitString(location.file_name(),"/");
  // int sz = ss.size();
  
  std::cout << " " << ss.back() << ":" << location.line()
    //            << location.function_name()
            << ": " << Message;
}

//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::beginRun(const art::Run& ArtRun) {

  if (_hist_booked == 0) {
                                        // make sure we do it only once
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
    _hist_booked = 1;
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

  _nstrawdigis   = 0;
  _nstrawhits    = 0;
  _ntimeclusters = 0;
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
    print_(std::format("ERROR: StrawDigiCollection:{:s} is not available. Bail out\n",
                       _sdCollTag.encode().data()));
    return -1;
  }

  ok =  ArtEvent.getByLabel(_sdCollTag,sdawfch);
  if (ok) { 
    _sdawfc = sdawfch.product();
  }
  else {
    print_(std::format("WARNING: StrawDigiADCWaveformCollection:{:s} is not available. Bail out\n",
                       _sdCollTag.encode().data()));
    return -1;
  }

  if (_makeTC != 0) {
    art::Handle<mu2e::TimeClusterCollection>             tcch;
    ok =  ArtEvent.getByLabel(_tcCollTag,tcch);
    if (ok) { 
      _tcc           = tcch.product();
      _ntimeclusters = _tcc->size();
    }
    else {
      print_(std::format("WARNING: TimeClusterCollection:{:s} is not available. Bail out\n",
                       _tcCollTag.encode().data()));
      return -1;
    }
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
int mu2e::MakeDigiNtuple::process_adc_waveform(float* Wf, WfParam_t* Wp) {
//-----------------------------------------------------------------------------
// waveform processing
// 1. determine the baseline
//-----------------------------------------------------------------------------
  Wp->bl = 0;
  for (int i=0; i<_nSamplesBL; ++i) {
    Wp->bl += Wf[i];
  }
  Wp->bl = Wp->bl/_nSamplesBL;
//-----------------------------------------------------------------------------
// 2. subtract the baseline and calculate the charge
//-----------------------------------------------------------------------------
  for (int i=0; i<_n_adc_samples; i++) {
    Wf[i] = Wf[i]-Wp->bl;
  }

  int   tail  = 0;
  Wp->fs = -1;
  Wp->q  = 0;
  Wp->qt = 0;
  Wp->ph = -1;
  for (int i=_nSamplesBL; i<_n_adc_samples; ++i) {
    if (Wf[i] > _minPulseHeight) {
      if (tail == 0) {
                                        // first sample above the threshold
        if (Wp->fs < 0) Wp->fs = i;

        Wp->q += Wf[i];
        if (Wf[i] > Wp->ph) {
          Wp->ph = Wf[i];
        }
      }
    }
    else if (Wf[i] < 0) {
      if (Wp->ph > 0) {
        tail  = 1;
      }
      if (tail == 1) Wp->qt -= Wf[i];
    }
  }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
  // if (Wp->q < 100) {
  //   TLOG(TLVL_DEBUG+1) << "event=" << _art_event->run() << ":"
  //                      << _art_event->subRun() << ":" << _art_event->event() 
  //                      << " Q=" << Wp->q;
  // }
  return 0;
}
  
//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSD() {

  _event->sd->Delete();
  
  for (int i=0; i<_nstrawdigis; i++) {
    const mu2e::StrawDigi*            sd    = &_sdc->at(i);
    const mu2e::StrawDigiADCWaveform* sdawf = &_sdawfc->at(i);
    int ns = sdawf->samples().size();
                                        // one-time initializatiion
    if (_n_adc_samples == -1) _n_adc_samples = ns;

    DaqStrawDigi* nt_sd = new ((*_event->sd)[i]) DaqStrawDigi(ns);
    nt_sd->sid          = sd->strawId().asUint16();
 
    int pln = sd->strawId().plane();
    int pnl = sd->strawId().panel();
    const TrkPanelMap_t* tpm = _panel_map[pln][pnl];

    nt_sd->mnid         = tpm->mnid;
    
    nt_sd->tdc0         = sd->TDC(mu2e::StrawEnd::cal);
    nt_sd->tdc1         = sd->TDC(mu2e::StrawEnd::hv );
    nt_sd->tot0         = sd->TOT(mu2e::StrawEnd::cal);
    nt_sd->tot1         = sd->TOT(mu2e::StrawEnd::hv );
    nt_sd->pmp          = sd->PMP();
    nt_sd->flag         = *((uint8_t*) &sd->digiFlag());
//-----------------------------------------------------------------------------
// store the waveform
//-----------------------------------------------------------------------------
    for (int is=0; is<ns; is++) {
      nt_sd->adc[is] = sdawf->samples()[is];
    }
//-----------------------------------------------------------------------------
// process the waveform and store teh waveform parameters
//-----------------------------------------------------------------------------
    float wf[100];
    for (int is=0; is<ns; is++) {
      wf[is] = sdawf->samples()[is];
    }

    WfParam_t wp;
    process_adc_waveform(wf,&wp);
    
    nt_sd->fs = wp.fs;
    nt_sd->bl = wp.bl;
    nt_sd->ph = wp.ph;

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

  for (int i=0; i<2; i++) {
    for (int k=0; k<6; k++) {
      _event->nsh[i][k] = 0;
    }
  }

  _event->sh->Clear();

  for (int i=0; i<_nstrawhits; i++) {
    const mu2e::StrawHit* sh = &_shc->at(i);
    int pln = sh->strawId().plane();
    int pnl = sh->strawId().panel();
    const TrkPanelMap_t* tpm = _panel_map[pln][pnl];

    int pcie_addr = tpm->dtc % 2; // convention
    _event->nsh[pcie_addr][tpm->link] += 1;
   
    DaqStrawHit* nt_sh = new ((*_event->sh)[i]) DaqStrawHit();
    nt_sh->sid         = sh->strawId().asUint16();
    nt_sh->zface       = tpm->zface;
    nt_sh->mnid        = tpm->mnid;
    nt_sh->time        = sh->time(mu2e::StrawEnd::cal);
    nt_sh->dt          = sh->dt();      // cal-hv
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

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillTC() {

  _event->tc->Clear();

  for (int i=0; i<_ntimeclusters; i++) {
    const mu2e::TimeCluster* tc = &_tcc->at(i);
   
    DaqTimeCluster* nt_tc = new ((*_event->tc)[i]) DaqTimeCluster();

    nt_tc->nsh  = tc->nStrawHits();
    nt_tc->nch  = tc->nhits();
    nt_tc->t0   = tc->t0().t0();
//-----------------------------------------------------------------------------
    // to calculate tmin and tmax need to loop over the hits, not now
    if (_debugMode  > 0) {
      if (_debugBit[2] != 0) {
        printf("%8i %5i %5i %12.2f %12.2f %12.2f\n",
               _event->evn,
               nt_tc->nsh,
               nt_tc->nch,
               nt_tc->t0,
               nt_tc->tmin,
               nt_tc->tmax);
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
  _event->run    = ArtEvent.run();
  _event->srn    = ArtEvent.subRun();
  _event->evn    = ArtEvent.event();
                                        // defined in getData()
  _event->nsd    = _nstrawdigis;
  _event->nshtot = _nstrawhits;
  _event->ntc    = _ntimeclusters;

  //   DaqEvent::fgSd.clear();
  if (_debugMode > 0) {
    print_(std::format("_nstrawdigis:{}",_nstrawdigis));
  }

  if (_makeSD) fillSD();
  if (_makeSH) fillSH();
  if (_makeTC) fillTC();

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
