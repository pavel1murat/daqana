// ======================================================================
// clang-format off
// MakeDigiNtuple:  assume that digis have been produced, make a hit ntuple
//                  for the cross-subsystem timing studies
// tracker       :  assume that the waveforms are made
// debug bits    :  1: for all hits, print hit times assuming contiguous timing
//                  2: (fillTC) parameters of the time cluster
//                  3: TrkSegment::fgDebugMode
//                  4: SegmentFit::fgDebugMode
//                  5: print HepTransform's for all panels
//                  6: print reconstructed segments in the end
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
// #include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"

#include "art/Framework/Principal/Handle.h"
// #include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"

// #include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
// #include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"


#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "daqana/obj/DaqEvent.hh"
#include "daqana/mod/WfParam_t.hh"
#include "daqana/obj/TrkSegment.hh"
#include "daqana/obj/SegmentFit.hh"

#include <ostream>
#include <regex>
// #include <ranges>
#include <string>
#include <vector>
#include <algorithm>

#include <iostream>
// #include <memory>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGeoMatrix.h"

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
    Atom<art::InputTag>   ksCollTag     {Name("ksCollTag"     ), Comment("KS coll tag"               ),""};
    Atom<int>             debugMode     {Name("debugMode"     ), Comment("debug mode"                )};
    Sequence<std::string> debugBits     {Name("debugBits"     ), Comment("debug bits"                )};
    Atom<std::string>     outputDir     {Name("outputDir"     ), Comment("output directory"          )};
    Atom<int>             saveWaveforms {Name("saveWaveforms" ), Comment("save StrawDigiADCWaveforms")};
    Atom<int>             makeSD        {Name("makeSD"        ), Comment("make straw digi branch"     ),1};
    Atom<int>             makeSH        {Name("makeSH"        ), Comment("make straw hit branch"      ),1};
    Atom<int>             makeCH        {Name("makeCH"        ), Comment("make combo hit branch"      ),1};
    Atom<int>             makeTC        {Name("makeTC"        ), Comment("make time cluster branch"   ),1};
    Atom<int>             makeSeg       {Name("makeSeg"       ), Comment("make segment branch"        ),1};
    Atom<int>             makeTrk       {Name("makeTrk"       ), Comment("make track branch"          ),1};
    Atom<int>             ewLength      {Name("ewLength"      ), Comment("event window length, in units of 25 ns"),1000};
    Atom<int>             nSamplesBL    {Name("nSamplesBL"    ), Comment("n(samples) to determine the BL"),6};
    Atom<float>           minPulseHeight{Name("minPulseHeight"), Comment("min height of the first non-BL sample"),5};
    Atom<int>             minNSegments  {Name("minNSegments"  ), Comment("min N(segments)"                     )};
    Atom<float>           vDrift        {Name("vDrift"        ), Comment("vDrift, um/ns")               };
    Atom<float>           tOffset       {Name("tOffset"       ), Comment("T0 offset, ns")               };
  };

  // --- C'tor/d'tor:
  explicit MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config);
  virtual ~MakeDigiNtuple() {}

  int      getData(const art::Event& ArtEvent);
  
  void     print_(const std::string&  Message, const std::source_location& location = std::source_location::current());

  int      process_adc_waveform(float* Wf, WfParam_t* Wp);

  int      calculateMissingTrkParameters();
  
  int      makeSegments();

  int      fillSD ();
  int      fillSH ();
  int      fillCH ();
  int      fillTC ();
  int      fillSeg();
  int      fillSegSh();
  int      fillTrk();
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDAnalyzer
//-----------------------------------------------------------------------------
  virtual  void beginRun(const art::Run&   r);
  virtual  void analyze (const art::Event& e);
  virtual  void endRun  (const art::Run&   r);
  virtual  void beginJob();
  virtual  void endJob  ();

  int                      _debugMode;
  std::vector<std::string> _debugBits;
  int                      _debugBit[100];
  art::InputTag            _sdCollTag;          // straw digi collection tag
  art::InputTag            _shCollTag;          // straw hit collection tag
  art::InputTag            _tcCollTag;          // time cluster collection tag
  art::InputTag            _ksCollTag;          // kalseed collection tag
  std::string              _outputDir;
  int                      _saveWaveforms;
  int                      _makeSD;
  int                      _makeSH;
  int                      _makeCH;
  int                      _makeTC;
  int                      _makeSeg;
  int                      _makeTrk;
  int                      _ewLength;           // it is up to the user to make sure it is set correctly
  int                      _nSamplesBL;
  int                      _minPulseHeight;
  int                      _minNSegments;
  float                    _vDrift;
  float                    _tOffset;
    
  
  int                      _n_adc_samples;
  double                   _tdc_bin;            // TDC bin, in us
  double                   _tdc_bin_ns;         // TDC bin, ns

  DaqEvent*                _event;
  const art::Event*        _art_event;
  int                      _last_run;

  int                     _nstrawdigis;
  int                     _nstrawhits;
  int                     _ncombohits;
  int                     _ntimeclusters;
  int                     _nsegments;
  int                     _ntracks;
  
  int                     _ncalodigis;
  int                     _ncrvdigis;
  int                     _nstmdigis;
                          
  TFile*                  _file;
  TTree*                  _tree;
  TBranch*                _branch;

  int                     _hist_booked;
   
  const mu2e::StrawDigiCollection*             _sdc;
  const mu2e::StrawDigiADCWaveformCollection*  _sdawfc;
  const mu2e::StrawHitCollection*              _shc;
  const mu2e::ComboHitCollection*              _chc;
  const mu2e::TimeClusterCollection*           _tcc;
  const mu2e::KalSeedCollection*               _ksc;

  ProditionsHandle<Tracker>                    _alignedTracker_h;
  const mu2e::Tracker*                         _tracker;
  
  ProditionsHandle<TrackerPanelMap>            _tpm_h;
  const TrackerPanelMap*                       _trkPanelMap;
  
  TrkSegment                                   _tseg[36][6]; // for now, assume just one segment per panel
  std::vector<TrkSegment*>                     _ptseg;
  int                                          _nseg;
  
}; // MakeDigiNtuple

// ======================================================================

mu2e::MakeDigiNtuple::MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _debugMode     (config().debugMode     ()),
    _debugBits     (config().debugBits     ()),
    _sdCollTag     (config().sdCollTag     ()),
    _shCollTag     (config().shCollTag     ()),
    _tcCollTag     (config().tcCollTag     ()),
    _ksCollTag     (config().ksCollTag     ()),
    _outputDir     (config().outputDir     ()),
    _saveWaveforms (config().saveWaveforms ()),
    _makeSD        (config().makeSD        ()),
    _makeSH        (config().makeSH        ()),
    _makeCH        (config().makeCH        ()),
    _makeTC        (config().makeTC        ()),
    _makeSeg       (config().makeSeg       ()),
    _makeTrk       (config().makeTrk       ()),
    _ewLength      (config().ewLength      ()),
    _nSamplesBL    (config().nSamplesBL    ()),
    _minPulseHeight(config().minPulseHeight()),
    _minNSegments  (config().minNSegments  ()),
    _vDrift        (config().vDrift        ()),
    _tOffset       (config().tOffset       ()),
    _art_event     (nullptr)
{
  _n_adc_samples = -1;

  _tdc_bin             = (5/256.*1e-3);       // TDC bin width (Richie), in us
  _tdc_bin_ns          = _tdc_bin*1e3;        // convert to ns
  _hist_booked         = 0;
  _last_run            = -1;
  _nseg                = 0;
//-----------------------------------------------------------------------------
// parse debug bits
//-----------------------------------------------------------------------------
  for (int i=0; i<100; i++) _debugBit[i] = 0;

  const char* key;
                                        // a flag is an integer!
  int nbits = _debugBits.size();
  for (int i=0; i<nbits; i++) {
    int index(0), value(0);
    key               = _debugBits[i].data();
    sscanf(key,"bit%i:%i",&index,&value);
    _debugBit[index]  = value;
    
    print_(std::format("...{}: bit={:4d} is set to {}\n",__func__,index,_debugBit[index]));
  }

  SegmentHit::SetVDrift (_vDrift );
  SegmentHit::SetTOffset(_tOffset);
  
  TrkSegment::fgDebugMode = _debugBit[3];
  SegmentFit::fgDebugMode = _debugBit[4];
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
// Message should be \n terminated , if needed
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
//-----------------------------------------------------------------------------
// for each panel, build a transformation matrix, do it once
//-----------------------------------------------------------------------------
  // mu2e::GeomHandle<mu2e::Tracker> handle;
  // _tracker = handle.get();
  art::EventID eid(ArtRun.run(),1,1);
  _alignedTracker_h = mu2e::ProditionsHandle<mu2e::Tracker>();
  _tracker = _alignedTracker_h.getPtr(eid).get();
  
  for (int ipln=0; ipln<36; ipln++) {
    for (int ipnl=0; ipnl<6; ipnl++) {
      TrkSegment* ts = &_tseg[ipln][ipnl];
      ts->fPlane = ipln;
      ts->fPanel = ipnl;
//-----------------------------------------------------------------------------
// 1.build transformation matrix
//-----------------------------------------------------------------------------
      const mu2e::Plane* pln = &_tracker->getPlane(ipln);
      const mu2e::Panel* pnl = &pln->getPanel(ipnl);
      ts->fTrkPanel = (mu2e::Panel*) pnl;

      if ((_debugMode != 0) and (_debugBit[5])) {
        print_(std::format("-- HepTransform for plane:{:2}:{}\n",ipln,ipnl));
        std::cout << ts->fTrkPanel->dsToPanel() << std::endl;
      }
    }
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
  art::Handle<mu2e::ComboHitCollection>             chch;

  _nstrawdigis   = 0;
  _nstrawhits    = 0;
  _ncombohits    = 0;
  _ntimeclusters = 0;
  // _ncalodigis  = 0;
  // _ncrvdigis   = 0;
  // _nstmdigis   = 0;

  _sdc         = nullptr;
  _sdawfc      = nullptr;
  _shc         = nullptr;
  _chc         = nullptr;
  
  bool ok = ArtEvent.getByLabel(_shCollTag,shch);
  if (ok) { 
    _shc         = shch.product();
    _nstrawhits = _shc->size();
  }
  else {
    print_(std::format("ERROR: StrawHitCollection:{:s} is not available. Bail out\n",_shCollTag.encode().data()));
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
    art::Handle<mu2e::TimeClusterCollection>          tcch;
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
//-----------------------------------------------------------------------------
// assume that chCollTag == _shCollTag
//-----------------------------------------------------------------------------
    art::Handle<mu2e::ComboHitCollection>             chch;
    ok =  ArtEvent.getByLabel(_shCollTag,chch);
    if (ok) { 
      _chc           = chch.product();
      _ncombohits    = _chc->size();
    }
    else {
      print_(std::format("WARNING: ComboHitCollection:{:s} is not available. Bail out\n",
                       _shCollTag.encode().data()));
      return -1;
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  if (_makeTrk != 0) {
    art::Handle<mu2e::KalSeedCollection>  ksch;
    ok =  ArtEvent.getByLabel(_ksCollTag,ksch);
    if (ok) { 
      _ksc     = ksch.product();
      _ntracks = _ksc->size();
    }
    else {
      print_(std::format("WARNING: KalSeedCollection:{:s} is not available. Bail out\n",
                         _ksCollTag.encode().data()));
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
    int ich = sd->strawId().straw();
    //    const TrkPanelMap_t* tpm = _panel_map[pln][pnl];
    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);

    // works for one station, not more...

    int pcie_addr = tpm->dtc() % 2;                   // convention
    int link      = tpm->link();
    _event->nsd[pcie_addr][link][ich] += 1;

    nt_sd->mnid         = tpm->mnid();
    
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

  if (_debugMode > 0) {
    if (_debugBit[1] != 0) {
      printf("evn    sid  pln  pnl mnid    time    dt   tot0 tot1   edep\n");
      printf("---------------------------------------------------------\n");
    }
  }
  for (int i=0; i<_nstrawhits; i++) {
    const mu2e::StrawHit* sh = &_shc->at(i);
    int pln = sh->strawId().plane();
    int pnl = sh->strawId().panel();
    //    const TrkPanelMap_t* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);
    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);

    int pcie_addr = tpm->dtc() % 2;                      // convention
    _event->nsh[pcie_addr][tpm->link()] += 1;
   
    DaqStrawHit* nt_sh = new ((*_event->sh)[i]) DaqStrawHit();
    nt_sh->sid         = sh->strawId().asUint16();
    nt_sh->zface       = tpm->zface();
    nt_sh->mnid        = tpm->mnid();
    nt_sh->time        = sh->time(mu2e::StrawEnd::cal);
    nt_sh->dt          = sh->dt();      // cal-hv
    nt_sh->tot0        = sh->TOT(mu2e::StrawEnd::cal);
    nt_sh->tot1        = sh->TOT(mu2e::StrawEnd::hv );
    nt_sh->edep        = sh->energyDep();
    if (sh->energyDep() > _event->maxEdep) _event->maxEdep = sh->energyDep();

    if (_debugMode  > 0) {
      if (_debugBit[1] != 0) {
        printf("%8i %5i %3i %3i %3i %12.4f %12.4f %6.1f %6.1f %7.4f\n",
               _event->evn,
               (int) nt_sh->sid, pln, pnl, nt_sh->mnid,
               nt_sh->time, nt_sh->dt,
               nt_sh->tot0, nt_sh->tot1,
               nt_sh->edep);
      }
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// comboo hits
//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCH() {

  if (_debugMode > 0) {
    if (_debugBit[1] != 0) {
      printf("evn    sid  pln  pnl mnid    time    dt   tot0 tot1   edep\n");
      printf("---------------------------------------------------------\n");
    }
  }
  for (int i=0; i<_ncombohits; i++) {
    const mu2e::ComboHit* ch = &_chc->at(i);
    int pln = ch->strawId().plane();
    int pnl = ch->strawId().panel();
    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);

    DaqComboHit* nt_ch = new ((*_event->ch)[i]) DaqComboHit();
    nt_ch->sid         = ch->strawId().asUint16();
    nt_ch->nsh         = ch->nStrawHits();
    nt_ch->zface       = tpm->zface();
    nt_ch->mnid        = tpm->mnid();
    nt_ch->time        = ch->correctedTime();
    nt_ch->dtime       = ch->driftTime();
    nt_ch->x           = ch->pos().x();
    nt_ch->y           = ch->pos().y();
    nt_ch->z           = ch->pos().z();
    nt_ch->ux          = ch->uDir().x();
    nt_ch->uy          = ch->uDir().y();
    nt_ch->ures        = ch->uRes();
    nt_ch->vres        = ch->vRes();

    if (_debugMode  > 0) {
      if (_debugBit[1] != 0) {
        // printf("%8i %5i %3i %3i %3i %12.4f %12.4f %6.1f %6.1f %7.4f\n",
        //        _event->evn,
        //        (int) nt_sh->sid, pln, pnl, nt_sh->mnid,
        //        nt_sh->time, nt_sh->dt,
        //        nt_sh->tot0, nt_sh->tot1,
        //        nt_sh->edep);
      }
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillTC() {

  for (int itc=0; itc<_ntimeclusters; itc++) {
    const mu2e::TimeCluster* tc = &_tcc->at(itc);
   
    DaqTimeCluster* nt_tc = new ((*_event->tc)[itc]) DaqTimeCluster();

    nt_tc->nsh     = tc->nStrawHits();
    nt_tc->nch     = tc->nhits();
    nt_tc->t0      = tc->t0().t0();
    
    LsqSums2 szy;

    for (int ih=0; ih<nt_tc->nch; ih++) {
      StrawHitIndex hit_index   = tc->hits().at(ih);
      const mu2e::ComboHit* hit = &_chc->at(hit_index);

      double y    = hit->pos().Y();
      double z    = hit->pos().Z();
      double sigy = hit->uRes()*hit->uDir().X();
      double w    = 1./sigy/sigy;
      szy.addPoint(z,y,w);
      
      int plane = hit->strawId().plane();
      int pln2  = plane % 2;
      int panel = hit->strawId().panel();
      //      const TrkPanelMap_t *tpm = _panel_map[plane][panel];
      const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(plane,panel);
      int zface = tpm->zface();

      int loc = 6*pln2+panel;
      
      if (nt_tc->_nh_panel[loc] == 0) nt_tc->npanels++;
      nt_tc->_nh_panel[loc]++;
      nt_tc->_time_panel[loc] += hit->correctedTime();
      nt_tc->_edep_panel[loc] += hit->energyDep();

      if (hit->energyDep() > nt_tc->edep_max) nt_tc->edep_max = hit->energyDep();

      if (nt_tc->_nhp[pln2] == 0) nt_tc->nplanes++;
      nt_tc->_nhp[pln2]++;
      nt_tc->_timep [pln2] += hit->correctedTime();

      if (nt_tc->_nhf[zface] == 0) {
        nt_tc->nfaces++;
        nt_tc->_mnid[zface] = tpm->mnid();
      }
      nt_tc->_nhf[zface]++;
      nt_tc->_timef [zface] += hit->correctedTime();

      if (hit->correctedTime() < nt_tc->tmin) nt_tc->tmin = hit->correctedTime();
      if (hit->correctedTime() > nt_tc->tmax) nt_tc->tmax = hit->correctedTime();
    }

    nt_tc->y0     = szy.y0();
    nt_tc->dydz   = szy.dydx();
    nt_tc->chi2yz = szy.chi2Dof();

    for (int ip=0; ip<2; ip++) {
      if (nt_tc->_nhp[ip] > 0) nt_tc->_timep[ip] = nt_tc->_timep[ip]/nt_tc->_nhp[ip];
    }

    for (int i=0; i<4; i++) {
      if (nt_tc->_nhf[i] > 0) nt_tc->_timef[i] = nt_tc->_timef[i]/nt_tc->_nhf[i];
    }
//-----------------------------------------------------------------------------
// average time and charge of the hits in a given panel
//-----------------------------------------------------------------------------
    for (int i=0; i<12; i++) {
      if (nt_tc->_nh_panel[i] > 0) {
        nt_tc->_time_panel[i] = nt_tc->_time_panel[i]/nt_tc->_nh_panel[i];
        nt_tc->_edep_panel[i] = nt_tc->_edep_panel[i]/nt_tc->_nh_panel[i];
      }

      if (nt_tc->_nh_panel[i] > nt_tc->max_nh_panel) {
        nt_tc->max_nh_panel = nt_tc->_nh_panel[i];
      }
    }
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

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::calculateMissingTrkParameters() {
  int rc(0);

  if (_debugMode != 0) std::cout << __func__ << " START" << std::endl;
  if (_debugMode != 0) std::cout << __func__ << " END"   << std::endl;
  return rc;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillTrk() {

  calculateMissingTrkParameters();

  for (int itrk=0; itrk<_ntracks; itrk++) {
    const mu2e::KalSeed* ks = &_ksc->at(itrk);
   
    DaqTrack* nt_trk = new ((*_event->trk)[itrk]) DaqTrack();

    nt_trk->nhits    = ks->nHits();
    nt_trk->chi2     = ks->chisquared();
    nt_trk->t0       = ks->t0().t0();
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::makeSegments() {
  int rc(0);

  if (_debugMode != 0) print_(std::format("{}  START: ntimeclusters:{}\n", __func__,_ntimeclusters));
//-----------------------------------------------------------------------------
// cleanup from the previous event, initially set _nseg to 0
//-----------------------------------------------------------------------------
  for (int i=0; i<_nseg; i++) {
    _ptseg[i]->Clear();
  }
  _ptseg.clear();
  _nseg = 0;
//-----------------------------------------------------------------------------
// for each time cluster make segments - see how it goes
// will need to handle the segment hits as well
//-----------------------------------------------------------------------------
  for (int itc=0; itc<_ntimeclusters; itc++) {
    const mu2e::TimeCluster* tc = &_tcc->at(itc);
    int nsh = tc->nStrawHits();
    if (_debugMode) {
      if (_debugBit[3] != 0) {
        print_(std::format(" -- PM: time cluster:{:2d} segment with {:2d} hits\n",itc, nsh));
      }
    }
    for (int ih=0; ih<nsh; ih++) {
      StrawHitIndex hit_index   = tc->hits().at(ih);
      const mu2e::ComboHit* ch = &_chc->at(hit_index);
      mu2e::StrawId sid = ch->strawId();
      int panel = sid.panel();
      int plane = sid.plane();
      if (_debugBit[3] != 0) {
        print_(std::format(" {}: hit number:{:3d} plane:{} panel:{} straw:{:2d}\n",
                           __func__,ih,plane,panel,sid.straw()));
      }
//-----------------------------------------------------------------------------
// cosmic track: assume one segment per panel, in principle, there could be more than one
//-----------------------------------------------------------------------------
      int pln2 = plane % 2;
      TrkSegment* ts = &_tseg[pln2][panel];
      if (ts->nHits() == 0) {
        _ptseg.push_back(ts);
        _nseg    += 1;
      }
                                        // and add the hit to the list
      SegmentHit sgh(ch);
      ts->fListOfHits.emplace_back(sgh);
    }
  }
//-----------------------------------------------------------------------------
// cleanup: loop over segments one more time and try to identify 'extra' hits
//-----------------------------------------------------------------------------
  struct SubSegment {
    int first_hit;                      // indes in the segment hit list
    int last_hit;
    int nhits() { return last_hit-first_hit+1; }
  };

  for (int i=0; i<_nseg; i++) {
    TrkSegment* ts = _ptseg[i];
//-----------------------------------------------------------------------------
// sort segment hits in the ascending wire order
// make sure that the standalone ROOT code also does that
//-----------------------------------------------------------------------------
    std::sort(ts->fListOfHits.begin(),ts->fListOfHits.end(),
              [] (const SegmentHit& a, const SegmentHit& b) {
                return a.ComboHit()->strawId().straw() < b.ComboHit()->strawId().straw();
              });

    int nhits = ts->nHits();
    if (_debugMode != 0) print_(std::format("{}  iseg:{} nhits:{}\n",__func__,i,nhits));
    if (nhits < 4) {
      ts->fMask |= 0x1 ; // not enough hits
                                                                  continue;
    }

    std::vector<SubSegment> list_of_subsegments;

    int last_layer  = -1;
    int last_straw  = -1;
    int first_hit   = -1;
    int last_hit    = -1;
    int best        = -1;
    int nmax        = -1;
    for (int i=0; i<nhits; ++i) {
      const mu2e::ComboHit* ch = ts->Hit(i)->ComboHit();
      int   layer = ch->strawId().getLayer();
      int   straw = ch->strawId().getStraw();
      //      float segment_length = (straw-first_straw)/2.;
      if (last_layer == -1) {
        last_layer = layer;
        last_straw = straw;
        first_hit  = i;
      }

//-----------------------------------------------------------------------------
// check the gap size may need different constants
//-----------------------------------------------------------------------------
      int gap = (straw-last_straw-2)/2;
      if (gap > 4) {
                                        // too large of a gap, make a sub-segment
        SubSegment sbs(first_hit,last_hit);
        list_of_subsegments.push_back(sbs);
        if (sbs.nhits() > nmax) {
          nmax = sbs.nhits();
          best = list_of_subsegments.size()-1;
        }
        first_hit = i;
      }
      if (layer != last_layer) {
        last_layer = layer;
      }
      last_hit   = i;
      last_straw = straw;
    }
//-----------------------------------------------------------------------------
// last subsegment
//-----------------------------------------------------------------------------
    if (first_hit != -1) {
      SubSegment sbs(first_hit,last_hit);
      list_of_subsegments.emplace_back(sbs);
      if (sbs.nhits() > nmax) {
        nmax = sbs.nhits();
        best = list_of_subsegments.size()-1;
      }
    }
//-----------------------------------------------------------------------------
// out of subsegments leave only the one with the largest number of hits
//-----------------------------------------------------------------------------
    int nsbs = list_of_subsegments.size();
    for (int isbs=0; isbs<nsbs; ++isbs) {
      if (isbs != best) {
                                        // wrong subsegment, mask its hits off
        SubSegment& sbs = list_of_subsegments[isbs];
        for (int ih=sbs.first_hit; ih<sbs.last_hit+1; ++ih) {
          ts->Hit(ih)->fMask |= TrkSegment::kSubsegmentBit;
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// loop over segments with 3 hits or more, fi each segment and compare
// the segment position ond slope with the track parameters transpated into the local
// coordinate system of the panel (or in the global ?)
// need add more segment parameters : n(layers) = ? - require N(transitions) = 1
// - number of good hits in each layer
// - number of 'holes' on each side - tomorrow? 
//-----------------------------------------------------------------------------

  TrkSegment::Par_t par;
  int niter(4);
  for (int i=0; i<_nseg; i++) {
    TrkSegment* ts = _ptseg[i];
    if (_debugMode != 0) print_(std::format("{}  iseg:{} nhits:{}\n",__func__,i,ts->nHits()));
    if (ts->nHits() < 4) {
      ts->fMask |= 0x1 ; // not enough hits
      continue;
    }
//-----------------------------------------------------------------------------
// list of hits already created, but SegmentHits need to be initialized from ComboHits
// this looks ugly, but OK for the purpose
//----------------------------------------------------------------------------- 
    int rc = ts->InitHits(nullptr);

    int converged(0);
    if (rc == 0) {
      SegmentFit sfitter(ts);
//-----------------------------------------------------------------------------
// use tangent line and the first and the last hits
// perform 4 fits, find the best
//-----------------------------------------------------------------------------
      sfitter.DefineDriftDirections();
//-----------------------------------------------------------------------------
// perform fit using all points and starting from ts.fPar4
//-----------------------------------------------------------------------------
      converged = sfitter.Fit(niter,0,nullptr,&par);
    }

    if (_debugMode) {
      std::cout << "rc:" << rc << " converged:" << converged;
    }
  }
//-----------------------------------------------------------------------------
// at this point, all track hits should be assigned to segments
// debug printout
//-----------------------------------------------------------------------------
  if (_debugMode and (_debugBit[6] != 0)) {
    std::cout << "nseg:" << _nseg << std::endl;

    for (int i=0; i<_nseg; i++) {
      TrkSegment* ts = _ptseg[i];
      ts->print();
    }
  }

  if (_debugMode != 0) std::cout << __func__ << ":END rc:" << rc << std::endl;
  return 0;
}


//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSegSh() {
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSeg() {

  makeSegments();

  int nsegsh = 0;
  int nseg4  = 0;
  for (int iseg=0; iseg<_nseg; iseg++) {
    TrkSegment* ts = _ptseg[iseg];
    //    if (ts->points.size() < 4) continue;
    nseg4 += 1;
   
    DaqSegment* nt_ts = new ((*_event->seg)[iseg]) DaqSegment();
    mu2e::StrawId sid = ts->Hit(0)->ComboHit()->strawId();

    nt_ts->sid     = sid.asUint16() & 0xff80;            // straw ID of the straw=0 of the panel
    nt_ts->nh      = ts->nHits();
    nt_ts->ntrans  = ts->fNTransitions;
    nt_ts->ngh     = ts->fNGoodHits;
    nt_ts->nghl[0] = ts->fNghl[0];
    nt_ts->nghl[1] = ts->fNghl[1];
    nt_ts->nmhl[0] = ts->fNmhl[0];                       // number of missing hits/layer
    nt_ts->nmhl[1] = ts->fNmhl[1];                       // number of missing hits/layer
    nt_ts->t0      = ts->T0()+ts->fTMean;
    nt_ts->chi2d   = ts->fPar.chi2dof;                   // chi2/dof
    nt_ts->z0      = -1.;
    if (ts->DyDx() != 0) nt_ts->y0 = -ts->Y0()/ts->DyDx();
    else                 nt_ts->y0 = -1.e6;             
    nt_ts->ymean   = ts->fXMean;                         // local on a panel
    nt_ts->dzdy    = ts->DyDx();
    if (_ntracks > 0) {
                                        // transform track parameters into a local coordinate system
                                        // assume one track...
      
      const mu2e::KalSeed*   ks = &_ksc->at(0);

      KinKal::KinematicLine   kline = ks->nearestSegment(ts->Hit(0)->ComboHit()->pos())->kinematicLine();
      ROOT::Math::XYZVector   pos   = kline.pos0();
      ROOT::Math::XYZVector   dir   = kline.direction();
     
      CLHEP::Hep3Vector dirm(dir.x(),dir.y(),dir.z());
      CLHEP::Hep3Vector dirl = ts->fTrkPanel->dsToPanel().rotation()*dirm;

      CLHEP::Hep3Vector posm(pos.x(),pos.y(),pos.z());
      CLHEP::Hep3Vector posl = ts->fTrkPanel->dsToPanel()*posm;

      nt_ts->y0t     = 0;                                     // at z=Z(mid panel) - to be figured
      nt_ts->dzdyt   = dirl[2]/dirl[1];                       // dzdy of the track (local coord system)
    }
    else {
      nt_ts->y0t     = 1.e6;                                     // at z=Z(mid panel) - to be figured
      nt_ts->dzdyt   = 1.e6;
    }
    if (_debugMode) {
      print_(std::format(" iseg:{} dz/dy(seg):{:12.5f} dz/dy(trk):{:12.5f}\n",iseg,nt_ts->dzdy,nt_ts->dzdyt));
    }
//-----------------------------------------------------------------------------
// fill segment straw hit branch
//-----------------------------------------------------------------------------
    int pln = sid.plane();
    int pnl = sid.panel();
    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);
    DaqTrkStrawHit* nt_tsh(nullptr);
    for (int ih=0; ih<nt_ts->nh; ih++) {
      const mu2e::StrawHit* sh = &_shc->at(ih);
      SegmentHit* sgh          = ts->Hit(ih);
      const mu2e::ComboHit* ch = sgh->ComboHit();

      int ihh            = nsegsh+ih;
      nt_tsh = new ((*_event->segsh)[ihh]) DaqTrkStrawHit();

      nt_tsh->sid     = ch->strawId().asUint16();           // hit id = sid | (segment #) << 16 | (track #) << 24
      nt_tsh->zface   = tpm->zface();                       // z-ordered face ... can be deduced from sid...
      nt_tsh->mnid    = tpm->mnid();                        // Minnesota panel ID 
      nt_tsh->time    = sh->time(mu2e::StrawEnd::cal); // 0:CAL
      nt_tsh->dt      = sh->dt();                           // cal-hv
      nt_tsh->tot0    = sh->TOT(mu2e::StrawEnd::cal);
      nt_tsh->tot1    = sh->TOT(mu2e::StrawEnd::hv );
      nt_tsh->edep    = ch->energyDep();

      nt_tsh->rdrift  = ts->R  (sgh);                    // drift distance
      nt_tsh->doca    = ts->Doca(sgh);                      // track-wire distance, signed 
      nt_tsh->dr      = ts->Dr(sgh);                        // (track-wire distance)-Rdrift*drift_sign
//-----------------------------------------------------------------------------
// if drho is positive, the hit drift radius needs to be increased, and the T0 - reduced...
//-----------------------------------------------------------------------------
      nt_tsh->drho    = fabs(nt_tsh->doca)-fabs(nt_tsh->rdrift); // unsigned residual

      nt_tsh->iseg    = iseg;
      nt_tsh->itrk    = -1;
      nt_tsh->ihit    = ih;
    }
    nsegsh += nt_ts->nh;
  }

  _event->nseg   = _nseg;
  _event->nsegsh = nsegsh;

  return 0;
}


//-----------------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::analyze(const art::Event& ArtEvent) {
  //  int const packet_size(16); // in bytes

  _art_event = &ArtEvent;

  if (_last_run != (int) _art_event->run()) {
    _trkPanelMap = &_tpm_h.get(ArtEvent.id());
    _last_run    = _art_event->run();
  }

  if (_debugMode > 0) {
    print_(std::format("-- START event:{}:{}:{}\n",ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event()));
  }

  int rc = getData(ArtEvent);
  if (rc < 0) return;
//-----------------------------------------------------------------------------
// clear event
//-----------------------------------------------------------------------------
  _event->Clear();

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
  _event->nsdtot  = _nstrawdigis;
  _event->nshtot  = _nstrawhits;
  _event->nch     = _ncombohits;
  _event->ntc     = _ntimeclusters;
  _event->ntrk    = _ntracks;
  _event->maxEdep = 0;

  if (_debugMode > 0) {
    print_(std::format("_nstrawdigis:{}\n",_nstrawdigis));
  }

  if (_makeSD ) fillSD ();
  if (_makeSH ) fillSH ();
  if (_makeCH ) fillCH ();
  if (_makeTC ) fillTC ();
  if (_makeSeg) {
    fillSeg  ();
    fillSegSh();
  }
  if (_makeTrk) fillTrk();

  if (_debugMode > 0) {
    print_(std::format("_event->strawdigis->GetEntries():{}\n",_event->nsdtot));
  }

  if (_event->nseg >= _minNSegments) {
    _tree->Fill();
  }

  if (_debugMode > 0) print_("-- END\n");
}



// ======================================================================

DEFINE_ART_MODULE(mu2e::MakeDigiNtuple)

// ======================================================================
