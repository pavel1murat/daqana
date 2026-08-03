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
//                 21: print number of CRV things
//                 31-40: makeSegments printout
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
// #include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"

#include "art/Framework/Principal/Handle.h"
#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_SubEventHeader.h"
#include "artdaq-core-mu2e/Overlays/DTC_Types/DTC_Subsystem.h"

#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
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
#include "TRACE/tracemf.h"
#define TRACE_NAME "MakeDigiNtuple"
#define LOG_STREAM "MakeDigiNtuple"

namespace mu2e {
  class MakeDigiNtuple;
}

using namespace fhicl;

// ======================================================================
class mu2e::MakeDigiNtuple : public art::EDAnalyzer {

public:

  enum {
    e_DEBUG   = 0,
    e_INFO    = 1,
    e_WARNING = 2,
    e_ERROR   = 3,
    e_SEVERE  = 4,
  };

  struct Config {
    
    Atom<art::InputTag>   caldCollTag   {Name("caldCollTag"   ), Comment("calorimeter digi coll tag"  )};
    Atom<art::InputTag>   calhCollTag   {Name("calhCollTag"   ), Comment("calorimeter hits coll tag"  )};
    Atom<art::InputTag>   calcCollTag   {Name("calcCollTag"   ), Comment("calo cluster coll tag"      )};
    Atom<art::InputTag>   crvdCollTag   {Name("crvdCollTag"   ), Comment("CRV digi coll tag"          )};
    Atom<art::InputTag>   crvpCollTag   {Name("crvpCollTag"   ), Comment("CRV reco pulse coll tag"    )};
    Atom<art::InputTag>   crvcCollTag   {Name("crvcCollTag"   ), Comment("CRV coins cluster coll tag" )};
    Atom<art::InputTag>   sdCollTag     {Name("sdCollTag"     ), Comment("straw digi coll tag"        )};
    Atom<art::InputTag>   shCollTag     {Name("shCollTag"     ), Comment("straw hit  coll tag"        )};
    Atom<art::InputTag>   tcCollTag     {Name("tcCollTag"     ), Comment("time cluster coll tag"      )};
    Atom<art::InputTag>   ksCollTag     {Name("ksCollTag"     ), Comment("KS coll tag"                )};  // ,""};
    Atom<int>             debugMode     {Name("debugMode"     ), Comment("debug mode"                 )};
    Sequence<std::string> debugBits     {Name("debugBits"     ), Comment("debug bits"                 )};
    Atom<std::string>     outputDir     {Name("outputDir"     ), Comment("output directory"           )};
    Atom<int>             saveWaveforms {Name("saveWaveforms" ), Comment("save StrawDigiADCWaveforms" )};
    Atom<int>             makeCalD      {Name("makeCalD"      ), Comment("make CAL digis"             )}; // ,1};
    Atom<int>             makeCalH      {Name("makeCalH"      ), Comment("make CAL hits"              )}; // ,1};
    Atom<int>             makeCalC      {Name("makeCalC"      ), Comment("make CAL clusters"          )}; // ,1};
    Atom<int>             makeCrvD      {Name("makeCrvD"      ), Comment("make CRV digis"             )}; // ,1};
    Atom<int>             makeCrvP      {Name("makeCrvP"      ), Comment("make CRV pulses"            )}; // ,1};
    Atom<int>             makeCrvC      {Name("makeCrvC"      ), Comment("make CRV cclusters"         )}; // ,1};
    Atom<int>             makeSD        {Name("makeSD"        ), Comment("make straw digi branch"     )}; // ,1};
    Atom<int>             makeSH        {Name("makeSH"        ), Comment("make straw hit branch"      )}; // ,1};
    Atom<int>             makeCH        {Name("makeCH"        ), Comment("make combo hit branch"      )}; // ,1};
    Atom<int>             makeFragments {Name("makeFragments" ), Comment("make artdaq branch"         )};       // ,1};
    Atom<int>             makeTC        {Name("makeTC"        ), Comment("make time cluster branch"   )};       // ,1};
    Atom<int>             makeSeg       {Name("makeSeg"       ), Comment("make segment branch"        )};        // ,1};
    Atom<int>             makeTrk       {Name("makeTrk"       ), Comment("make track branch"          )};        // ,1};
    Atom<int>             ewLength      {Name("ewLength"      ), Comment("event window length, in units of 25 ns"),1000};
    Atom<int>             nSamplesBL    {Name("nSamplesBL"    ), Comment("n(samples) to determine the BL"),6};
    Atom<float>           minPulseHeight{Name("minPulseHeight"), Comment("min height of the first non-BL sample"),5};
    Atom<float>           minSDPHToSave {Name("minSDPHToSave" ), Comment("min PH of the SD to save")   };
    Atom<int>             minNSegments  {Name("minNSegments"  ), Comment("min N(segments)")            };
    Atom<float>           vDrift        {Name("vDrift"        ), Comment("vDrift, um/ns")              };
    Atom<float>           tOffset       {Name("tOffset"       ), Comment("T0 offset, ns")              };
  };

  // --- C'tor/d'tor:
  explicit MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config);
  virtual ~MakeDigiNtuple() {}

  int      getData(const art::Event& ArtEvent);
  
  void         print_(int Level, const std::string&  Message,
                      const std::source_location& location = std::source_location::current());

  int      process_adc_waveform(float* Wf, WfParam_t* Wp);

  int      calculateMissingTrkParameters();
  
  int      makeSegments();

  int      fillFragments();

  int      fillCalD();
  int      fillCalH();
  int      fillCalC();

  int      fillCrvD();
  int      fillCrvP();
  int      fillCrvC();              // CRV coincidence clusters

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
  std::vector<std::string> _sDebugBits;
  int                      _debugBit[100];
  art::InputTag            _caldCollTag;        // CAL digi collection tag
  art::InputTag            _calhCollTag;        // CAL hit  collection tag
  art::InputTag            _calcCollTag;        // CAL cluster collection tag
  art::InputTag            _crvdCollTag;        // CRV digi collection tag
  art::InputTag            _crvpCollTag;        // CRV reco pulse collection tag
  art::InputTag            _crvcCollTag;        // CRV CC collection tag
  art::InputTag            _sdCollTag;          // straw digi collection tag
  art::InputTag            _shCollTag;          // straw hit collection tag
  art::InputTag            _tcCollTag;          // time cluster collection tag
  art::InputTag            _ksCollTag;          // kalseed collection tag
  std::string              _outputDir;
  int                      _saveWaveforms;
  int                      _makeCalD;
  int                      _makeCalH;
  int                      _makeCalC;
  int                      _makeCrvD;
  int                      _makeCrvP;
  int                      _makeCrvC;
  int                      _makeFragments;
  int                      _makeSD;
  int                      _makeSH;
  int                      _makeCH;
  int                      _makeTC;
  int                      _makeSeg;
  int                      _makeTrk;
  int                      _ewLength;           // it is up to the user to make sure it is set correctly
  int                      _nSamplesBL;
  int                      _minNSegments;
  float                    _minPulseHeight;
  float                    _minSDPHToSave;
  float                    _vDrift;
  float                    _tOffset;
    
  
  int                      _n_adc_samples;
  double                   _tdc_bin;            // TDC bin, in us
  double                   _tdc_bin_ns;         // TDC bin, ns

  DaqEvent*                _event;
  const art::Event*        _art_event;
  int                      _last_run;

  int                     _ncald;
  int                     _ncalh;
  int                     _ncalc;
  int                     _nstrawdigis;
  int                     _nstrawhits;
  int                     _ncombohits;
  int                     _ntimeclusters;
  int                     _nsegments;
  int                     _ntracks;
  
  int                     _ncalod;      // N(calo digis)

  int                     _ncrvd;       // N(CRV digis)
  int                     _ncrvp;       // N(CRV pulses)
  int                     _ncrvc;       // N(CRV coincidence clusters - muon stub candidates)

  int                     _nstmd;       // STM digis ?
                          
  TFile*                  _file;
  TTree*                  _tree;
  TBranch*                _branch;

  int                     _hist_booked;
   
  const mu2e::CaloDigiCollection*              _caldc;
  const mu2e::CaloHitCollection*               _calhc;
  const mu2e::CaloClusterCollection*           _calcc;
  const mu2e::CrvDigiCollection*               _crvdc;
  const mu2e::CrvRecoPulseCollection*          _crvpc;
  const mu2e::CrvCoincidenceClusterCollection* _crvcc;

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

//-----------------------------------------------------------------------------
// Level:     0: debug
//            1: info
//            2: warning
//            3: error
//------------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::print_(int Level, const std::string& Message, const std::source_location& location) {

  struct split {
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
  } xx;
  
  std::string s;
  if (_art_event) s = std::format("event: {}:{}:{} ",_art_event->run(),_art_event->subRun(),_art_event->event());

  std::vector<std::string> ss   = xx.splitString(location.file_name()    ,"/");
  std::vector<std::string> func = xx.splitString(location.function_name(),":");

  if (Level == e_DEBUG) {
    MF_LOG_VERBATIM(LOG_STREAM) << s << ss.back() << ":" << location.line() << ":" << func.back() << " : " << Message;
  } 
  else if (Level == e_INFO) {
    MF_LOG_PRINT(LOG_STREAM) << s << ss.back() << ":" << location.line() << ":" << func.back() << " : " << Message;
  }
  else if (Level == e_WARNING) {                // warning
    MF_LOG_PRINT(LOG_STREAM) << "WARNING: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }

  else if (Level == e_ERROR) {                // 
    MF_LOG_PROBLEM(LOG_STREAM) << "ERROR: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }

  else if (Level == e_SEVERE) {                // 
    MF_LOG_ABSOLUTE(LOG_STREAM) << "SEVERE: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }
}

// ======================================================================

mu2e::MakeDigiNtuple::MakeDigiNtuple(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _debugMode     (config().debugMode     ()),
    _sDebugBits    (config().debugBits     ()),

    _caldCollTag   (config().caldCollTag   ()),
    _calhCollTag   (config().calhCollTag   ()),
    _calcCollTag   (config().calcCollTag   ()),

    _crvdCollTag   (config().crvdCollTag   ()),
    _crvpCollTag   (config().crvpCollTag   ()),
    _crvcCollTag   (config().crvcCollTag   ()),

    _sdCollTag     (config().sdCollTag     ()),
    _shCollTag     (config().shCollTag     ()),
    _tcCollTag     (config().tcCollTag     ()),
    _ksCollTag     (config().ksCollTag     ()),
    _outputDir     (config().outputDir     ()),
    _saveWaveforms (config().saveWaveforms ()),
    _makeCalD      (config().makeCalD      ()),
    _makeCalH      (config().makeCalH      ()),
    _makeCalC      (config().makeCalC      ()),
    _makeCrvD      (config().makeCrvD      ()),
    _makeCrvP      (config().makeCrvP      ()),
    _makeCrvC      (config().makeCrvC      ()),
    _makeFragments (config().makeFragments ()),
    _makeSD        (config().makeSD        ()),
    _makeSH        (config().makeSH        ()),
    _makeCH        (config().makeCH        ()),
    _makeTC        (config().makeTC        ()),
    _makeSeg       (config().makeSeg       ()),
    _makeTrk       (config().makeTrk       ()),
    _ewLength      (config().ewLength      ()),
    _nSamplesBL    (config().nSamplesBL    ()),
    _minNSegments  (config().minNSegments  ()),
    _minPulseHeight(config().minPulseHeight()),
    _minSDPHToSave (config().minSDPHToSave ()),
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
  int nbits = _sDebugBits.size();
  for (int i=0; i<nbits; i++) {
    int index(0), value(0);
    key               = _sDebugBits[i].data();
    sscanf(key,"bit%i:%i",&index,&value);
    _debugBit[index]  = value;
    
    print_(e_INFO,std::format("...{}: bit={:4d} is set to {}\n",__func__,index,_debugBit[index]));
  }

  SegmentHit::SetVDrift (_vDrift );
  SegmentHit::SetTOffset(_tOffset);
  
  TrkSegment::fgDebugMode = _debugBit[3];
  SegmentFit::fgDebugMode = _debugBit[4];
}


//-----------------------------------------------------------------------------
void mu2e::MakeDigiNtuple::beginRun(const art::Run& ArtRun) {
  
  if (_hist_booked == 0) {
                                        // make sure we do it only once
    art::ServiceHandle<art::TFileService> tfs;
    TH1::AddDirectory(kFALSE);

  
    //  _file = new TFile(Form("%s/make_digi_ntuple_%06i.root",outputDir_.data(),ArtRun.run()),"RECREATE");
    TTree::SetMaxTreeSize(4000000000LL);

    _tree = tfs->make<TTree>("digis","digis");

    _event = new DaqEvent();

    _branch = _tree->Branch("evt","DaqEvent",_event,64000,99);
    _branch->SetAutoDelete(kFALSE);
    
    if (_branch) { 
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
        print_(e_DEBUG,std::format("-- HepTransform for plane:{:2}:{}\n",ipln,ipnl));
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
  bool ok;

  _sdc         = nullptr;
  _sdawfc      = nullptr;
  _nstrawdigis = 0;

  if (_makeSD) {
    art::Handle<mu2e::StrawDigiCollection>            sdch;
    art::Handle<mu2e::StrawDigiADCWaveformCollection> sdawfch;
    
    ok = ArtEvent.getByLabel(_sdCollTag,sdch);
    if (ok) { 
      _sdc         = sdch.product();
      _nstrawdigis = _sdc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: StrawDigiCollection:{:s} is not available. Bail out\n",
                                 _sdCollTag.encode().data()));
      return -1;
    }

    ok =  ArtEvent.getByLabel(_sdCollTag,sdawfch);
    if (ok) { 
      _sdawfc = sdawfch.product();
    }
    else {
      print_(e_WARNING,std::format("WARNING: StrawDigiADCWaveformCollection:{:s} is not available. Bail out\n",
                                   _sdCollTag.encode().data()));
      return -1;
    }
  }

  _shc = nullptr;
  _nstrawhits = 0;
  if (_makeSH) {
    art::Handle<mu2e::StrawHitCollection>             shch;
    ok = ArtEvent.getByLabel(_shCollTag,shch);
    if (ok) { 
      _shc         = shch.product();
      _nstrawhits = _shc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: StrawHitCollection:{:s} is not available. Bail out\n",_shCollTag.encode().data()));
      return -1;
    }
  }

  _tcc = nullptr;
  _ntimeclusters = 0;
  if (_makeTC != 0) {
    art::Handle<mu2e::TimeClusterCollection>          tcch;
    ok =  ArtEvent.getByLabel(_tcCollTag,tcch);
    if (ok) { 
      _tcc           = tcch.product();
      _ntimeclusters = _tcc->size();
    }
    else {
      print_(e_WARNING,std::format("WARNING: TimeClusterCollection:{:s} is not available. Bail out\n",
                       _tcCollTag.encode().data()));
      return -1;
    }
//-----------------------------------------------------------------------------
// assume that chCollTag == _shCollTag
//-----------------------------------------------------------------------------
    _chc        = nullptr;
    _ncombohits = 0;
    if (_makeCH != 0) {
      art::Handle<mu2e::ComboHitCollection> chch;
      ok =  ArtEvent.getByLabel(_shCollTag,chch);
      if (ok) { 
        _chc           = chch.product();
        _ncombohits    = _chc->size();
      }
      else {
        print_(e_WARNING,std::format("WARNING: ComboHitCollection:{:s} is not available. Bail out\n",
                                     _shCollTag.encode().data()));
        return -1;
      }
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  _ksc     = nullptr;
  _ntracks = 0;
  if (_makeTrk != 0) {
    art::Handle<mu2e::KalSeedCollection>  ksch;
    ok =  ArtEvent.getByLabel(_ksCollTag,ksch);
    if (ok) { 
      _ksc     = ksch.product();
      _ntracks = _ksc->size();
    }
    else {
      print_(e_WARNING,std::format("WARNING: KalSeedCollection:{:s} is not available. Bail out\n",
                         _ksCollTag.encode().data()));
      return -1;
    }
  }
//-----------------------------------------------------------------------------
// calorimeter
//-----------------------------------------------------------------------------
  _caldc = nullptr;
  _ncald = 0;
  if (_makeCalD != 0) {
    art::Handle<mu2e::CaloDigiCollection>     caldch;
    ok = ArtEvent.getByLabel(_caldCollTag,caldch);
    if (ok) { 
      _caldc = caldch.product();
      _ncald = _caldc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CaloDigiCollection:{:s} not found.\n",
                         _caldCollTag.encode().data()));
    }
  }
  
  _calhc = nullptr;
  _ncalh = 0;
  if (_makeCalH != 0) {
    art::Handle<mu2e::CaloHitCollection>     calhch;
    ok = ArtEvent.getByLabel(_calhCollTag,calhch);
    if (ok) { 
      _calhc = calhch.product();
      _ncalh = _calhc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CaloHitCollection:{:s} not found.\n",
                         _calhCollTag.encode().data()));
    }
  }
  
  _calcc = nullptr;
  _ncalc = 0;
  if (_makeCalC != 0) {
    art::Handle<mu2e::CaloClusterCollection>     calcch;
    ok = ArtEvent.getByLabel(_calcCollTag,calcch);
    if (ok) { 
      _calcc = calcch.product();
      _ncalc = _calcc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CaloClusterCollection:{:s} not found.\n",
                         _calcCollTag.encode().data()));
    }
  }
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  _crvdc = nullptr;
  _ncrvd = 0;
  if (_makeCrvD != 0) {
    art::Handle<mu2e::CrvDigiCollection>     crvdch;
    ok = ArtEvent.getByLabel(_crvdCollTag,crvdch);
    if (ok) { 
      _crvdc = crvdch.product();
      _ncrvd = _crvdc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CrvDigiCollection:{:s} not found.\n",
                         _crvdCollTag.encode().data()));
    }
  }
  
  _crvpc = nullptr;
  _ncrvp = 0;
  if (_makeCrvP != 0) {
    art::Handle<mu2e::CrvRecoPulseCollection>      crvpch;
    ok = ArtEvent.getByLabel(_crvpCollTag,crvpch);
    if (ok) { 
      _crvpc = crvpch.product();
      _ncrvp = _crvpc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CrvRecoPulseCollection:{:s} not found\n",
                         _crvpCollTag.encode().data()));
    }
  }

  _crvcc = nullptr;
  _ncrvc = 0;
  if (_makeCrvC != 0) {
    art::Handle<mu2e::CrvCoincidenceClusterCollection>     crvcch;
    ok = ArtEvent.getByLabel(_crvcCollTag,crvcch);
    if (ok) { 
      _crvcc = crvcch.product();
      _ncrvc = _crvcc->size();
    }
    else {
      print_(e_ERROR,std::format("ERROR: CrvCoincidenceClusterCollection:{:s} not found\n",
                         _crvcCollTag.encode().data()));
    }
  }
//-----------------------------------------------------------------------------
// bit_21: print CRV inputs:
//-----------------------------------------------------------------------------
  if ((_debugMode != 0) and (_debugBit[21])) {
    print_(e_DEBUG,std::format("bit21: _ncrvd:{} _ncrvp:{} _ncrvc:{}\n",_ncrvd,_ncrvp,_ncrvc));
  }
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
int mu2e::MakeDigiNtuple::fillCalD() {
  for (int i=0; i<_ncald; i++) {
    const mu2e::CaloDigi* cald = &_caldc->at(i);
    int ns = cald->waveform().size();

    DaqCaloDigi* nt_cald = (DaqCaloDigi*) _event->cald->ConstructedAt(i);
    nt_cald->Init(ns);
    
    nt_cald->sipmid       = cald->SiPMID();
    nt_cald->t0           = cald->t0();
    nt_cald->ppos         = cald->peakpos();
//-----------------------------------------------------------------------------
// store the waveform
//-----------------------------------------------------------------------------
    for (int is=0; is<ns; is++) {
      nt_cald->wf[is] = cald->waveform()[is];
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCalH() {
  _event->edisk[0] = 0;
  _event->edisk[1] = 0;
  
  for (int i=0; i<_ncalh; i++) {
    const mu2e::CaloHit* calh = &_calhc->at(i);

    DaqCaloHit* nt_calh = (DaqCaloHit*) _event->calh->ConstructedAt(i);
    nt_calh->cid        = calh->crystalID();
    nt_calh->nsipms     = calh->nSiPMs();
    nt_calh->time       = calh->time();
    nt_calh->sigt       = calh->timeErr();
    nt_calh->edep       = calh->energyDep();
    nt_calh->sige       = calh->energyDepErr();

    int disk = 0;
    if (calh->crystalID() >= 674) disk = 1;
    _event->edisk[disk] += nt_calh->edep;
  }
  
  _event->ecal = _event->edisk[0]+_event->edisk[1];
  
  return 0;
}

//-----------------------------------------------------------------------------
// calorimeter clusters
//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCalC() {
  
  for (int i=0; i<_ncalc; i++) {
    const mu2e::CaloCluster* calc = &_calcc->at(i);

    DaqCaloCluster* nt_calc = (DaqCaloCluster*) _event->calc->ConstructedAt(i);
    nt_calc->disk       = calc->diskID();
    nt_calc->size       = calc->size();
    nt_calc->time       = calc->time();
    nt_calc->sigt       = calc->timeErr();
    nt_calc->edep       = calc->energyDep();
    nt_calc->sige       = calc->energyDepErr();
    nt_calc->split      = calc->isSplit();
    nt_calc->x          = calc->cog3Vector().x();
    nt_calc->y          = calc->cog3Vector().y();
    nt_calc->z          = calc->cog3Vector().z();
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCrvD() {
  for (int i=0; i<_ncrvd; i++) {
    const mu2e::CrvDigi* crvd = &_crvdc->at(i);
    int ns = crvd->GetADCs().size();

    DaqCrvDigi* nt_crvd = (DaqCrvDigi*) _event->crvd->ConstructedAt(i);
    nt_crvd->Init(ns);
    
    nt_crvd->sbid         = crvd->GetScintillatorBarIndex().asInt();
    nt_crvd->tdc          = crvd->GetStartTDC();
    nt_crvd->nzs          = crvd->IsNZS();
    nt_crvd->odd_ts       = crvd->HasOddTimestamp();
    nt_crvd->sipm         = crvd->GetSiPMNumber();
    nt_crvd->roc          = crvd->GetROC();
    nt_crvd->feb          = crvd->GetFEB();
    nt_crvd->ch           = crvd->GetFEBchannel();
//-----------------------------------------------------------------------------
// store the waveform
//-----------------------------------------------------------------------------
    // for (int is=0; is<ns; is++) {
    //   nt_crvd->adc[is] = crvd->GetADCs()[is];
    // }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCrvP() {
  for (int i=0; i<_ncrvp; i++) {
    const mu2e::CrvRecoPulse* crvp = &_crvpc->at(i);

    DaqCrvRecoPulse* nt_crvp = (DaqCrvRecoPulse*)_event->crvp->ConstructedAt(i);
    nt_crvp->Init();
    
    nt_crvp->npes            = crvp->GetPEs();
    
    nt_crvp->pes_ph          = crvp->GetPEsPulseHeight();

    nt_crvp->time            = crvp->GetPulseTime();
    nt_crvp->ph              = crvp->GetPulseHeight();
    nt_crvp->ped             = crvp->GetPedestal();
    nt_crvp->beta            = crvp->GetPulseBeta();
    nt_crvp->chi2            = crvp->GetPulseFitChi2();
    nt_crvp->le_time         = crvp->GetLEtime();
    nt_crvp->flags           = crvp->GetRecoPulseFlags().to_ulong();
    
    nt_crvp->npes_nofit      = crvp->GetPEsNoFit();
    nt_crvp->time_nofit      = crvp->GetPulseTimeNoFit();
    nt_crvp->tstart          = crvp->GetPulseStart();
    nt_crvp->tend            = crvp->GetPulseEnd();

    nt_crvp->sbid            = crvp->GetScintillatorBarIndex().asInt();
    nt_crvp->sipm            = crvp->GetSiPMNumber();
    nt_crvp->roc             = crvp->GetROC();
    nt_crvp->feb             = crvp->GetFEB();
    nt_crvp->ch              = crvp->GetFEBchannel();
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillCrvC() {
  for (int i=0; i<_ncrvc; i++) {
    const mu2e::CrvCoincidenceCluster* crvc = &_crvcc->at(i);

    DaqCrvCoincidenceCluster* nt_crvc = (DaqCrvCoincidenceCluster*) _event->crvc->ConstructedAt(i);
    nt_crvc->Init();

    nt_crvc->stype   = crvc->GetCrvSectorType();
    nt_crvc->tstart  = crvc->GetStartTime();
    nt_crvc->tend    = crvc->GetEndTime();
    nt_crvc->pes     = crvc->GetPEs();
    nt_crvc->time    = crvc->GetAvgHitTime();
    nt_crvc->x       = crvc->GetAvgHitPos().x();
    nt_crvc->y       = crvc->GetAvgHitPos().y();
    nt_crvc->z       = crvc->GetAvgHitPos().z();
    nt_crvc->nlayers = crvc->GetLayers().size();
    nt_crvc->nsides  = (crvc->HasTwoReadoutSides() == 0) ? 1 : 2;
  }
  return 0;
}

//-----------------------------------------------------------------------------
// called always
//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillFragments() {
  int rc(0);
  print_(e_DEBUG,"--START:");
  artdaq::Fragments    fragments;
  artdaq::FragmentPtrs containerFragments;

  auto fragmentHandles = _art_event->getMany<std::vector<artdaq::Fragment>>();
  _event->nfrag = 0;
  
  if (_debugMode > 0) {
    std::string msg = std::format("n_fragment_collections):{}",fragmentHandles.size());
    print_(e_DEBUG,msg);
  }

  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->empty())     continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      // not sure what this is....
      for (const auto& cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        for (size_t ii = 0; ii < contf.block_count(); ++ii) {
          containerFragments.push_back(contf[ii]);
          fragments.push_back(*containerFragments.back());
        }
      }
    }
    else {
      // 
      int n_fragments = handle->size();

      _event->nfrag += n_fragments;
      
      if (_debugMode) {
        print_(e_DEBUG,std::format("-- next fragment collection with n_fragments:{}",n_fragments));
      }
      
      for (int ifrag=0; ifrag<n_fragments; ifrag++) {
        const artdaq::Fragment* frag = &handle->at(ifrag);

        if (_debugMode and (_debugBit[0] > 0)) {
          print_(e_DEBUG,std::format("-- fragment number:{} version:{} timestamp:{} data_size:{} type:{} DTC_SubEventHeader.size:{}",
                             ifrag,frag->version(),frag->timestamp(),frag->dataSizeBytes(),
                             frag->typeString(),sizeof(DTCLib::DTC_SubEventHeader)));
          //          print_fragment(frag);
        }
//-----------------------------------------------------------------------------
// skip CFO fragment (type = 12)
//-----------------------------------------------------------------------------
        if (frag->type() == mu2e::FragmentType::CFO)        continue;
        uint8_t* fdata = (uint8_t*) (frag->dataBegin());
//-----------------------------------------------------------------------------
// skip fragments with the payload size less than the DTC header size
// do it only for the current data format (runs > 107236)
//-----------------------------------------------------------------------------
        if (frag->dataSizeBytes() <= sizeof(DTCLib::DTC_SubEventHeader)) {
          std::string msg = std::format("ERROR: fragment:{} data size:{} < DTC_SubEventHeader.size:{}. SKIP FRAGMENT",
                                        ifrag,frag->dataSizeBytes(),sizeof(DTCLib::DTC_SubEventHeader));
          print_(e_DEBUG,msg);
          continue;
        }
        fdata += sizeof(DTCLib::DTC_EventHeader);
//-----------------------------------------------------------------------------
// skip non-tracker fragments
// after a recent format change, a DTC fragment may contain ROC data from different
// subdetectors, make sure that at least one of them is the tracker ROC
//-----------------------------------------------------------------------------
        DTCLib::DTC_SubEventHeader* seh = (DTCLib::DTC_SubEventHeader*) fdata;
        if ((seh->link0_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link1_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link2_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link3_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link4_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link5_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker)     )
                                                            continue;
        // so far store only tracker fragments
        DaqFragment* nt_fr    = _event->NewFragment(ifrag);
        nt_fr->nbytes         = seh->inclusive_subevent_byte_count;
        nt_fr->ewtag          = ((uint64_t) seh->event_tag_low) | (((uint64_t) seh->event_tag_high) << 32);
        nt_fr->nrocs          = seh->num_rocs;
        nt_fr->event_mode     = seh->event_mode;
        nt_fr->mac_addr       = seh->dtc_mac;
        nt_fr->partition      = seh->partition_id;
        nt_fr->evb_mode       = seh->evb_mode;
        nt_fr->dtc_id         = seh->source_dtc_id;
        nt_fr->link_ssid[0]   = seh->link0_subsystem;
        nt_fr->link_ssid[1]   = seh->link1_subsystem;
        nt_fr->link_ssid[2]   = seh->link2_subsystem;
        nt_fr->link_ssid[3]   = seh->link3_subsystem;
        nt_fr->link_ssid[4]   = seh->link4_subsystem;
        nt_fr->link_ssid[5]   = seh->link5_subsystem;
        nt_fr->link_status[0] = seh->link0_status;
        nt_fr->link_status[1] = seh->link1_status;
        nt_fr->link_status[2] = seh->link2_status;
        nt_fr->link_status[3] = seh->link3_status;
        nt_fr->link_status[4] = seh->link4_status;
        nt_fr->link_status[5] = seh->link5_status;
        nt_fr->version        = seh->subevent_format_version;
        nt_fr->emtdc          = seh->emtdc;
        nt_fr->latency[0]     = seh->link0_drp_rx_latency;
        nt_fr->latency[1]     = seh->link1_drp_rx_latency;
        nt_fr->latency[2]     = seh->link2_drp_rx_latency;
        nt_fr->latency[3]     = seh->link3_drp_rx_latency;
        nt_fr->latency[4]     = seh->link4_drp_rx_latency;
        nt_fr->latency[5]     = seh->link5_drp_rx_latency;
//-----------------------------------------------------------------------------
// this is a tracker DTC fragment, loop over the ROCs
//-----------------------------------------------------------------------------
        // ushort*  buf          = (ushort*) fdata;
        // int      nbytes       = buf[0];
        uint8_t* roc_data     = fdata+sizeof(*seh);
        // uint8_t* last_address = fdata+nbytes;
        for (int lnk=0; lnk<6; lnk++) {
          RocDataHeaderPacket_t* rdh = (RocDataHeaderPacket_t*) roc_data;
          DaqRocData* nt_rd   = &nt_fr->roc[lnk];
          nt_rd->nbytes       = rdh->byteCount;
          nt_rd->ewtag        = ((uint64_t) rdh->eventTag[0]) | (((uint64_t) rdh->eventTag[1]) << 16) | (((uint64_t) rdh->eventTag) << 32);
          nt_rd->packet_type = rdh->packetType;
          nt_rd->link        = rdh->linkID;
          nt_rd->err         = rdh->DtcErrors;
          nt_rd->valid       = rdh->valid;
          nt_rd->npackets    = rdh->packetCount;
          nt_rd->ssid        = rdh->subsystemID;
          nt_rd->status      = rdh->status;
          nt_rd->version     = rdh->version;
          nt_rd->dtc_id      = rdh->dtcID;
          nt_rd->onspill     = rdh->onSpill;
          nt_rd->subrun      = rdh->subrun;
          nt_rd->event_mode  = rdh->eventMode;
          roc_data           = roc_data+rdh->byteCount;
        }
      }
    }
  }
  print_(e_DEBUG,std::format("--END: rc:{}",rc));
  return rc;
}


//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSD() {

  int idigi=0;
  for (int i=0; i<_nstrawdigis; i++) {
    const mu2e::StrawDigi*            sd    = &_sdc->at(i);
    const mu2e::StrawDigiADCWaveform* sdawf = &_sdawfc->at(i);
    int ns = sdawf->samples().size();
                                        // one-time initializatiion
    if (_n_adc_samples == -1) _n_adc_samples = ns;
//-----------------------------------------------------------------------------
// process the waveform and store the waveform parameters
//-----------------------------------------------------------------------------
    float wf[100];
    for (int is=0; is<ns; is++) {
      wf[is] = sdawf->samples()[is];
    }

    WfParam_t wp;
    process_adc_waveform(wf,&wp);
//-----------------------------------------------------------------------------
// for calibrations, write out a special small size ntuple with the cut on the SD pulse height
// by default, _minSdPulseHeight is -1.
//-----------------------------------------------------------------------------
    if (wp.ph < _minSDPHToSave)                             continue;
    DaqStrawDigi* nt_sd = _event->NewSD(idigi);
    idigi++;
    
    nt_sd->InitSD(ns);
    
    nt_sd->sid          = sd->strawId().asUint16();
 
    int pln = sd->strawId().plane();
    int pnl = sd->strawId().panel();
    // int ich = sd->strawId().straw();

    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);

    int dtc_id    = tpm->dtc();

    nt_sd->mnid         = tpm->mnid();
    
    nt_sd->tdc0         = sd->TDC(mu2e::StrawEnd::cal);
    nt_sd->tdc1         = sd->TDC(mu2e::StrawEnd::hv );
    nt_sd->tot0         = sd->TOT(mu2e::StrawEnd::cal);
    nt_sd->tot1         = sd->TOT(mu2e::StrawEnd::hv );
    nt_sd->pmp          = sd->PMP();
    // dtc_id runs from 1 to 36
    if (_event->pmp[dtc_id-1] == -1) {
      _event->pmp[dtc_id-1] = nt_sd->pmp;
    }
    else if (_event->pmp[dtc_id-1] !=  nt_sd->pmp) {
      print_(e_ERROR,std::format("dtc_id:{} _event->pmp[dtc_id]:{} hit_pmp{}",dtc_id,_event->pmp[dtc_id-1],nt_sd->pmp));
    }
    
    nt_sd->flag         = *((uint8_t*) &sd->digiFlag());
//-----------------------------------------------------------------------------
// store the waveform
//-----------------------------------------------------------------------------
    for (int is=0; is<ns; is++) {
      nt_sd->adc[is] = sdawf->samples()[is];
    }
    
    nt_sd->fs = wp.fs;
    nt_sd->bl = wp.bl;
    nt_sd->ph = wp.ph;

    if ((_debugMode  > 0) and (_debugBit[1] != 0)) {
//-----------------------------------------------------------------------------
// for all hits, print hit times assuming contiguous timing
//-----------------------------------------------------------------------------
      double t0_offset  = _event->evn*_ewLength*25;                  // in ns
      double t0         = t0_offset + nt_sd->tdc0*_tdc_bin_ns;
      double t1         = t0_offset + nt_sd->tdc1*_tdc_bin_ns;
      printf("%8i %5i MN%03i %8i %8i %12.4lf %12.4lf %6i %6i %6i 0x%04x\n",
             _event->evn,
             (int) nt_sd->sid, nt_sd->mnid, 
             nt_sd->tdc0, nt_sd->tdc1,
             t0, t1,
             nt_sd->tot0, nt_sd->tot1,
             nt_sd->pmp, nt_sd->flag);
    }
  }
  _nstrawdigis = idigi;
  return 0;
}

//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillSH() {

  if ((_debugMode > 0) and (_debugBit[1] != 0)) {
    printf("evn    sid  pln  pnl mnid    time    dt   tot0 tot1   edep\n");
    printf("---------------------------------------------------------\n");
  }
  for (int i=0; i<_nstrawhits; i++) {
    const mu2e::StrawHit* sh = &_shc->at(i);
    int pln = sh->strawId().plane();
    int pnl = sh->strawId().panel();
    //    const TrkPanelMap_t* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);
    const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(pln,pnl);

    int dtc_id = tpm->dtc();
    //    int pcie_addr = dtc_id % 2;                      // convention
    _event->nsh[dtc_id-1][tpm->link()] += 1;         // assume DTC_ID runs from 1 to 36 (tracker)
   
    DaqStrawHit* nt_sh = new ((*_event->sh)[i]) DaqStrawHit();
    nt_sh->sid         = sh->strawId().asUint16();
    nt_sh->zface       = tpm->zface();
    nt_sh->mnid        = tpm->mnid();
    nt_sh->time        = sh->time(mu2e::StrawEnd::cal);
    nt_sh->dt          = sh->dt();      // cal-hv
    nt_sh->tot0        = sh->TOT(mu2e::StrawEnd::cal);
    nt_sh->tot1        = sh->TOT(mu2e::StrawEnd::hv );
    nt_sh->edep        = sh->energyDep();
    nt_sh->dped        = sh->digitalPedestal();
    nt_sh->dpmp        = sh->digitalPulseHeight();
    if (sh->energyDep() > _event->maxEdep) _event->maxEdep = sh->energyDep();

    if ((_debugMode  > 0) and (_debugBit[1] != 0)) {
      printf("%8i %5i %3i %3i %3i %12.4f %12.4f %6.1f %6.1f %7.4f\n",
             _event->evn,
             (int) nt_sh->sid, pln, pnl, nt_sh->mnid,
             nt_sh->time, nt_sh->dt,
             nt_sh->tot0, nt_sh->tot1,
             nt_sh->edep);
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
    nt_ch->edep        = ch->energyDep();

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
// time clusters are filled using offline indices
//-----------------------------------------------------------------------------
int mu2e::MakeDigiNtuple::fillTC() {

  for (int itc=0; itc<_ntimeclusters; itc++) {
    const mu2e::TimeCluster* tc = &_tcc->at(itc);
   
    DaqTimeCluster* nt_tc = new ((*_event->tc)[itc]) DaqTimeCluster();

    nt_tc->nsh     = tc->nStrawHits();
    nt_tc->nch     = tc->nhits();
    nt_tc->t0      = tc->t0().t0();
    nt_tc->ngh     = 0;
    
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
      int stn   = plane / 2;
      int pln2  = plane % 2;            // 0 or 1
      int panel = hit->strawId().panel();
      const TrkPanelMap::Row* tpm = _trkPanelMap->panel_map_by_offline_ind(plane,panel);
      int zface = tpm->zface();

      int loc = 6*pln2+panel;
      
      if (nt_tc->_nh_panel[stn][loc] == 0) nt_tc->npanels++;
      nt_tc->_nh_panel  [stn][loc]++;
      nt_tc->_time_panel[stn][loc] += hit->correctedTime();
      nt_tc->_edep_panel[stn][loc] += hit->energyDep();

      if (hit->energyDep() > nt_tc->edep_max) nt_tc->edep_max = hit->energyDep();

      if (hit->energyDep() > 0.0005) {
        nt_tc->ngh += 1;
      }

      if (nt_tc->_nhp[stn][pln2] == 0) nt_tc->nplanes++;
      nt_tc->_nhp  [stn][pln2]++;
      nt_tc->_timep[stn][pln2] += hit->correctedTime();

      if (nt_tc->_nhf[stn][zface] == 0) {
        nt_tc->nfaces++;
        nt_tc->_mnid[stn][zface] = tpm->mnid();
      }
      nt_tc->_nhf  [stn][zface]++;
      nt_tc->_timef[stn][zface] += hit->correctedTime();

      if (hit->correctedTime() < nt_tc->tmin) nt_tc->tmin = hit->correctedTime();
      if (hit->correctedTime() > nt_tc->tmax) nt_tc->tmax = hit->correctedTime();
    }
//-----------------------------------------------------------------------------
// done looping over the time cluster hits,
//-----------------------------------------------------------------------------
    nt_tc->y0     = szy.y0();
    nt_tc->dydz   = szy.dydx();
    nt_tc->chi2yz = szy.chi2Dof();

    for (int stn=0; stn<18; stn++) {
      for (int ip=0; ip<2; ip++) {
        if (nt_tc->_nhp[stn][ip] > 0) nt_tc->_timep[stn][ip] = nt_tc->_timep[stn][ip]/nt_tc->_nhp[stn][ip];
      }

      for (int i=0; i<4; i++) {
        if (nt_tc->_nhf[stn][i] > 0) nt_tc->_timef[stn][i] = nt_tc->_timef[stn][i]/nt_tc->_nhf[stn][i];
      }
//-----------------------------------------------------------------------------
// average time and charge of the hits in a given panel
//-----------------------------------------------------------------------------
      for (int i=0; i<12; i++) {
        if (nt_tc->_nh_panel[stn][i] > 0) {
          nt_tc->_time_panel[stn][i] = nt_tc->_time_panel[stn][i]/nt_tc->_nh_panel[stn][i];
          nt_tc->_edep_panel[stn][i] = nt_tc->_edep_panel[stn][i]/nt_tc->_nh_panel[stn][i];
        }

        if (nt_tc->_nh_panel[stn][i] > nt_tc->max_nh_panel) {
          nt_tc->max_nh_panel = nt_tc->_nh_panel[stn][i];
        }
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

  if ((_debugMode != 0) and (_debugBit[31] != 0)) {
    print_(e_DEBUG,std::format("{}  START: ntimeclusters:{}\n", __func__,_ntimeclusters));
  }
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
    if (_debugMode and (_debugBit[31] != 0)) {
      print_(e_DEBUG,std::format(" -- PM: time cluster:{:2d} segment with {:2d} hits\n",itc, nsh));
    }
    for (int ih=0; ih<nsh; ih++) {
      StrawHitIndex hit_index   = tc->hits().at(ih);
      const mu2e::ComboHit* ch = &_chc->at(hit_index);
      mu2e::StrawId sid = ch->strawId();
      int panel = sid.panel();
      int plane = sid.plane();
      if (_debugMode and (_debugBit[31] != 0)) {
        print_(e_DEBUG,std::format("hit number:{:3d} plane:{} panel:{} straw:{:2d}\n",
                                   ih,plane,panel,sid.straw()));
      }
//-----------------------------------------------------------------------------
// cosmic track: assume one segment per panel, in principle, there could be more than one
//-----------------------------------------------------------------------------
      TrkSegment* ts = &_tseg[plane][panel];
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
    int first_hit;                      // index in the segment hit list
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
    
    if ((_debugMode != 0) and (_debugBit[31] != 0)) {
      print_(e_DEBUG,std::format("{}  iseg:{} nhits:{}\n",__func__,i,nhits));
    }

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
    if ((_debugMode != 0) and (_debugBit[31] != 0)) {
      print_(e_DEBUG,std::format("{}  iseg:{} nhits:{}\n",__func__,i,ts->nHits()));
    }
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

    if (_debugMode and (_debugBit[31] != 0)) {
      std::cout << "rc:" << rc << " converged:" << converged;
    }
  }
//-----------------------------------------------------------------------------
// at this point, all track hits should be assigned to segments
// debug printout
//-----------------------------------------------------------------------------
  if (_debugMode and (_debugBit[31] != 0)) {
    std::cout << "nseg:" << _nseg << std::endl;

    for (int i=0; i<_nseg; i++) {
      TrkSegment* ts = _ptseg[i];
      ts->print();
    }
  }

  if ((_debugMode != 0) and (_debugBit[31] != 0)) {
    print_(e_DEBUG,std::format("-- END rc:{}\n",rc));
  }
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
    if (_debugMode and (_debugBit[31] != 0)) {
      print_(e_DEBUG,std::format(" iseg:{} dz/dy(seg):{:12.5f} dz/dy(trk):{:12.5f}\n",iseg,nt_ts->dzdy,nt_ts->dzdyt));
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
    print_(e_DEBUG,std::format("-- START event:{}:{}:{}\n",ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event()));
  }

  int rc = getData(ArtEvent);
  if (rc < 0) return;
//-----------------------------------------------------------------------------
// initialize event to "empty" values
//-----------------------------------------------------------------------------
  _event->Clear();
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
  _event->run     = ArtEvent.run();
  _event->srn     = ArtEvent.subRun();
  _event->evn     = ArtEvent.event();
  _event->ewtag   = -1;                 // to be figured separately from either 
                                        // defined in getData()
  _event->nsdtot  = _nstrawdigis;
  _event->nshtot  = _nstrawhits;
  _event->nch     = _ncombohits;
  _event->ntc     = _ntimeclusters;
  _event->ntrk    = _ntracks;

  _event->ncald   = _ncald;
  _event->ncalh   = _ncalh;
  _event->ncalc   = _ncalc;

  _event->ncrvd   = _ncrvd;
  _event->ncrvp   = _ncrvp;
  _event->ncrvc   = _ncrvc;
  //  _event->nfrag   = _nfrag;
  
  _event->maxEdep = 0;

  if ((_debugMode > 0) and (_debugBit[11] != 0)) {
    print_(e_DEBUG,std::format("_nstrawdigis:{}\n",_nstrawdigis));
  }

  if (_makeFragments) fillFragments();
  
  if (_makeSD       ) fillSD ();
  if (_makeSH       ) fillSH ();
  if (_makeCH       ) fillCH ();
  if (_makeTC       ) fillTC ();

  if (_makeCalD     ) fillCalD ();
  if (_makeCalH     ) fillCalH ();
  if (_makeCalC     ) fillCalC ();

  if (_makeCrvD     ) fillCrvD ();
  if (_makeCrvP     ) fillCrvP ();
  if (_makeCrvC     ) fillCrvC ();

  if (_makeSeg) {
    fillSeg  ();
    fillSegSh();
  }
  if (_makeTrk      ) fillTrk();

  if (_event->nseg >= _minNSegments) {
    _tree->Fill();
  }

  if (_debugMode > 0) print_(e_DEBUG,"-- END\n");
}



// ======================================================================

DEFINE_ART_MODULE(mu2e::MakeDigiNtuple)

// ======================================================================
