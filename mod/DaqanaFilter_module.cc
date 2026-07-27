//
//  PM: straw hit filter, reject the noise
//  use single straw combo hits
//  StrawHitReco module always produces single-straw combohit collection - use that
//
//  debugBit[0] : print all WFs
//  debugBit[1] : print good WFs
#include "TH1F.h"
//#include "TFolder.h"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
// mu2e
// data
// #include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "TRACE/tracemf.h"
#define TRACE_NAME "DaqAnaFilter"

#define LOG_STREAM "DaqAnaFilter"

// c++
#include <iostream>
#include <memory>
#include <map>

namespace mu2e {
  class DaqanaFilter : public art::EDFilter {
  public:

    enum {
      e_DEBUG   = 0,
      e_INFO    = 1,
      e_WARNING = 2,
      e_ERROR   = 3,
      e_SEVERE  = 4,
    };

    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool>            debugMode        {Name("debugMode"        ), Comment("Debug Mode, default:false")};
      fhicl::Sequence<std::string> debugBits        {Name("debugBits"        ), Comment("debug bits")};
      fhicl::Atom<art::InputTag>   chCollTag        {Name("chCollTag"        ), Comment("ComboHitColl tag") };
      fhicl::Atom<art::InputTag>   crvcCollTag      {Name("crvcCollTag"      ), Comment("CRV CC coll tag") };
      fhicl::Atom<art::InputTag>   tcCollTag        {Name("tcCollTag"        ), Comment("TimeClusterCollection tag") };
      fhicl::Atom<bool>            requireCrvc      {Name("requireCrvc"      ), Comment("1: need CRVC, 0:ignore absense")};
      
      fhicl::Atom<float>           maxDt            {Name("maxDt"            ), Comment("max abs(DT)")};
      fhicl::Atom<float>           maxTcCrvcDt      {Name("maxTcCrvcDt"      ), Comment("max TC-CRVC deltaT")};
      fhicl::Atom<float>           minEDep          {Name("minEDep"          ), Comment("min EDep")};
      fhicl::Atom<int>             minNGoodCh       {Name("minNGoodCh"       ), Comment("min N good combo hits")};
      fhicl::Atom<int>             minNTimeClusters {Name("minNTimeClusters" ), Comment("min N time clusters")};
      fhicl::Atom<int>             minNTcCrvcMatches{Name("minNTcCrvcMatches"), Comment("min N TC-CRVC matches")};
      fhicl::Atom<bool>            fillHistograms   {Name("fillHistograms"   ), Comment("fill histogrms, default:false")};
    };

    enum {
      kMaxEvtHistSets = 10,
      kMaxChHistSets  = 10,
      kMaxTcHistSets  = 10,
    };

    struct EventHist_t {
      TH1F* h_evt;
      TH1F* h_ncht;
      TH1F* h_nchg;
      TH1F* h_ntc;
      TH1F* h_ntcg;
    };

    struct ComboHitHist_t {
      TH1F* h_dt10;   // HV-CAL
      TH1F* h_edep;
    };

    struct TimeClusterHist_t {
      TH1F* h_nh;
      TH1F* h_nhg;
      TH1F* h_edep;
      TH1F* h_dt;
    };

    struct Hist_t {
      ComboHitHist_t*    ch [kMaxChHistSets];
      TimeClusterHist_t* tc [kMaxTcHistSets];
      EventHist_t*       evt[kMaxEvtHistSets];
    } _hist;
    
    struct Index_t {
      int rn;
      int slot;                           // 0-17
      int plane;                          // offline
      int panel;                          // offline
      int pnl12;                          // panel index within the station (0-11)
      int mnid;
      int ch;
    };
    
    using Parameters = art::EDFilter::Table<Config>;

    explicit DaqanaFilter(const Parameters& config);

    
    int book_evt_histograms(EventHist_t*       Hist, Index_t* Index, art::TFileDirectory* Dir);
    int book_ch_histograms (ComboHitHist_t*    Hist, Index_t* Index, art::TFileDirectory* Dir);
    int book_tc_histograms (TimeClusterHist_t* Hist, Index_t* Index, art::TFileDirectory* Dir);
    int book_histograms(int RunNumber);
    
    int fill_histograms(Hist_t* Hist);
    
    void      print_(int                         Level,
                     const std::string&          Message,
                     const std::source_location& location = std::source_location::current());

  private:
    bool filter  (art::Event& ArtEvent) override;
    bool beginRun(art::Run&   Run     ) override;
    bool endRun  (art::Run&   Run     ) override;

    art::InputTag              _chCollTag;
    art::InputTag              _crvcCollTag;
    art::InputTag              _tcCollTag;

    bool                     _debugMode;
    std::vector<std::string> _debugBits;
    int                      _debugBit[100];
    bool                     _requireCrvc;
    
    float                    _maxDt;
    float                    _minEDep;
    int                      _minNGoodCh;
    int                      _minNTimeClusters;
    int                      _minNTcCrvcMatches;

    int                      _fillHistograms;
    
    int                      _nevt;
    int                      _nevp;
    int                      _ncht;
    int                      _nchg;
    int                      _ntc;
    int                      _ntcg;     // N(time clusters with at least one hit above minEDep)
    int                      _ncrvc;
    int                      _n_tc_crvc_matches;

    std::vector<int>         _tc_mask;

    const mu2e::ComboHitCollection*    _chc;
    const mu2e::TimeClusterCollection* _tcc;
    
    // const mu2e::CrvDigiCollection*               _crvdc;
    // const mu2e::CrvRecoPulseCollection*          _crvpc;
    const mu2e::CrvCoincidenceClusterCollection* _crvcc;

    const art::Event* _art_event;
    int               _rn;
    bool              _run_initialized;
  };

//-----------------------------------------------------------------------------
  DaqanaFilter::DaqanaFilter(const Parameters& conf)
    : art::EDFilter{conf},
      _chCollTag         (conf().chCollTag()),
      _crvcCollTag       (conf().crvcCollTag()),
      _tcCollTag         (conf().tcCollTag()),
      _debugMode         (conf().debugMode()),
      _debugBits         (conf().debugBits()),
      _requireCrvc       (conf().requireCrvc()),
      _maxDt             (conf().maxDt    ()),
      _minEDep           (conf().minEDep  ()),
      _minNGoodCh        (conf().minNGoodCh()),
      _minNTimeClusters  (conf().minNTimeClusters()),
      _minNTcCrvcMatches (conf().minNTcCrvcMatches()),
      _fillHistograms    (conf().fillHistograms())
  {
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
    }
    _run_initialized = false;
  }

//-----------------------------------------------------------------------------
  void mu2e::DaqanaFilter::print_(int Level, const std::string& Message, const std::source_location& location) {

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


//-----------------------------------------------------------------------------
  int DaqanaFilter::book_evt_histograms(EventHist_t* Hist, Index_t* Index, art::TFileDirectory* Dir) {

    std::string prefix = std::format("run:{:6d}",Index->rn);
    std::string name, title;
     
    title = std::format("{} event number",prefix);
    Hist->h_evt  = Dir->make<TH1F>("evt" , title.data(), 10000,   0., 1.e7 );

    title = std::format("{} nch total",prefix);
    Hist->h_ncht = Dir->make<TH1F>("ncht", title.data(),   500, -0.5, 499.5);

    title = std::format("{} nchg",prefix);
    Hist->h_nchg = Dir->make<TH1F>("nchg", title.data(),   100, -0.5,  99.5);

    title = std::format("{} N(time clusters)",prefix);
    Hist->h_ntc  = Dir->make<TH1F>("ntc" , title.data(),   100, -0.5,  99.5);
    
    title = std::format("{} N(time clusters with at least {} good hits )",prefix,_minNGoodCh);
    Hist->h_ntcg = Dir->make<TH1F>("ntcg", title.data(),   100, -0.5,  99.5);
    
    return 0;
  }

//-----------------------------------------------------------------------------
  int DaqanaFilter::book_ch_histograms(ComboHitHist_t* Hist, Index_t* Index, art::TFileDirectory* Dir) {

    std::string prefix = std::format("run:{:6d}",Index->rn);
    std::string name, title;
     
    title = std::format("{} dt10 = T(HV)-T(CAL)",prefix);
    Hist->h_dt10 = Dir->make<TH1F>("dt10", title.data(),   500,  -25,   25. );

    title = std::format("{} edep, keV",prefix);
    Hist->h_edep = Dir->make<TH1F>("edep", title.data(),   100,   -2,   8. );

    return 0;
  }
  
//-----------------------------------------------------------------------------
  int DaqanaFilter::book_tc_histograms(TimeClusterHist_t* Hist, Index_t* Index, art::TFileDirectory* Dir) {
    std::string prefix = std::format("run:{:6d}",Index->rn);
    std::string name, title;
     
    title = std::format("{} nh",prefix);
    Hist->h_nh   = Dir->make<TH1F>("nh"  , title.data(),   200, -0.5, 199.5 );
    
    title = std::format("{} nhg",prefix);
    Hist->h_nhg  = Dir->make<TH1F>("nhg" , title.data(),   200, -0.5, 199.5 );
    
    title = std::format("{} dt = T(hit)-T0",prefix);
    Hist->h_dt   = Dir->make<TH1F>("dt"  , title.data(),   500, -250., 250. );
    
    title = std::format("{} edep, keV",prefix);
    Hist->h_edep = Dir->make<TH1F>("edep", title.data(),   100,  -2.,   8. );
    return 0;
  }
  
//-----------------------------------------------------------------------------
  int DaqanaFilter::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;

     TH1::AddDirectory(kFALSE);
   
     Index_t index;
     index.rn = RunNumber;

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
     int book_evt_histset[kMaxEvtHistSets];
     for (int i=0; i<kMaxEvtHistSets; i++) { book_evt_histset[i] = 0; }

     book_evt_histset[0] = 1;
     book_evt_histset[1] = 1;

     for (int i=0; i<kMaxEvtHistSets; i++) {
       if (book_evt_histset[i] == 0) continue;
       std::string folder_name = std::format("evt_{:02d}",i);
       art::TFileDirectory dir = tfs->mkdir(folder_name);
       _hist.evt[i]  = new EventHist_t;
       book_evt_histograms(_hist.evt[i],&index,&dir);
     }
//-----------------------------------------------------------------------------
// combo hit histograms
//-----------------------------------------------------------------------------
     int book_ch_histset[kMaxChHistSets];
     for (int i=0; i<kMaxChHistSets; i++) { book_ch_histset[i] = 0; }

     book_ch_histset[0] = 1;
     book_ch_histset[1] = 1;

     for (int i=0; i<kMaxChHistSets; i++) {
       if (book_ch_histset[i] == 0) continue;
       std::string folder_name = std::format("ch_{:02d}",i);
       art::TFileDirectory dir = tfs->mkdir(folder_name);
       _hist.ch[i]  = new ComboHitHist_t;
       book_ch_histograms(_hist.ch[i],&index,&dir);
     }
//-----------------------------------------------------------------------------
// time cluster histograms
//-----------------------------------------------------------------------------
     int book_tc_histset[kMaxTcHistSets];
     for (int i=0; i<kMaxTcHistSets; i++) { book_tc_histset[i] = 0; }

     book_tc_histset[0] = 1;

     for (int i=0; i<kMaxTcHistSets; i++) {
       if (book_tc_histset[i] == 0) continue;
       std::string folder_name = std::format("tc_{:02d}",i);
       art::TFileDirectory dir = tfs->mkdir(folder_name);
       _hist.tc[i]  = new TimeClusterHist_t;
       book_tc_histograms(_hist.tc[i],&index,&dir);
     }

    return 0;
  }
  
//-----------------------------------------------------------------------------
  int DaqanaFilter::fill_histograms(Hist_t* Hist) {

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
    Hist->evt[0]->h_evt->Fill(_art_event->event());
    Hist->evt[0]->h_ncht->Fill(_ncht);
    Hist->evt[0]->h_nchg->Fill(_nchg);
    Hist->evt[0]->h_ntc->Fill(_ntc);
    Hist->evt[0]->h_ntcg->Fill(_ntcg);
    
//-----------------------------------------------------------------------------
// combo hit histograms
//-----------------------------------------------------------------------------
    for (int i = 0; i<_ncht; ++i) {
      const mu2e::ComboHit* ch = &_chc->at(i);

      float dt10 = ch->endTime(StrawEnd::hv)-ch->endTime(StrawEnd::cal);
      
      ComboHitHist_t* hr = Hist->ch[0];
      hr->h_dt10->Fill(dt10);
      hr->h_edep->Fill(ch->energyDep()*1.e3);

      if (ch->energyDep() > 0.0005) {
        ComboHitHist_t* hist = Hist->ch[1];
        hist->h_dt10->Fill(dt10);
        hist->h_edep->Fill(ch->energyDep()*1.e3);
      }
    }
  
//-----------------------------------------------------------------------------
// time cluster histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<_ntc; i++) {
      const TimeCluster* tc = &_tcc->at(i);
      float t0      = tc->t0().t0();
      int nh        = tc->nhits();
      int nhg       = 0;

      TimeClusterHist_t* hist = Hist->tc[0];
      
      for (int ih=0; ih<nh; ih++) {
        StrawHitIndex hit_index   = tc->hits().at(ih);
        const mu2e::ComboHit* hit = &_chc->at(hit_index);
        if (hit->energyDep() > _minEDep) {
          nhg++;
        }

        float dt = hit->correctedTime()-t0;
        hist->h_dt->Fill(dt);
        hist->h_edep->Fill(hit->energyDep()*1.e3);
      }
      
      hist->h_nh->Fill(nh);
      hist->h_nhg->Fill(nhg);
    }
                     
    return 0;
  }

//-----------------------------------------------------------------------------
  bool DaqanaFilter::beginRun(art::Run& ArtRun) {
    _rn = ArtRun.run();

    if (not _run_initialized) {
      if (_fillHistograms) {
                                        // make sure the job doesn't crash
        book_histograms(_rn);
      }
      _run_initialized = true;
    }
    return true;
  }
  
//-----------------------------------------------------------------------------
  bool DaqanaFilter::endRun(art::Run& run) {
    const float rate = (_nevt > 0) ? float(_nevp)/float(_nevt) : 0.f;
    TLOG(TLVL_DEBUG) << "passed:" << _nevp << " events out of:" << _nevt << " for a ratio of:" << rate;
    return true;
  }


  //-----------------------------------------------------------------------------
  bool DaqanaFilter::filter(art::Event& ArtEvent) {
    
    _art_event = &ArtEvent;         // should always be the first line
    
    if (_debugMode) print_(e_DEBUG,"-- START");
    
    ++_nevt;
//-----------------------------------------------------------------------------
// straw hits
//-----------------------------------------------------------------------------    
    art::Handle<mu2e::ComboHitCollection> chch;
    if (!ArtEvent.getByLabel(_chCollTag, chch)) {
      TLOG(TLVL_ERROR) << std::format("No straw hit collection tag:{}",_chCollTag.encode());
      _chc = nullptr;
    }
    else {
      _chc   =  chch.product();
      _ncht  = _chc->size();
    }
//-----------------------------------------------------------------------------
// time clusters
//-----------------------------------------------------------------------------    
    art::Handle<mu2e::TimeClusterCollection> tcch;
    if (!ArtEvent.getByLabel(_tcCollTag, tcch)) {
      TLOG(TLVL_ERROR) << std::format("No time cluster collection tag:{}",_tcCollTag.encode());
      _tcc = nullptr;
    }
    else {
      _tcc   =  tcch.product();
      _ntc   = _tcc->size();
    }
//-----------------------------------------------------------------------------
// process waveforms, count good hits
//-----------------------------------------------------------------------------
    _nchg          = 0;
    
    for (int i = 0; i<_ncht; ++i) {
      const mu2e::ComboHit* ch = &_chc->at(i);
      float dt = ch->endTime(StrawEnd::hv)-ch->endTime(StrawEnd::cal);
      if (fabs(dt)  > _maxDt  ) continue;
      if (ch->energyDep() < _minEDep) continue;
      _nchg++;
    }
//-----------------------------------------------------------------------------
// count number of good time clusters
// good time cluster: a TC with at least _minNGoodHits above the _minEDep threshold
// _minEDep threshold is the same for good hits and good hits in time clusters
//-----------------------------------------------------------------------------
    _ntcg = 0;
    _tc_mask.clear();
    
    for (int i=0; i<_ntc; i++) {
      const TimeCluster* tc = &_tcc->at(i);
      int n_good_ch = 0;
      int tc_bitmask = 0;
      int nch        = tc->nhits();
      
      for (int ih=0; ih<nch; ih++) {
        StrawHitIndex hit_index   = tc->hits().at(ih);
        const mu2e::ComboHit* hit = &_chc->at(hit_index);
        if (hit->energyDep() > _minEDep) {
          n_good_ch++;
        }
      }
      
      if (nch < 10) tc_bitmask |= 0x1;

      if (n_good_ch < _minNGoodCh) {
        tc_bitmask |= 0x02;
      }
      else {
        _ntcg++;
      }
     
      _tc_mask.push_back(tc_bitmask);
    }
//------------------------------------------------------------------------------
    _ncrvc = 0;
    art::Handle<mu2e::CrvCoincidenceClusterCollection>     crvcch;
    bool ok = ArtEvent.getByLabel(_crvcCollTag,crvcch);
    if (ok) { 
      _crvcc = crvcch.product();
      _ncrvc = _crvcc->size();
    }
    else {
      if (_requireCrvc) {
        print_(e_ERROR,std::format("CrvCoincidenceClusterCollection:{:s} not found\n",
                                   _crvcCollTag.encode().data()));
      }
    }
//-----------------------------------------------------------------------------
// matches between the CRVC and the tracker time clusters
//-----------------------------------------------------------------------------
    _n_tc_crvc_matches = 0;
    for (int itc=0; itc<_ntc; itc++) {
      const TimeCluster* tc = &_tcc->at(itc);
      if (_tc_mask[itc] != 0) continue;
      float tc_t0 = tc->t0().t0();
      
      for (int k=0; k<_ncrvc; k++) {
        const mu2e::CrvCoincidenceCluster* crvc = &_crvcc->at(k);

        float crvc_time = crvc->GetAvgHitTime();
        if (fabs(crvc_time-tc_t0) < 2000) {
          _n_tc_crvc_matches += 1;
        }
      }
    }
    
   
    if (_fillHistograms) fill_histograms(&_hist);
  
    if (_debugMode) print_(e_DEBUG,std::format("-- END, n good combo hits:{}",_nchg));
    
    bool passed = false;
    if ((_nchg              >= _minNGoodCh       ) and
        (_ntcg              >= _minNTimeClusters ) and
        (_n_tc_crvc_matches >= _minNTcCrvcMatches)     ) {
      
      passed = true;
      _nevp++;
      // if (_fillHistograms) fill_histograms(&_hist[1]);
    }
  
    return passed;
  }
}

DEFINE_ART_MODULE(mu2e::DaqanaFilter)
