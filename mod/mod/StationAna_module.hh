// -*- buffer-read-only:t -*- 
// ======================================================================
//
// ======================================================================
#ifndef __daqana_mod_StationAna_hh__
#define __daqana_mod_StationAna_hh__

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

#ifndef __CLING__
#include "fhiclcpp/ParameterSet.h"
#else
namespace fhiclcpp {
  class ParameterSet;
};
#endif

// #include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
// #include "artdaq-core/Data/Fragment.hh"

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

  class StationAna : public THistModule {

    enum {
      kNStraws           = 96,
      kMaxNLinks         =  6,
      kMaxNHWfPerChannel = 10,
      kMaxStations       =  2       // for now, just make it an array
    };

    enum {
      kNEventHistSets =  10,
      kNPanelHistSets =  10         // for now, just make it an array
    };

  public:
                                        // per-ROC histograms
    struct EventHist_t {
      TH1F* evt;
      TH1F* nsht;
      TH1F* nshg;
      TH1F* nshdt;
      TH1F* nsht_vs_evt;
      TH1F* max_edep;
    };

    struct StrawHist_t {
      TH1F* nsht;
      TH1F* dtch;                       // ch : CAL-HV
      TH1F* edep;
      TH1F* tcal;
      TH1F* nshg;
      TH1F* dtchg;                       // ch : CAL-HV
      TH1F* edepg;
      TH1F* tcalg;
    };

    struct PanelHist_t {
      TH1F*         nsht;
      TH1F*         nsht_vs_evt;
      TH1F*         occup;
      TH1F*         tcal;
      TH1F*         dtch;
      TH1F*         edep;
      TH1F*         nshg;                // n good hits
      TH1F*         occupg;
      TH1F*         tcalg;
      TH1F*         dtchg;
      TH1F*         edepg;
      StrawHist_t   straw[96];
    };

    struct PanelHistSet_t {
      PanelHist_t   panel[12];
    };

    struct Hist_t {
      EventHist_t*      event    [kNEventHistSets];
      PanelHistSet_t    panel_set[kNPanelHistSets];
    };
    
    struct StrawData_t {
      std::vector<const mu2e::StrawHit*> list_of_hits;
      int  nsht() { return list_of_hits.size(); }
      
      std::vector<const mu2e::StrawHit*> list_of_good_hits;
      int  nshg() { return list_of_good_hits.size(); }
    };

    struct PanelData_t {
      int         mnid;
      int         nsht;
      int         nshg;
      StrawData_t straw_data[96];
    };
                                        // pointer to the raw event data
    struct EventData_t {
      const art::Event*  event;
      int                nsht;
      int                nshg;
      int                nshdt;
      float              max_edep;
      int                run_number;
      int                srn_number;
      int                evt_number;
      PanelData_t        panel[12];
    } _edata;

//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    int              _debugMode;
    int              _debugBit[100];
    art::InputTag    _shCollTag;
    float            _maxDt;
    float            _minEDep;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    const mu2e::StrawHitCollection* _shc;
    int                             _n_straw_hits;
    Hist_t                          _hist;
    int                             _initialized;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    explicit StationAna(fhicl::ParameterSet const& pset);
    // explicit StationAna(const art::EDAnalyzer::Table<Config>& config);
    virtual ~StationAna() {}
    
    virtual void analyze (const art::Event& ArtEvent) override;
    virtual void beginRun(const art::Run&   ArtRun  ) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist);
    void         book_panel_histograms  (art::TFileDirectory* Dir, int RunNumber, PanelHist_t* Hist, int Mnid);
    void         book_straw_histograms  (art::TFileDirectory* Dir, int RunNumber, StrawHist_t* Hist, int Mnid, int Is);

    void         book_histograms        (int RunNumber);
    void         debug                  (const art::Event& event);
  
    void         fill_straw_histograms  (StrawHist_t* Hist, StrawData_t* Data);
    void         fill_panel_histograms  (PanelHist_t* Hist, PanelData_t* Data);
    void         fill_event_histograms  (EventHist_t* Hist, EventData_t* Data);

    int          fill_histograms();

    int          getData(const art::Event& ArtEvent);

    int          init_event       (const art::Event& AnEvent);
    void         print_(const std::string& Message,
                        const std::source_location& location = std::source_location::current());
  };
}
#endif
