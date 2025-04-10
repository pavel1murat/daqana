// ======================================================================
//
// ======================================================================
#ifndef __daqana_mod_CalFragmentAna_hh__
#define __daqana_mod_CalFragmentAna_hh__

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"

#include "Stntuple/mod/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

namespace mu2e {

  class CalFragmentAna : public THistModule {

    enum {
      kNChannels         = 96,
      kMaxDtcs           = 10,
      kMaxNLinks         =  6,
    };

    enum {
      kNEventHistSets    =  10,
      kNDtcHistSets      = 100,      // for now, just make it an array
      kNRocHistSets      = 100       // for now, just make it an array
    };

  public:

    struct WfParam_t {
      int   fs;     // first sample above _minPulseHeight
      float bl;     // baseline
      float ph;     // pulse height
      float q;      // Q(positive)
      float qt;     // Q(tail)
    };

    //-----------------------------------------------------------------------------  
    // Subevent header structure
    //-----------------------------------------------------------------------------
    struct SubEventHeader_t {
      uint32_t       byteCount   : 24;    
      uint16_t       unused      :  8;
      uint16_t       eventTag[3];
      uint16_t       numRocs     :  8;
      uint8_t        eventMode[5];
      // Packet #2
      uint8_t        dtcMacByte0;
      uint8_t        dtcMacByte2;
      uint8_t        evbMode;
      uint8_t        dtcID;
      uint8_t        subsystemID :  3;
      uint32_t       unused2     : 29;
      uint8_t        linkStatus[6];
      uint8_t        version;
      uint8_t        rf0Tdc;
      // Packet #3
      uint16_t       link4Lat;
      uint16_t       link5Lat;
      uint16_t       unused3[2];
      uint16_t       link0Lat;
      uint16_t       link1Lat;
      uint16_t       link2Lat;
      uint16_t       link3Lat;
      uint16_t       who_knows_what[12];
    };

    struct DtcDMAPacket_t {              // 2 16-byte words
      uint16_t       byteCount   : 16;
      uint8_t        subsystemID :  4;
      uint16_t       packetType  :  4;
      uint16_t       ROCID       :  3;
      uint16_t       unused      :  4;
      bool           valid       :  1;
    };

    struct RocDataHeaderPacket_t {            // 8 16-byte words in total
      // Word 0
      uint16_t            byteCount;
      // Word 1
      uint16_t            unused       : 4;
      uint16_t            packetType   : 4;
      uint16_t            linkID       : 3;
      uint16_t            DtcErrors    : 4;
      uint16_t            valid        : 1;
      // Word 2
      uint16_t            packetCount  : 11;
      uint16_t            unused2      : 2;
      uint16_t            subsystemID  : 3;
      // Words 3-5
      uint16_t            eventTag[3];
      // Word 6
      uint8_t             status       : 8;
      uint8_t             version      : 8;
      // Word 7
      uint8_t             dtcID        : 8;
      uint8_t             onSpill      : 1;
      uint8_t             subrun       : 2;
      uint8_t             eventMode    : 5;
      // Decoding status
      int                 empty()      { return (status & 0x01) == 0; }
      int                 invalid_dr() { return (status & 0x02); }
      int                 corrupt()    { return (status & 0x04); }
      int                 timeout()    { return (status & 0x08); }
      int                 overflow()   { return (status & 0x10); }
      int                 error_code() { return (status & 0x1e); }
    };

    struct DtcDataBlock_t : public RocDataHeaderPacket_t {
      uint16_t            hitData[10000];
    };

    struct EventHist_t {
      TH1F* nbtot;
      TH1F* nfrag;
      TH1F* nhits;
      TH1F* fsize;
      TH1F* n_empty;
      TH1F* n_invalid_dr;
      TH1F* n_corrupt;
      TH1F* n_timeouts;
      TH1F* n_overflows;
      TH1F* error_code;
      TH1F* nerr_tot;
      TH1F* valid;
      TH1F* eflg_vs_evt;
      TH1F* nerr_vs_evt;
      std::vector<TH1F*> fragment_size_per_link;
      TH1F* type_dist;
      TH1F* id_dist;
      TH1F* hits_per_channel;
      std::vector<TH1F*> channel_dist_per_board;
      TH1F* adc_values;
      TH1F* hitn_adc_value;      
    };

    struct Hist_t {
      EventHist_t*   event  [kNEventHistSets];
    } _hist;
    
    // pointer to the raw event data
    struct FragmentData_t {
      int nbytes;
    };

    struct RocData_t {
      int       link;
      int       size;
      int       nhits;
      int       nbytes;
      int       npackets;
      int       valid;
      int       dtc_id;
      
      int       n_empty;
      int       n_invalid_dr;
      int       n_corrupt;
      int       n_timeouts;
      int       n_overflows;
      int       error_code;
      int       nerr_tot;
      int       type;
      int       id;
      uint8_t   channel;
      std::vector<int> hit_channel;
      std::vector<uint16_t> footerData;
      std::vector<uint16_t> adc_values;
    };

    struct CalDtcData_t {
      int       dtcid;
      int       nbytes;
      int       nhits;
      int       error;
      RocData_t roc[6];
      
      int       n_empty;
      int       n_invalid_dr;
      int       n_corrupt;
      int       n_timeouts;
      int       n_overflows;
      int       type;
      int       id;
      uint8_t   channel;
    };

    struct EventData_t {
      const art::Event*  _event;
      int   nbtot;
      int   nhtot;
      int   nfrag;
      int   valid;
      
      int   n_nb_errors;     // wrong event size
      int   n_nwfs_errors;   // hit reported too many waveform samples
      int   n_linkid_errors; // hit reported wrong link ID
      int   n_chid_errors;   // hit reported wrong channel ID
      int   n_nchh_errors;   // too many hits in one channel
      
      int   n_empty;
      int   n_invalid_dr;
      int   n_corrupt;
      int   n_timeouts;
      int   n_overflows;
      int   error_code;
      int   nerr_tot;
     
      CalDtcData_t    caldtc[kMaxDtcs];
      std::vector<FragmentData_t> fragments;
    } _edata;

    // Data part: in reality, this is the fragment data.
    DtcDataBlock_t*  _calFragment;

    // Talk-to parameters
    int              _diagLevel;
    int              _minNBytes;
    int              _maxNBytes;
    int              _dataHeaderOffset;
    std::vector<int> _activeLinks_0;
    std::vector<int> _activeLinks_1;
    art::InputTag    _calCollTag;
    int              _analyzeFragments;
    int              _timeWindow;               // time window in ns
    int              _minNErrors;               // minimum number of errors for printout
    int              _errorCode;                // error code to print
    int              _validateAdcPatterns;
    int              _fillHistograms;           // <=0: don't fill
    int              _rocDataFormat;            // digis, patterns, etc
    int              _initialized;

    // Reference channels for different DTCs
    std::vector<int>* _activeLinks[2];         // active links - connected ROCs
    int              _nActiveLinks[2];
    
    double           _gen_offset[kNChannels];
    int              _node;

    std::unordered_map<std::string, TH1F*> _hitHistograms;
    std::vector< std::vector< std::vector< std::vector<art::TFileDirectory*> > > > _calDir;

    // NEW: Configuration flag to switch between old and new analysis methods.
    bool _useNewAnalysis;

    explicit CalFragmentAna(fhicl::ParameterSet const& pset);
    virtual ~CalFragmentAna() {}

    virtual void beginRun(const art::Run& ARun) override;
    virtual void beginJob() override;
    virtual void endJob() override;
    virtual void analyze(const art::Event& e) override;
    
    void book_calorimeter_directories();
    void analyze_dtc_fragment(const art::Event& e, const artdaq::Fragment* Fragment);
    void analyze_roc_data(RocDataHeaderPacket_t* Dh, RocData_t* Rd);
    void book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist);
    void book_histograms(int RunNumber);
    void debug(const art::Event& event);
    void fill_event_histograms(EventHist_t* Hist, EventData_t* Data);
    int fill_histograms();
    int init_event(const art::Event& AnEvent);
    void print_fragment(const artdaq::Fragment* Fragment, int NWords);
    void print_message(const char* Message);
  };
}

#endif
