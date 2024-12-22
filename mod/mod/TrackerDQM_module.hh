#ifndef __TrackerDQM_module_hh__
#define __TrackerDQM_module_hh__

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic pop

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "fhiclcpp/types/OptionalAtom.h"

#include "TBufferFile.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TBrowser.h"

// #include "Offline/RecoDataProducts/inc/StrawDigi.hh"
// #include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
// #include "Offline/DataProducts/inc/TrkTypes.hh"

#include "artdaq-core/Data/Fragment.hh"

#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
//-----------------------------------------------------------------------------
class TrackerDQM : public art::EDAnalyzer {
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<int>             diagLevel              {Name("diagLevel"          )    , Comment("diag level"                 ) };
      fhicl::Atom<art::InputTag>   trkfCollTag            {Name("trkfCollTag"        )    , Comment("track fragment coll tag"    ) };
      fhicl::Atom<int>             minNBytes              {Name("minNBytes"          )    , Comment("min N(bytes)"               ) };
      fhicl::Atom<int>             maxNBytes              {Name("maxNBytes"          )    , Comment("max N(bytes)"               ) };
      fhicl::Sequence<int>         activeLinks_0          {Name("activeLinks_0"      )    , Comment("DTC0 active Links"          ) };
      fhicl::Sequence<int>         activeLinks_1          {Name("activeLinks_1"      )    , Comment("DTC1 active Links"          ) };
      fhicl::Sequence<int>         refChCal               {Name("refChCal"           )    , Comment("reference channel CAL side" ) };
      fhicl::Sequence<int>         refChHV                {Name("refChHV"            )    , Comment("reference channel HV  side" ) };
      fhicl::Atom<int>             dumpDTCRegisters       {Name("dumpDTCRegisters"   )    , Comment("1:dumpDTCRegisters"         ) };
      fhicl::Atom<int>             analyzeFragments       {Name("analyzeFragments"   )    , Comment("1:analyzeFragments"         ) };
      fhicl::Atom<int>             maxFragmentSize        {Name("maxFragmentSize"    )    , Comment("max fragment size"          ) };
      fhicl::Atom<int>             pulserFrequency        {Name("pulserFrequency"    )    , Comment("pulser frequency"           ) };
      fhicl::Atom<int>             nADCPackets            {Name("nADCPackets"        )    , Comment("N(ADC packets/hit)"         ) };
      fhicl::Atom<int>             port                   {Name("port"               )    , Comment("port"                       ) };
      fhicl::Atom<int>             timeWindow             {Name("timeWindow"         )    , Comment("time window, 25 ns ticks"   ) };
      fhicl::Atom<int>             nSamplesBL             {Name("nSamplesBL"         )    , Comment("N(baseline samples)"        ) };
      fhicl::Atom<float>           minPulseHeight         {Name("minPulseHeight"     )    , Comment("min pulse over threshold"   ) };
      fhicl::Atom<float>           minNErrors             {Name("minNErrors"         )    , Comment("min N(errors) to print"     ) };
      fhicl::Atom<float>           errorCode              {Name("errorCode"          )    , Comment("error code to print"        ) };
      fhicl::Atom<float>           validateADCPatterns    {Name("validateADCPatterns")    , Comment("validate ADC patterns"      ) };
      fhicl::Atom<float>           fillHistograms         {Name("fillHistograms"     )    , Comment("1:fill histograms"          ) };
      fhicl::Sequence<int>         plotWaveforms          {Name("plotWaveforms"      )    , Comment("[link, channel]"            ) };

      fhicl::Sequence<std::string> debugBits              {Name("debugBits"          )    , Comment("debug bits"                 ) };
    };

                                        // TODO use constants from StrawID
  enum {
    kNStations         = 18,
    kNPlanesPerStation =  2,
    kNPanelsPerPlane   =  6,
    kNChannels         = 96,
    kMaxNLinks         =  6,
    kMaxNHWfPerChannel = 10,
    kMaxNSamples       = 30
  };

  enum {
    kNEventHistSets    =  10,
    kNStationHistSets  = 100       // for now, just make it an array
  };

  enum {
    kNBytesErrorBit     = 0x0001,
    kNWfsErrorBit       = 0x0002,
    kLinkIDErrorBit     = 0x0004,  // link ID error bit
    kChIDErrorBit       = 0x0008,
    kNChHitsErrorBit    = 0x0010,
    kHitErrorBit        = 0x0100,  // hit error reported by the digi FPGA
    kAdcPatternErrorBit = 0x0200,  // wrong ADC pattern
      
    kNErrorBits         = 7
  };

  enum {
    kRocPattern1        = 1,
    kDigiCheckerBoard   = 4
  };

  struct WfParam_t {
    int   fs;                           // first sample above _minPulseHeight
    float bl;                           // baseline
    float ph;                           // pulse height
    float q;                            // Q(positive)
    float qt;                           // Q(tail)
    float q_x_i;
    float ns;                           // nsamples in the charge integration
    float tm;                           // q_x_i/ns
  };

  struct PlotWaveform_t {
    int station;
    int link;
    int channel;
  };

//-----------------------------------------------------------------------------
// subevent header is 3 packets
//-----------------------------------------------------------------------------
  struct SubEventHeader_t {
    // packet #1 
    uint32_t       byteCount   : 24;    
    uint16_t       unused      :  8;
    uint16_t       eventTag[3]     ;
    uint16_t       numRocs     :  8;
    uint8_t        eventMode[5]    ;
    // packet #2
    uint8_t        dtcMacByte0     ;
    uint8_t        dtcMacByte2     ;
    uint8_t        evbMode         ;
    uint8_t        dtcID           ;
    uint8_t        subsystemID :  3;
    uint32_t       unused2     : 29;
    uint8_t        linkStatus[6]   ;
    uint8_t        version         ;
    uint8_t        rf0Tdc          ;
    // packet #3
    uint16_t       link4Lat        ;
    uint16_t       link5Lat        ;
    uint16_t       unused3[2]      ;
    uint16_t       link0Lat        ;
    uint16_t       link1Lat        ;
    uint16_t       link2Lat        ;
    uint16_t       link3Lat        ;
  };

  struct DtcDMAPacket_t {              // 2 16-byte words
    uint16_t       byteCount   : 16;   ///< Byte count of current block
    uint8_t        subsystemID :  4;   ///< Hop count
    uint16_t       packetType  :  4;   ///< Packet type               DTC_PacketType
    uint16_t       ROCID       :  3;   ///< Link identifier of packet DTC_LinkID
    uint16_t       unused      :  4;
    bool           valid       :  1;   ///< Whether the DTC believes the packet to be valid
  };

  struct DtcDataHeaderPacket_t : public DtcDMAPacket_t {  // 8 16-byte words in total
    uint16_t            nPackets     : 11;
    uint16_t            unused       :  5;
    uint16_t            eventTag[3];             // DTC_EventWindowTag
    uint8_t             status;                  // DTC_DataStatus
    uint8_t             version;
    uint8_t             DTCID;
    uint8_t             EVBMode;
  };

  struct DtcDataBlock_t : public DtcDataHeaderPacket_t {
    uint16_t            hitData[10000];
  };

  struct RocDataHeaderPacket_t {        // 8 16-byte words in total
                                        // 16-bit word 0
    uint16_t            byteCount    : 16;
                                        // 16-bit word 1
    uint16_t            unused       : 4;
    uint16_t            packetType   : 4;
    uint16_t            linkID       : 3;
    uint16_t            DtcErrors    : 4;
    uint16_t            valid        : 1;
                                        // 16-bit word 2
    uint16_t            packetCount  : 11;
    uint16_t            unused2      : 2;
    uint16_t            subsystemID  : 3;
                                        // 16-bit words 3-5
    uint16_t            eventTag[3];
                                        // 16-bit word 6
    uint8_t             status       : 8;
    uint8_t             version      : 8;
                                        // 16-bit word 7
    uint8_t             dtcID        : 8;
    uint8_t             onSpill      : 1;
    uint8_t             subrun       : 2;
    uint8_t             eventMode    : 5;
                                        // decoding status
      
    int                 empty     () { return (status & 0x01) == 0; }
    int                 invalid_dr() { return (status & 0x02); }
    int                 corrupt   () { return (status & 0x04); }
    int                 timeout   () { return (status & 0x08); }
    int                 overflow  () { return (status & 0x10); }
      
    int                 error_code() { return (status & 0x1e); }
  };

//-----------------------------------------------------------------------------
// per-channel histograms
//-----------------------------------------------------------------------------
  struct ChannelHist_t {
    TH1F*         nhits;
    TH1F*         time[2];
    TH1F*         t0  [2];            // early times in ns
    TH1F*         t1  [2];            // late  times in ns
    TH1F*         tot [2];
    TH1F*         pmp;
    TH1F*         dt01[2];            // T0-T1 for each hit, with different binning
    TH1F*         dt0;                // T0 distance between the two consequtive pulses
    TH1F*         dt1;                // T1 distance between the two consequtive pulses
    TH1F*         dt2;                // T2 = (dt1+dt2)/2
    TH1F*         dt0r;               // T0(ich,0)-T0(ref,0)
    TH1F*         dt1r;               // T1(ich,0)-distance between the two pulses (if more than one)

    TH1F*         fsample;
    TH1F*         bline;
    TH1F*         pheight;
    TH1F*         q;                  // waveform charge Q
    TH1F*         qt;                 // tail charge Qt
    TH1F*         qtq;                // Qt/Q
    
    TH1F*         raw_wf[kMaxNHWfPerChannel];
    TH1F*         wf    [kMaxNHWfPerChannel];
  };
//-----------------------------------------------------------------------------
// per-event histograms
//-----------------------------------------------------------------------------
  struct EventHist_t {
    TH1F*         nbtot;
    TH1F*         nfrag;
    TH1F*         nhits;
    TH1F*         fsize;
    TH1F*         n_nb_errors;
    TH1F*         n_nwfs_errors;
    TH1F*         n_linkid_errors;
    TH1F*         n_chid_errors;
    TH1F*         n_nchh_errors;
    TH1F*         valid;

    TH1F*         error_code;
    TH1F*         nerr_tot;

    TH1F*         n_empty;
    TH1F*         n_invalid_dr;
    TH1F*         n_corrupt;
    TH1F*         n_timeouts;
    TH1F*         n_overflows;

    TH1F*         nerr_vs_evt;
    TH1F*         eflg_vs_evt;
  };
//-----------------------------------------------------------------------------
// per-ROC histograms (or per-panel) histograms
//-----------------------------------------------------------------------------
  struct RocHist_t {
    TH1F*         nbytes;
    TH1F*         npackets;
    TH1F*         nhits;
    TH1F*         valid;
    TH1F*         error_code;  // ROC error code

    TH1F*         n_empty;
    TH1F*         n_invalid_dr;
    TH1F*         n_corrupt;
    TH1F*         n_timeouts;
    TH1F*         n_overflows;

    TH1F*         nerr_tot;
    TH1F*         nerr_vs_evt;
    TH1F*         eflg_vs_evt;

    TH2F*         nh_vs_ch;
    TH2F*         nh_vs_adc0;

    TH2F*         dt0r_vs_ch;
    TH2F*         dt1r_vs_ch;
                                        // time difference between the two reference channels,
                                        // each TDC separately
    TH1F*         dt0r01;
    TH1F*         dt1r01;

    TH2F*         dt0rc_vs_ch[2];
    TH2F*         dt1rc_vs_ch[2];

    TH2F*         dt0rc_vs_adc[2];
    TH2F*         dt1rc_vs_adc[2];

    TH1F*         nhits_vs_ich;
    TH1F*         nhits_vs_adc[2];

    TProfile*     fs_vs_ich;
    TProfile*     bl_vs_ich;
    TProfile*     ph_vs_ich;
    TProfile*     q_vs_ich;
    TProfile*     qt_vs_ich;
    TProfile*     qtq_vs_ich;

    ChannelHist_t channel[kNChannels];
  };

//-----------------------------------------------------------------------------
// forgetting, for now, about multiple DTC's
//-----------------------------------------------------------------------------
  struct DtcHist_t {
    RocHist_t     roc[kMaxNLinks];
  };

  struct StationHist_t {
    DtcHist_t     dtc[2];
  };

  struct Hist_t {
    EventHist_t*   event  [kNEventHistSets  ];
    StationHist_t* station[kNStationHistSets];
  };
    
  struct ChannelData_t {
    int      error;                    // 1:filled up
    float    dt0r;                     // time dist btw this channel and an FPGA reference channel, TDC0, ns
    float    dt1r;                     // time dist btw this channel and an FPGA reference channel, TDC1, ns
    float    dt0r_c;                   // the same, corrected for the FPGA-specific generator time offset
    float    dt1r_c;
    
    std::vector<mu2e::TrackerDataDecoder::TrackerDataPacket*> hit;
    std::vector<WfParam_t>                              wp;

    int nhits() { return hit.size(); }
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
    int       error_code;             // anything but not-empty a byte..

    int       nerr_tot;
      
    ChannelData_t  channel[kNChannels];
    ChannelData_t* ref_ch [2];

    float     dt0r01;                 // time difference between the two reference channels, TDC0, ns
    float     dt1r01;                 // time difference between the two reference channels, TDC1, ns
  };

  struct StationData_t {
    int       dtcid[2];
    int       nbytes[2];
    int       nhits[2];
    int       error[2];
    RocData_t roc[2][6];
      
    int       n_empty;
    int       n_invalid_dr;
    int       n_corrupt;
    int       n_timeouts;
    int       n_overflows;
  };
                                        // pointer to the raw event data
  struct FragmentData_t {
    int       nbytes;
  };

  struct EventData_t {
    const art::Event*  _event;

    int           nbtot;                  // total nbytes
    int           nhtot;
    int           nfrag;
    int           valid;
      
    int           n_nb_errors;     // 0x01 : wrong event size
    int           n_nwfs_errors;   // 0x02 : hit reported too many wafeform samples
    int           n_linkid_errors; // 0x04 : hit reported wrond channel ID
    int           n_chid_errors;   // 0x08 : hit reported wrond channel ID
    int           n_nchh_errors;   // 0x10 : too many hits in one channel
      
    int           n_empty;
    int           n_invalid_dr;
    int           n_corrupt;
    int           n_timeouts;
    int           n_overflows;
    int           error_code;
    int           nerr_tot;

    StationData_t station[kNStations];

    std::vector<FragmentData_t> fragments;
  } _edata;

//-----------------------------------------------------------------------------
// in reality, this is the fragment data, an event can contain multiple fragments
//-----------------------------------------------------------------------------
  DtcDataBlock_t*  _trkFragment;
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
  int              _diagLevel;
  art::InputTag    _trkfCollTag;
  int              _minNBytes;
  int              _maxNBytes;
  //  int              _dataHeaderOffset;
  std::vector<int>  _activeLinks_0;
  std::vector<int>  _activeLinks_1;
  std::vector<int>*_activeLinks[2];         // active links - connected ROCs
  std::vector<int> _refChCal;           // reference channel on CAL side FPGA
  std::vector<int> _refChHV;            // reference channel on HV  side FPGA
  int              _dumpDTCRegisters;
  int              _analyzeFragments;
  int              _maxFragmentSize;
  int              _pulserFrequency;    // in kHz, either 60 or 250
  
  int              _timeWindow;         // time window (spacing between the two EWMs for a given run) 
  int              _nADCPackets;        // number of waveform packets
  int              _nSamplesBL;
  float            _minPulseHeight;
  int              _nStations;
  int              _minNErrors;               // min number of errros for printout on bit3
  int              _errorCode;                // errorCode to print
  int              _validateADCPatterns;      //
  int              _fillHistograms;           // <=0 : don't
  int              _fillWaveformHistograms;   // <=0 : don't
  int              _rocDataFormat;            // digis, patterns, etc

  std::vector<std::string> _debugBits;
  
  int              _debugBit[100];

  int              _port;                     // http://localhost:port serves histograms
  std::vector<int> _plotWaveforms;
//-----------------------------------------------------------------------------
// the rest
//-----------------------------------------------------------------------------
  int              _nActiveLinks[2];
  int              _referenceChannel[kMaxNLinks][2];
    
  int              _adc_index_0 [kNChannels]; // seq num of the channel 'i' in the readout sequence
  //  int              _adc_index_1 [kNChannels]; // fixed map, seq num of the channel 'i' in the readout sequence
  double           _gen_offset  [kNChannels];

  double           _freq;           // generator frequency, defined by the run number
  double           _dt;             // expected distance between the two pulses
  double           _tdc_bin;        // 
  double           _tdc_bin_ns;     // TDC bin, in nanoseconds
  int              _initialized;    // histograms are booked in beginRun, protect ...

  int              _station;
  int              _plane;

  Hist_t           _hist;
  PlotWaveform_t   _plot_wf;
  
  art::ServiceHandle<art::TFileService> tfs;

  TApplication*      _app;
  TCanvas*           _canvas[100];
  TBrowser*          _browser;
//
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
public:
  
  explicit     TrackerDQM(art::EDAnalyzer::Table<Config> const& conf);
  
  void         analyze (art::Event const& anEvent) override;
  void         beginJob()                          override;
  void         beginRun(art::Run   const& aRun   ) override;
  void         endJob  ()                          override;
  void         endRun  (art::Run   const& aRun   ) override;


  void         analyze_dtc_fragment(const art::Event& e, const artdaq::Fragment* Fragment);
  void         analyze_roc_data    (RocDataHeaderPacket_t* Dh, RocData_t* Rd);
  void         analyze_roc_patterns(RocDataHeaderPacket_t* Dh, RocData_t* Rd);

  void         book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist,
                                       int Link, int Ich);
  void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t*   Hist);
  void         book_dtc_histograms    (art::TFileDirectory* Dir, int RunNumber, DtcHist_t*     Hist,
                                       int IStation, int IDtc);
  void         book_roc_histograms    (art::TFileDirectory* Dir, int RunNumber, RocHist_t*     Hist,
                                       int IStation, int IDtc, int Link);
    
  void         book_histograms        (int RunNumber);
  void         debug                  (const art::Event& event);
  int          DebugBit(int I)        {return _debugBit[I]; }
  int          dtcIndex               (int DtcID);
  
  void         fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data);
  void         fill_dtc_histograms    (DtcHist_t*     Hist, StationData_t* Data, int IDtc);
  void         fill_event_histograms  (EventHist_t*   Hist, EventData_t*   Data);
  void         fill_roc_histograms    (RocHist_t*     Hist, RocData_t*     Data);
  int          fill_station_histograms(StationHist_t* Hist, EventData_t*   Data);

  // returns -1 if in trouble
  int          fill_histograms        ();

  int          init_event             (const art::Event& AnEvent);
    
  // NWords: number of 2-byte words
  void         print_fragment     (const artdaq::Fragment* Fragment, int NWords);
  void         print_hit          (const mu2e::TrackerDataDecoder::TrackerDataPacket* Hit);
  void         print_message      (const char* Message);

  void         unpack_adc_waveform(mu2e::TrackerDataDecoder::TrackerDataPacket* Hit, float* Wf, WfParam_t* Wp);

};

#endif
