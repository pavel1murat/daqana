#ifndef __TrackerDQM_module_hh__
#define __TrackerDQM_module_hh__

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic pop

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

      fhicl::Atom<int>             diagLevel              {Name("diagLevel"         )    , Comment("diag level"                 ) };
      fhicl::Atom<int>             minNBytes              {Name("minNBytes"         )    , Comment("min N(bytes)"               ) };
      fhicl::Atom<int>             maxNBytes              {Name("maxNBytes"         )    , Comment("max N(bytes)"               ) };
      fhicl::Atom<int>             dataHeaderOffset       {Name("dataHeaderOffset"  )    , Comment("data header offset"         ) };
      fhicl::Sequence<int>         activeLinks            {Name("activeLinks"       )    , Comment("active Links"               ) };
      fhicl::Sequence<int>         refChCal               {Name("refChCal"          )    , Comment("reference channel CAL side" ) };
      fhicl::Sequence<int>         refChHV                {Name("refChHV"           )    , Comment("reference channel HV  side" ) };
      fhicl::Atom<int>             dumpDTCRegisters       {Name("dumpDTCRegisters"  )    , Comment("1:dumpDTCRegisters"         ) };
      fhicl::Atom<int>             analyzeFragments       {Name("analyzeFragments"  )    , Comment("1:analyzeFragments"         ) };
      fhicl::Atom<int>             maxFragmentSize        {Name("maxFragmentSize"   )    , Comment("max fragment size"          ) };
      fhicl::Atom<int>             pulserFrequency        {Name("pulserFrequency"   )    , Comment("pulser frequency"           ) };
      fhicl::Atom<int>             port                   {Name("port"              )    , Comment("port"                       ) };
      fhicl::Atom<int>             timeWindow             {Name("timeWindow"        )    , Comment("time window, 25 ns ticks"   ) };
      fhicl::Atom<int>             nSamplesBL             {Name("nSamplesBL"        )    , Comment("N(baseline samples)"        ) };
      fhicl::Atom<float>           minPulseHeight         {Name("minPulseHeight"    )    , Comment("min pulse over threshold"   ) };
      fhicl::Sequence<int>         plotWaveforms          {Name("plotWaveforms"     )    , Comment("[link, channel]"            ) };
    };

                                        // TODO use constants from StrawID
  enum {
    kNStations          = 18,
    kNPlanesPerStation  =  2,
    kNPanelsPerPlane    =  6,
    kNChannels          = 96,
    kMaxNLinks          =  6,
    kMaxNHitsPerChannel = 10,
    kMaxNSamples        = 30
  };

  enum {
    kNBytesErrorBit     = 0x0001,
    kNWfsErrorBit       = 0x0002,
    kLinkIDErrorBit     = 0x0004,  // link ID error bit
    kChIDErrorBit       = 0x0008,
    kNChHitsErrorBit    = 0x0010,

    kNErrorBits         = 5
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
    TH1F*         dt0;                // T0 distance between the two consequtive pulses
    TH1F*         dt1;                // T1 distance between the two consequtive pulses
    TH1F*         dt2;                // T2 = (dt1+dt2)/2
    TH1F*         dt0r;               // T0(ich,0)-T0(ref,0)
    TH1F*         dt1r;               // T1(ich,0)-distance between the two pulses (if more than one)

    TH1F*         fsample;
    TH1F*         tmean;              // mean time
    TH1F*         bline;
    TH1F*         pheight;
    TH1F*         q;                  // waveform charge Q
    TH1F*         qt;                 // tail charge Qt
    TH1F*         qtq;                // Qt/Q

    TH1F*         raw_wf[kMaxNHitsPerChannel];
    TH1F*         wf    [kMaxNHitsPerChannel];
  };

//-----------------------------------------------------------------------------
// per-event histograms
//-----------------------------------------------------------------------------
  struct EventHist_t {
    TH1F*         nbtot;
    TH1F*         nfrag;
    TH1F*         nhits;
    TH1F*         fsize;
    TH1F*         error;
    TH1F*         error_rate;
    TH1F*         n_nb_errors;
    TH1F*         n_nwfs_errors;
    TH1F*         n_linkid_errors;
    TH1F*         n_chid_errors;
    TH1F*         n_nchh_errors;
    TH1F*         valid;
  };

//-----------------------------------------------------------------------------
// per-ROC histograms (or per-panel) histograms
//-----------------------------------------------------------------------------
  struct RocHist_t {
    TH1F*         nbytes;
    TH1F*         npackets;
    TH1F*         nhits;
    TH1F*         valid;
    TH2F*         nh_vs_ch;
    TH2F*         nh_vs_adc1;
    
    TH2F*         dt0r_vs_ch;
    TH2F*         dt1r_vs_ch;
                                        // time difference between the two reference channels,
                                        // each TDC separately
    TH1F*         dt0r01;
    TH1F*         dt1r01;

    TH1F*         sum_error_vs_ch;

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

//-----------------------------------------------------------------------------
// forgetting, for now, about multiple DTC's
//-----------------------------------------------------------------------------
  struct Hist_t {
    EventHist_t     event;
    RocHist_t       roc[kNStations][kNPlanesPerStation][kNPanelsPerPlane];
  };

  struct ChannelData_t {
    int      error;                    // 1:filled up
    int      nhits;
    float    dt0r;                     // time dist btw this channel and an FPGA reference channel, TDC0, ns
    float    dt1r;                     // time dist btw this channel and an FPGA reference channel, TDC1, ns
    float    dt0r_c;                   // the same, corrected for the FPGA-specific generator time offset
    float    dt1r_c;
    
    mu2e::TrackerDataDecoder::TrackerDataPacket* hit[kMaxNHitsPerChannel];
    WfParam_t  wp[kMaxNHitsPerChannel];
  };
  
  struct RocData_t {
    int       size;
    int       nhits;
    int       nbytes;
    int       npackets;
    int       valid;
    
    ChannelData_t  channel[kNChannels];
    ChannelData_t* ref_ch [2];
    
    float     dt0r01;                 // time difference between the two reference channels, TDC0, ns
    float     dt1r01;                 // time difference between the two reference channels, TDC1, ns
  };

  struct FragmentData_t {
    int       nbytes;
  };

  struct EventData_t {
    int       nbtot;                  // total nbytes
    int       nhtot;
    int       nfrag;
    int       valid;
    int       error;
    int       n_nb_errors;     // 0x01 : wrong event size
    int       n_nwfs_errors;   // 0x02 : hit reported too many wafeform samples
    int       n_linkid_errors; // 0x04 : hit reported wrond channel ID
    int       n_chid_errors;   // 0x08 : hit reported wrond channel ID
    int       n_nchh_errors;   // 0x10 : too many hits in one channel

    RocData_t rdata[kNStations][kNPlanesPerStation][kMaxNLinks];
    
    std::vector<FragmentData_t> fragments;
  } _event_data;

//-----------------------------------------------------------------------------
// in reality, this is the fragment data, an event can contain multiple fragments
//-----------------------------------------------------------------------------
  DtcDataBlock_t*  _trkFragment;
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
  int              _diagLevel;
  int              _minNBytes;
  int              _maxNBytes;
  int              _dataHeaderOffset;
  std::vector<int> _activeLinks;        // active links - connected ROCs
  std::vector<int> _refChCal;           // reference channel on CAL side FPGA
  std::vector<int> _refChHV;            // reference channel on HV  side FPGA
  art::InputTag    _trkfCollTag;
  int              _dumpDTCRegisters;
  int              _analyzeFragments;
  int              _maxFragmentSize;
  int              _pulserFrequency;    // in kHz, either 60 or 250
  
  int              _timeWindow;         // time window (spacing between the two EWMs for a given run) 
  int              _nSamplesBL;
  float            _minPulseHeight;
  int              _port;               // http://localhost:port serves histograms
  std::vector<int> _plotWaveforms;
//-----------------------------------------------------------------------------
// the rest
//-----------------------------------------------------------------------------
  int              _nActiveLinks;
  int              _referenceChannel[kMaxNLinks][2];
  int              _station;
  int              _plane;
    
  int              _adc_index_0 [kNChannels]; // seq num of the channel 'i' in the readout sequence
  int              _adc_index_1 [kNChannels]; // fixed map, seq num of the channel 'i' in the readout sequence
  double           _gen_offset  [kNChannels];

  double           _freq;           // generator frequency, defined by the run number
  double           _dt;             // expected distance between the two pulses
  double           _tdc_bin;        // 
  double           _tdc_bin_ns;     // TDC bin, in nanoseconds
  int              _initialized;    // histograms are booked in beginRun, protect ...

  Hist_t           _hist;
  PlotWaveform_t   _plot_wf;
  
  art::ServiceHandle<art::TFileService> tfs;

  TApplication*      _app;
  TCanvas*           _canvas[100];
  TBrowser*          _browser;
  const art::Event*  _event;
//
public:
  
  explicit     TrackerDQM(art::EDAnalyzer::Table<Config> const& conf);
  
  void         analyze (art::Event const& anEvent) override;
  void         beginJob()                          override;
  void         beginRun(art::Run   const& aRun   ) override;
  void         endJob  ()                          override;
  void         endRun  (art::Run   const& aRun   ) override;


  void         analyze_fragment       (const art::Event& e, const artdaq::Fragment* Fragment);

  void         book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int Ich);
  void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t*   Hist);
  void         book_roc_histograms    (art::TFileDirectory* Dir, int RunNumber, RocHist_t*     Hist, int Link);
  void         book_histograms        (int RunNumber);
  
  //  void         fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data);
  void         fill_event_histograms  (EventHist_t*   Hist, EventData_t*   Data);
  void         fill_roc_histograms    (RocHist_t*     Hist, RocData_t*     Data);

  // returns -1 if in trouble
  int          fill_histograms        ();

  // NWords: number of 2-byte words
  void         printFragment   (const artdaq::Fragment* Fragment, int NWords);
  void         unpack_adc_waveform(mu2e::TrackerDataDecoder::TrackerDataPacket* Hit, float* Wf, WfParam_t* Wp);

};

#endif
