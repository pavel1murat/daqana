// ======================================================================
//
// ======================================================================
#ifndef __daqana_mod_TrkFragmentAna_hh__
#define __daqana_mod_TrkFragmentAna_hh__

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

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

  class TrkFragmentAna : public THistModule {

    enum { kNChannels         = 96,
           kMaxNLinks         =  6,
           kMaxNHWfPerChannel = 10,
           kMaxStations       = 2       // for now, just make it an array
    };

    enum {
      kNBytesErrorBit     = 0x0001,
      kNWfsErrorBit       = 0x0002,
      kLinkIDErrorBit     = 0x0004,  // link ID error bit
      kChIDErrorBit       = 0x0008,
      kNChHitsErrorBit    = 0x0010,
      
      kNErrorBits         = 5
    };

  public:

    struct WfParam_t {
      int   fs;                         // first sample above _minPulseHeight
      float bl;                         // baseline
      float ph;                         // pulse height
      float q;                          // Q(positive)
      float qt;                         // Q(tail)
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
                                        // per-ROC histograms
    struct EventHist_t {
      TH1F*         nbtot;
      TH1F*         nfrag;
      TH1F*         nhits;
      TH1F*         fsize;
      TH1F*         error;
      TH1F*         n_nb_errors;
      TH1F*         n_nwfs_errors;
      TH1F*         n_linkid_errors;
      TH1F*         n_chid_errors;
      TH1F*         n_nchh_errors;
      TH1F*         valid;
    };

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
// multiple DTC's
//-----------------------------------------------------------------------------
    struct DtcHist_t {
      RocHist_t     roc[kMaxNLinks];
    };

    struct StationHist_t {
      DtcHist_t     dtc[2];
    };

    struct Hist_t {
      EventHist_t   event;
      StationHist_t stn[2];
    } _Hist;

    struct ChannelData_t {
      int      nhits;
      float    dt0r;                     // time dist btw this channel and an FPGA reference channel, TDC0, ns
      float    dt1r;                     // time dist btw this channel and an FPGA reference channel, TDC1, ns
      float    dt0r_c;                   // the same, corrected for the FPGA-specific generator time offset
      float    dt1r_c;
      
      std::vector<TrackerDataDecoder::TrackerDataPacket*> hit;
      std::vector<WfParam_t>                              wp; 
    };

    struct RocData_t {
      int       link;
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

    struct StationData_t {
      int        dtcid[2];
      int        nbytes[2];
      int        nhits[2];
      int        error[2];
      RocData_t  roc[2][6];
    };
                                        // pointer to the raw event data
    struct FragmentData_t {
      int       nbytes;
    };

    struct EventData_t {
      int           nbtot;                  // total nbytes
      int           nhtot;
      int           nfrag;
      int           valid;
      int           error;
      int           n_nb_errors;     // 0x01 : wrong event size
      int           n_nwfs_errors;   // 0x02 : hit reported too many wafeform samples
      int           n_linkid_errors; // 0x04 : hit reported wrond channel ID
      int           n_chid_errors;   // 0x08 : hit reported wrond channel ID
      int           n_nchh_errors;   // 0x10 : too many hits in one channel
      StationData_t station[kMaxStations];

      std::vector<FragmentData_t> fragments;
    } _edata;

    // struct Config {
    //   fhicl::Atom<int>               diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    //   fhicl::Atom<int>               parseTRK {fhicl::Name("parseTRK" ), fhicl::Comment("parseTRK"        )};
    //   fhicl::Atom<art::InputTag>     trkTag   {fhicl::Name("trkTag"   ), fhicl::Comment("trkTag"          )};
    //   fhicl::Table<TAnaDump::Config> tAnaDump {fhicl::Name("tAnaDump" ), fhichl::Comment("Diag plugin"    )};
    // };
//-----------------------------------------------------------------------------
// data part
//-----------------------------------------------------------------------------
                                        // in reality, this is the fragment data, 
                                        // an event can contain multiple fragments
    DtcDataBlock_t*  _trkFragment;
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    int              _diagLevel;
    int              _minNBytes;
    int              _maxNBytes;
    int              _dataHeaderOffset;
    std::vector<int>  _activeLinks_0;
    std::vector<int>  _activeLinks_1;
    std::vector<int>*_activeLinks[2];         // active links - connected ROCs
    std::vector<int> _refChCal;               // reference channel on CAL side FPGA
    std::vector<int> _refChHV;                // reference channel on HV  side FPGA
    art::InputTag    _trkfCollTag;
    int              _dumpDTCRegisters;
    int              _analyzeFragments;
    int              _maxFragmentSize;
    int              _pulserFrequency;          // in kHz, either 60 or 250

    int              _timeWindow;               // time window (spacing between the two EWMs for a given run)
    int              _nADCPackets;              // number of waveform packets
    int              _nSamplesBL;               // number of first samples used to determine the baseline
    float            _minPulseHeight;           // threshold for the charge integration;
    int              _nStations;
//-----------------------------------------------------------------------------
// the rest, use the same reference channels for different DTCs - the ROC FW is the same
//-----------------------------------------------------------------------------
    int              _nActiveLinks[2];
    int              _referenceChannel[kMaxNLinks][2];
    
    int              _adc_index_0 [kNChannels]; // seq num of the channel 'i' in the readout sequence
    int              _adc_index_1 [kNChannels]; // fixed map, seq num of the channel 'i' in the readout sequence
    double           _gen_offset  [kNChannels];

    double           _freq;           // generator frequency, defined by the run number
    double           _dt;             // expected distance between the two pulses
    double           _tdc_bin;        // 
    double           _tdc_bin_ns;     // TDC bin, in nanoseconds
    int              _initialized;    // histograms are booked in beginRun, protect ...
    int              _idtc;           // fragment number, today - proxy to the DTC ID
    int              _station;

    const art::Event*  _event;

    explicit TrkFragmentAna(fhicl::ParameterSet const& pset);
    // explicit TrkFragmentAna(const art::EDAnalyzer::Table<Config>& config);
    virtual ~TrkFragmentAna() {}
    
    virtual void beginRun(const art::Run& ARun) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    virtual void analyze         (const art::Event& e) override;
    void         analyze_fragment(const art::Event& e, const artdaq::Fragment* Fragment);

    void         book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int Link, int Ich);

    void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t*   Hist);
    
    void         book_dtc_histograms    (art::TFileDirectory* Dir, int RunNumber, DtcHist_t*     Hist,
                                         int IStation, int IDtc);
    void         book_roc_histograms    (art::TFileDirectory* Dir, int RunNumber, RocHist_t*     Hist,
                                         int IStation, int IDtc, int Link);
    
    void         book_histograms        (int RunNumber);
  
    void         fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data);
    void         fill_dtc_histograms    (DtcHist_t*     Hist, StationData_t* Data, int IDtc);
    void         fill_event_histograms  (EventHist_t*   Hist, EventData_t*   Data);
    void         fill_roc_histograms    (RocHist_t*     Hist, RocData_t*     Data);

                                        // returns -1 if in trouble
    int          fill_histograms        ();

    // NWords: number of 2-byte words
    void         printFragment      (const artdaq::Fragment* Fragment, int NWords);
    int          unpack_adc_waveform(TrackerDataDecoder::TrackerDataPacket* Hit, float* Wf, WfParam_t* Wp);
  };
}
#endif
