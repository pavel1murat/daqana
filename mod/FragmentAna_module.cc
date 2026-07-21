// ======================================================================
// clang-format off
// PM
// FragmentAna:  DTC fragment-level analysis of the raw data - mostly printouts
// _debugMode > 0 : enables diagnostic printouts
// _debugBit[0]: print raw fragments
// _debugBit[1]: tracker digis
// _debugBit[2]: tracker waveforms
// ======================================================================
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/Decoders/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>
#include <format>

#include <regex>
#include <string>
#include <format>

#include <map>
#include <memory>

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

// #define TRACEMF_USE_VERBATIM 1

#include "TRACE/tracemf.h"
#define TRACE_NAME "FragmentAna"
                                         // message facility stream
#define MF_STREAM "FRAGMENT_ANA"


namespace mu2e {
  class FragmentAna;
}
// ======================================================================

class mu2e::FragmentAna : public art::EDAnalyzer {

public:

  enum {
    e_DEBUG   = 0,
    e_INFO    = 1,
    e_WARNING = 2,
    e_ERROR   = 3,
    e_SEVERE  = 4,
  };

  struct Config {
    fhicl::Atom<int>             debugMode       {fhicl::Name("debugMode"       ), fhicl::Comment("2026-04-25 PM:")};
    fhicl::Sequence<std::string> debugBits       {fhicl::Name("debugBits"       ), fhicl::Comment("debug bits"    )};
    fhicl::Atom<int>             analyzeFragments{fhicl::Name("analyzeFragments"), fhicl::Comment("2026-04-25 PM:")};
  };

  // --- C'tor/d'tor:
  explicit FragmentAna(const art::EDAnalyzer::Table<Config>& config);
  virtual ~FragmentAna() {}

  void         print_(int Level, const std::string&  Message,
                      const std::source_location& location = std::source_location::current());

  void         print_fragment(const artdaq::Fragment* Frag);

    // --- overloaded functions of the art producer
  virtual void analyze (const art::Event& ArtEvent) override;
  virtual void beginRun(const art::Run&   ArtRun  ) override;

  enum {kNDebugBits = 100};

private:
                                        // talk-to parameters
  int                      _debugMode ;
  std::vector<std::string> _debugBits;
  int                      _debugBit[kNDebugBits];
  int       _analyzeFragments;
                                        // the rest
  int       nADCPackets_{-1};           // N(ADC packets per hit)
  int       nSamples_   {-1};           // N(ADC samples per hit)
  int       np_per_hit_ {-1};           // N(data packets per hits)

  const art::Event*                 _event;
  int                               _last_run;

  ProditionsHandle<TrackerPanelMap> _tpm_h;
  const TrackerPanelMap*            _trackerPanelMap;
};

// ======================================================================
mu2e::FragmentAna::FragmentAna(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer   {config},
    _debugMode        (config().debugMode       ()),
    _debugBits        (config().debugBits       ()),
    _analyzeFragments (config().analyzeFragments()),
    _event            (nullptr)
{
//-----------------------------------------------------------------------------
// initialize debug bits  : debugBits: [ "bit0:1" , "bit14:1" ]
//-----------------------------------------------------------------------------
  for (int i=0; i<kNDebugBits; ++i) _debugBit[i] = 0;

  const char* key;
  int nbits = _debugBits.size(); // from FCL

  for (int i=0; i<nbits; i++) {
    int index(0), value(0);
    key               = _debugBits[i].data();
    sscanf(key,"bit%i:%i",&index,&value);
    _debugBit[index]  = value;

    print_(e_INFO,std::format("FragmentAna: bit={:4d} is set to {}",index,_debugBit[index]));
  }

  _last_run = -1;

}

//-----------------------------------------------------------------------------
// Level:     0: debug
//            1: info
//            2: warning
//            3: error
//------------------------------------------------------------------------------
void mu2e::FragmentAna::print_(int Level, const std::string& Message, const std::source_location& location) {

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
  if (_event) s = std::format("event: {}:{}:{} ",_event->run(),_event->subRun(),_event->event());

  std::vector<std::string> ss = xx.splitString(location.file_name(),"/");

  if (Level == e_DEBUG) {
    MF_LOG_TRACE("FRAGMENT_ANA") << s << ss.back() << ":" << location.line() << " : " << Message;
  } 
  else if (Level == e_INFO) {
    MF_LOG_VERBATIM("FRAGMENT_ANA") << s << ss.back() << ":" << location.line() << ":" << location.function_name() << ": " << Message;
  }
  else if (Level == e_WARNING) {
    MF_LOG_PRINT(MF_STREAM) << "WARNING: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }
  else if (Level == e_ERROR) {
    MF_LOG_PROBLEM(MF_STREAM) << "ERROR: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }
  else if (Level == e_SEVERE) {
    MF_LOG_ABSOLUTE(MF_STREAM) << "SEVERE: " << s << ss.back() << ":" << location.line() << " : " << Message;
  }
}

//-----------------------------------------------------------------------------
// HEX print of a fragment, the Mu2e data come in 2-byte words
//-----------------------------------------------------------------------------
void mu2e::FragmentAna::print_fragment(const artdaq::Fragment* Frag) {
  ushort* buf = (ushort*) (Frag->dataBegin());
  int nw      = Frag->dataSizeBytes()/2;
  int loc     = 0;

  for (int i=0; i<nw; i++) {
    if (loc == 0) printf(" 0x%08x: ",i*2);

    ushort  word = buf[i];
    printf("0x%04x ",word);

    loc += 1;
    if (loc == 8) {
      printf("\n");
      loc = 0;
    }
  }

  if (loc != 0) printf("\n");
}

//-----------------------------------------------------------------------------
void mu2e::FragmentAna::beginRun(const art::Run&  ArtRun) {
}

// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::FragmentAna::analyze(const art::Event& event) {
  int const packet_size(16); // in bytes

  if (_debugMode > 0) print_(e_INFO,"-- START");

  _event = &event;                      // cache to print events

  _trackerPanelMap = &_tpm_h.get(event.id());
//-----------------------------------------------------------------------------
// defined by the first hit
//-----------------------------------------------------------------------------
  artdaq::Fragments    fragments;
  artdaq::FragmentPtrs containerFragments;

  auto fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

  if (_debugMode > 0) {
    std::string msg = std::format("n_fragment_collections):{}",fragmentHandles.size());
    print_(e_INFO,msg);
  }

  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->empty())     continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      for (const auto& cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        for (size_t ii = 0; ii < contf.block_count(); ++ii) {
          containerFragments.push_back(contf[ii]);
          fragments.push_back(*containerFragments.back());
        }
      }
    }
    else {
//-----------------------------------------------------------------------------
// the 'handle' handles a list of artdaq fragments
// each artdaq fragment corresponds to a single DTC, or a plane
// loop over them
//-----------------------------------------------------------------------------
      int n_fragments = handle->size();
      
      if (_debugMode) {
        print_(e_INFO,std::format("-- next fragment collection,  n_fragments:{}",n_fragments));
      }
      
      for (int ifrag=0; ifrag<n_fragments; ifrag++) {
        const artdaq::Fragment* frag = &handle->at(ifrag);

        if (_debugMode and (_debugBit[0] > 0)) {
          print_(e_INFO,
                 std::format("-- fragment number:{} type:{} version:{} ID:{} timestamp:{} data_size:{} type:{} DTC_SubEventHeader.size:{}",
                             ifrag,frag->type(),frag->version(),frag->fragmentID(),frag->timestamp(),frag->dataSizeBytes(),
                             frag->typeString(),sizeof(DTCLib::DTC_SubEventHeader)));
          print_fragment(frag);
        }
        if (not _analyzeFragments) continue;

        if (frag->type() == mu2e::FragmentType::CFO) {
//-----------------------------------------------------------------------------
// skip CFO fragment (type = 12)
//-----------------------------------------------------------------------------
          continue;
        }
//-----------------------------------------------------------------------------
// skip fragments with the payload size less than the DTC header size
// do it only for the current data format (runs > 107236)
//-----------------------------------------------------------------------------
        if (frag->dataSizeBytes() <= sizeof(DTCLib::DTC_SubEventHeader)) {
          std::string msg = std::format("fragment:{} data size:{} < DTC_SubEventHeader.size:{}. SKIP FRAGMENT",
                                        ifrag,frag->dataSizeBytes(),sizeof(DTCLib::DTC_SubEventHeader));
          print_(e_ERROR,msg);
          continue;
        }
//-----------------------------------------------------------------------------
// skip data header
//-----------------------------------------------------------------------------
        uint8_t* fdata = (uint8_t*) (frag->dataBegin()) + sizeof(DTCLib::DTC_EventHeader);
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
            (seh->link5_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker)     ) {
          
          continue;
        }
        uint32_t dtc_id = seh->source_dtc_id;
//-----------------------------------------------------------------------------
// this is a tracker DTC fragment, loop over the ROCs
//-----------------------------------------------------------------------------
        ushort*  buf          = (ushort*) fdata;
        int      nbytes       = buf[0];             // frag.dataSizeBytes() includes extra 0x20
        uint8_t* roc_data     = fdata+sizeof(*seh);
        uint8_t* last_address = fdata+nbytes;

        while (roc_data < last_address) {
          RocDataHeaderPacket_t* rdh = (RocDataHeaderPacket_t*) roc_data;
          int nhits = 0;
          // int header_printed = 0;
//------------------------------------------------------------------------------
// skip empty ROC blocks
//------------------------------------------------------------------------------
          if (rdh->packetCount > 1) {
            if (nADCPackets_ < 0) {
//-----------------------------------------------------------------------------
// take the number of ADC packets per hit from the hit data,
// trust that but watch if it changes
// so far, any corruptions we saw were contained withing the ROC payload, and nhits
// was a reliable number
// 2026-03-20: don't store waveforms if only one packet per hit
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* h0;
              h0           = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+packet_size);
//-----------------------------------------------------------------------------
// don't expect zero ADC packets, the error seen was that one of the ROCs just
// doesn't send the second packet at all... work under that assumption
//-----------------------------------------------------------------------------
              if (h0->NumADCPackets == 0) {
                print_(e_ERROR,std::format("dtc_id:{} link_id:{} N(ADC packets) = 0, skip ROC data",
                                   dtc_id,(int) rdh->linkID));
                
                roc_data += (rdh->packetCount+1)*packet_size;
                continue;
              }
              nADCPackets_ = h0->NumADCPackets;
              nSamples_    = 3+12*nADCPackets_;
              np_per_hit_  = nADCPackets_+1;
            }
            uint32_t link_id = rdh->linkID;
            nhits            = rdh->packetCount/(nADCPackets_+1);
//-----------------------------------------------------------------------------
// there should not be more than 255 hits per ROC, if nhits>=255 it is a corruption,
// stop processing of the event 
//-----------------------------------------------------------------------------
            if (nhits > 255) {
              print_(e_ERROR,std::format("nhits:{}, skip event",nhits));
              break;
            }
            
            const TrkPanelMap::Row* tpmd = _trackerPanelMap->panel_map_by_online_ind(dtc_id,link_id);
            if (tpmd == nullptr) {
//-----------------------------------------------------------------------------
// either DTC ID or link ID are corrupted. Haven't seen that so far, switch to the next ROC anyway
//-----------------------------------------------------------------------------
              print_(e_ERROR,std::format("either dtc_id:{} or link_id:{} is corrupted, skip ROC data",
                                 dtc_id,link_id));

              roc_data += (nhits*np_per_hit_+1)*packet_size;
              continue;
            }

            if (_debugMode) {
              print_(e_INFO,std::format("-- DTC:{} ROC:{} nhits:{}",dtc_id,link_id,nhits));
            }

            for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* hit_data ;

              int offset = (ihit*np_per_hit_+1)*packet_size;   // in bytes
              hit_data   = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+offset);
              if (roc_data+offset >= last_address) {
                print_(e_ERROR,std::format("dtc_id:{} link_id:{} roc_data:{} offset:{} last_address:{} , SKIPPING",
                                   dtc_id, link_id, (void*) roc_data, offset, (void*) last_address));
                break;
              }
//-----------------------------------------------------------------------------
// at this point, check consistency between the channel_id, dtc_id and link_id for a given run
// panel ID is a derivative of the DTC ID and the link iD
// mn_id - 'MinnesotaID' of the panel
//-----------------------------------------------------------------------------
              mu2e::StrawDigiFlag digi_flag;
              uint16_t channel = static_cast<uint16_t>(hit_data->StrawIndex);
              uint16_t chid    = mu2e::StrawId(channel).straw(); // channel ID within the panel

              if (chid >= StrawId::_nstraws) {
                if (_debugBit[52] == 0) {
                  print_(e_ERROR,std::format("hit with corrupted chid:{:04x} : straw:{} / dtc_id:{} link_id:{}, SKIPPING",
                                     hit_data->StrawIndex, chid, dtc_id, link_id));
                }
                continue;
              }

              uint16_t mnid    = channel >> mu2e::StrawId::_panelsft;

              if (hit_data->NumADCPackets != nADCPackets_) {
                int np = hit_data->NumADCPackets;
                print_(e_ERROR,std::format("wrong NADCpackets:{} , expected:{}, GO TO THE NEXT ROC",
                                   np,nADCPackets_));
                break;
              }
            }
          }
//-----------------------------------------------------------------------------
// end fo ROC data processing, on to the next one
//-----------------------------------------------------------------------------
          roc_data += (nhits*np_per_hit_+1)*packet_size;
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// Store the straw digis in the event
//-----------------------------------------------------------------------------

  if (_debugMode) print_(e_INFO,"-- END");
}


// //-----------------------------------------------------------------------------
// uint16_t mu2e::FragmentAna::parse_minnesota_label(std::string label){
//     if ((label.size() != 5) || (label[0] != 'M') || (label[1] != 'N')){
//         std::string msg = "invalid minnesota label: " + label;
//         throw cet::exception("FragmentAna") << msg << std::endl;
//     }
//     std::string substr = label.substr(2, 3);
//     unsigned int parsed;
//     int scanned = sscanf(substr.c_str(), "%u", &parsed);
//     if (scanned != 1){
//       std::string msg = "failed to parse minnesota label: " + label;
//       throw cet::exception("FragmentAna") << msg << std::endl;
//     }
//     uint16_t rv = static_cast<uint16_t>(parsed);
//     return rv;
// }

// ======================================================================

DEFINE_ART_MODULE(mu2e::FragmentAna)

// ======================================================================
