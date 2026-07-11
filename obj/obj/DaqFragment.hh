///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqFragment_hh__
#define __daqana_obj_DaqFragment_hh__

#include "TObject.h"
                                        // for ntuple
struct RocData_t {
  uint32_t  nbytes;                     // forgot if that includes the ROC header
  uint64_t  ewtag;
  uint16_t  packet_type;
  uint16_t  link;
  uint16_t  err;
  uint16_t  valid;
  uint16_t  npackets;
  uint8_t   ssid;                    // subsystem ID
  uint16_t  status;                     //
  uint16_t  version;                     //
  uint16_t  dtc_id;
  uint16_t  onspill;
  uint16_t  subrun;
  uint16_t  event_mode;
};



class DaqFragment : public TObject {
public:
  uint32_t  nbytes;
  uint64_t  ewtag;
  uint16_t  nrocs;
  uint64_t  event_mode;
  uint16_t  mac_addr;
  uint16_t  partition;
  uint16_t  evb_mode;
  uint16_t  dtc_id;
  uint32_t  link_ssid[6];                  // 3 bits per subsystem
  uint8_t   link_status[6];
  uint8_t   version;
  uint8_t   emtdc;
  uint16_t  latency[6];
  
  RocData_t roc[6];

  DaqFragment();

  virtual ~DaqFragment();

  //  ClassDefOverride(DaqFragment,1);
};

#endif
