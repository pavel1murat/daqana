#ifndef __daqana_mod_WfParam_t__
#define __daqana_mod_WfParam_t__

struct WfParam_t {
  int   fs;                         // first sample above _minPulseHeight
  float bl;                         // baseline
  float ph;                         // pulse height
  float q;                          // Q(positive)
  float qt;                         // Q(tail)
};

#endif
