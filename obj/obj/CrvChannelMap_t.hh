///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_mod_CrvChannelMap_t_hh__
#define __daqana_mod_CrvChannelMap_t_hh__

#include<deque>

class CrvChannelMap_t {
public:
  enum {
    kNRocs           = 18,
    kNFebs           = 24,
    kNChannelsFeb    = 64,
    kNChannelsMax    = kNRocs*kNFebs*kNChannelsFeb,
  };
  
  struct Data_t {
    int  och;                           // offline channel
    int  roc;
    int  feb;
    int  fch;                           // FEB channel
  };

  static CrvChannelMap_t* fgInstance;
                                        // there are 6000 smth bars < 25000 channels
  
  std::deque<Data_t> _data;
  
  int      fRunNumber;
  int      fNChannels;
  
  Data_t* _ch_data_by_offline   [kNChannelsMax];     // [mnid]
  Data_t* _ch_data_by_online    [kNRocs][kNFebs][kNChannelsFeb];  //

private:
  CrvChannelMap_t(int RunNumber);
  ~CrvChannelMap_t();

public:
  static CrvChannelMap_t* Instance(int RunNumber);

  Data_t* ch_data_by_online (int roc, int feb, int fch) { return _ch_data_by_online [roc][feb][fch]; }
  Data_t* ch_data_by_offline(int och)                   { return _ch_data_by_offline[och]; }

};

#endif
