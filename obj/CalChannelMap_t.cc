//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <format>
#include <deque>
#include <vector>
#include "daqana/obj/CalChannelMap_t.hh"

CalChannelMap_t* CalChannelMap_t::fgInstance(nullptr);

//-----------------------------------------------------------------------------
CalChannelMap_t::CalChannelMap_t() {
  _data.resize(kMaxNChannels);
}

//-----------------------------------------------------------------------------
int CalChannelMap_t::Init(int RunNumber) {
  fRunNumber = RunNumber;
  // initialize
                                        // assume using spack
  
  std::string fn = std::format("{}/Offline/CRVConditions/data/extracted_v04.txt",getenv("SPACK_ENV"));
    
  std::ifstream input(fn);
  if (!input) {
    std::cerr << "Cannot open " << fn << '\n';
    return -1;
  }

  std::string line;
  std::size_t line_number = 0;

  while (std::getline(input, line)) {
    ++line_number;
    
    if (line.empty() || line[0] == '#') continue;
    
    std::istringstream stream(line);
    
    Data_t dat;
    
    // This also skips the header:
    // Channel ROC FEB FEBchannel
    if (!(stream >> dat.cid >> dat.sipm >> dat.disk >> dat.crate >> dat.board)) {
      dat.channel = dat.cid*2+dat.sipm;
      if (line_number == 1)
        continue;
      
      std::cerr << "Invalid line " << line_number
                << ": " << line << '\n';
      return -2;
    }

    // if (!valid_indices(dat.roc,dat.feb,dat.fch)) {
    //   std::cerr << "Indices out of range on line "
    //             << line_number << ": " << line << '\n';
    //   return;
    // }
    
    _data[dat.channel] = dat;

    // Data_t* p = &_data.back();
    // _ch_data_by_offline[dat.och]                   = p;
    // _ch_data_by_online [dat.roc][dat.feb][dat.fch] = p;
  }

  std::cout << "Read " << _data.size() << " channel mappings\n";

  return 0;

}

//-----------------------------------------------------------------------------
CalChannelMap_t::~CalChannelMap_t() {
}

//-----------------------------------------------------------------------------
CalChannelMap_t* CalChannelMap_t::Instance() {
  if  (fgInstance == nullptr) {
    fgInstance = new CalChannelMap_t();
  }
  return fgInstance;
}
