//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <format>
#include <deque>
#include "daqana/obj/CrvChannelMap_t.hh"

CrvChannelMap_t* CrvChannelMap_t::fgInstance(nullptr);

//-----------------------------------------------------------------------------
CrvChannelMap_t::CrvChannelMap_t(int RunNumber) {
  fRunNumber = RunNumber;
  // initialize
                                        // assume using spack
  
  std::string fn = std::format("{}/Offline/CRVConditions/data/extracted_v04.txt",getenv("SPACK_ENV"));
    
  std::ifstream input(fn);
  if (!input) {
    std::cerr << "Cannot open " << fn << '\n';
    return;
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
    if (!(stream >> dat.och >> dat.roc >> dat.feb >> dat.fch)) {
      if (line_number == 1)
        continue;
      
      std::cerr << "Invalid line " << line_number
                << ": " << line << '\n';
      return;
    }

    // if (!valid_indices(dat.roc,dat.feb,dat.fch)) {
    //   std::cerr << "Indices out of range on line "
    //             << line_number << ": " << line << '\n';
    //   return;
    // }
    
    _data.push_back(dat);

    Data_t* p = &_data.back();
    _ch_data_by_offline[dat.och]                   = p;
    _ch_data_by_online [dat.roc][dat.feb][dat.fch] = p;
  }

  std::cout << "Read " << _data.size() << " channel mappings\n";

}

//-----------------------------------------------------------------------------
CrvChannelMap_t::~CrvChannelMap_t() {
}

//-----------------------------------------------------------------------------
CrvChannelMap_t* CrvChannelMap_t::Instance(int RunNumber) {
  if  (fgInstance == nullptr) {
    fgInstance = new CrvChannelMap_t(RunNumber);
  }
  return fgInstance;
}
