//
#include <iostream>
#include <string>
#include <format>
#include "toml++/toml.hpp"
#include "daqana/obj/RunDb.hh"

RunDb* RunDb::fgInstance(nullptr);

//-----------------------------------------------------------------------------
RunDb::RunDb() {
  // initialize
}

//-----------------------------------------------------------------------------
RunDb::~RunDb() {
}

//-----------------------------------------------------------------------------
RunDb* RunDb::Instance() {
  static RunDb* run_db(nullptr);

  if  (run_db == nullptr) {
    run_db     = new RunDb();
    fgInstance = run_db;
  }
  return fgInstance;
}

//-----------------------------------------------------------------------------
int RunDb::GetRunInfo(int RunNumber, RunDb::Data_t* RunInfo) {
  int rc(0);

  RunInfo->run_number = -1;
  
  int subdir = (RunNumber / 1000 )*1000;
  
  std::string fn = std::format("{}/daqana/rundb/{:06d}/{:06d}_rundb.toml",getenv("SPACK_ENV"),subdir,subdir);
  toml::table tbl = toml::parse_file(fn);
    
  auto maps = tbl["rundata"].as_array();
  //  int n_run_ranges = maps->size();

  for (auto&& node : *maps) {
    toml::table& tbl = *node.as_table();

    int rn = tbl["run"].value_or(-1);
    if (rn == RunNumber) {
      RunInfo->run_number  = RunNumber;        // for kicks
      auto arr = tbl["pulsed_channels"].as_array();
      for (int i=0; i<12; i++) {
        auto v = arr->get_as<int64_t>(i);
        RunInfo->pulsed_channel[i] = v->value_or(-1);
      }
      
      break;
    }
  }

  RunInfo->run_type    = 1;                // for kicks
  return rc;
}
