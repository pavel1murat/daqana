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
  
  std::string fn = std::format("{}/daqana/rundb/{:06d}/{:06d}_runinfo.toml",getenv("SPACK_ENV"),subdir,subdir);
  toml::table tbl = toml::parse_file(fn);
    
  auto maps = tbl["runinfo"].as_array();
  //  int n_run_ranges = maps->size();

  for (auto&& node : *maps) {
    toml::table& tbl = *node.as_table();

    int rn = tbl["run_number"].value_or(-1);
    if (rn == RunNumber) {
      RunInfo->run_number   = RunNumber;                            // for kicks
      RunInfo->run_type     = tbl["run_type"     ].value_or(-1);
      RunInfo->event_mode   = tbl["event_mode"   ].value_or(1); // by default, on-spill 
      RunInfo->cfo_rate     = tbl["cfo_rate"     ].value_or(0.0); // to mark undefined
      RunInfo->ew_length    = tbl["ew_length"    ].value_or(0.0); // to mark undefined
      RunInfo->trigger_rate = tbl["trigger_rate" ].value_or(0.0); // to mark undefined

      if (RunInfo->run_type == e_PULSE_INJECTION) {
        // for pulse injection runs need to know pulsed channels, normally, one out of eight
        auto arr = tbl["pulsed_channels"].as_array();
        for (int i=0; i<12; i++) {
          auto v = arr->get_as<int64_t>(i);
          RunInfo->pulsed_channel[i] = v->value_or(-1);
        }
      }
      
      break;
    }
  }

  RunInfo->run_type    = 1;                // for kicks
  return rc;
}
