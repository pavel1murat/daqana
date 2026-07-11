//
#include <iostream>
#include <string>
#include <format>
#include "toml++/toml.hpp"
#include "daqana/obj/RunInfoDb.hh"

RunInfoDb* RunInfoDb::fgInstance(nullptr);

//-----------------------------------------------------------------------------
RunInfoDb::RunInfoDb() {
    // initialize
}

//-----------------------------------------------------------------------------
RunInfoDb::~RunInfoDb() {
}

//-----------------------------------------------------------------------------
RunInfoDb* RunInfoDb::Instance() {
  static RunInfoDb* run_db(nullptr);

  if  (run_db == nullptr) {
    run_db     = new RunInfoDb();
    fgInstance = run_db;
  }
  return fgInstance;
}

//-----------------------------------------------------------------------------
int RunInfoDb::GetRunInfo(int RunNumber, RunInfoDb::Data_t* RunInfo) {
  int rc(0);

  RunInfo->run_number  = RunNumber;        // for kicks
  RunInfo->run_type    = 1;                // for kicks
  RunInfo->ref_channel = 21;               // OK for run=122655

  for (int i=0; i<36; i++) {
    RunInfo->plane_flag[i] = 1;
  }

  RunInfo->plane_flag[ 0] = 0;
  RunInfo->plane_flag[ 1] = 0;

  RunInfo->plane_flag[ 8] = 0;
  RunInfo->plane_flag[ 9] = 0;
  
  RunInfo->plane_flag[10] = 0;
  RunInfo->plane_flag[11] = 0;
  
  RunInfo->plane_flag[12] = 0;
  RunInfo->plane_flag[13] = 0;
  
  RunInfo->plane_flag[14] = 0;
  RunInfo->plane_flag[15] = 0;
  
  RunInfo->plane_flag[16] = 0;
  RunInfo->plane_flag[17] = 0;
  
                                        // assume using spack
//  std::string fn = std::format("{}/daqana/calibrations/TrkPanelMap.toml",getenv("SPACK_ENV"));
//  toml::table tbl = toml::parse_file(fn);
//    
//  auto maps = tbl["rundata"].as_array();
//  //  int n_run_ranges = maps->size();
//
//  for (auto&& node : *maps) {
//    toml::table& entry_table = *node.as_table();
//
//    auto* range = entry_table["run_range"].as_array();
//
//    int min_run = range->at(0).as_integer()->get();
//    int max_run = range->at(1).as_integer()->get();
//
//    if ((RunNumber >= min_run) and (RunNumber <= max_run)) {
//      // found run range - panel data - array of records
//
//      auto panel_data_array = entry_table["panel_data"].as_array();
//      for (auto&& item : *panel_data_array) {
//        toml::table& p = *item.as_table();
//        
//        int mnid   = p["mnid" ].value_or(-1);
//        Data_t* r = &_data[mnid];
//        r->mnid = mnid;
//        r->dtc_id = p["dtc_id"].value_or(-1);
//        r->link   = p["link"  ].value_or(-1);
//        r->plane  = p["plane" ].value_or(-1);
//        r->ppid   = p["ppid"  ].value_or(-1);
//        r->panel  = p["panel" ].value_or(-1);
//        r->zface  = p["zface" ].value_or(-1);
//
//        _panel_data_by_mnid[mnid] = r;
//        _panel_data_by_online[r->dtc_id][r->link] = r;
//        _panel_data_by_offline[r->plane][r->panel] = r;
//      }
//    }
//  }
//
  return rc;
}
