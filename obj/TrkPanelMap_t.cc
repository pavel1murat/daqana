//
#include <iostream>
#include <string>
#include <format>
#include "toml++/toml.hpp"
#include "daqana/obj/TrkPanelMap_t.hh"

TrkPanelMap_t* TrkPanelMap_t::fgInstance(nullptr);

//-----------------------------------------------------------------------------
TrkPanelMap_t::TrkPanelMap_t(int RunNumber) {
  fRunNumber = RunNumber;
    // initialize
                                        // assume using spack
  std::string fn = std::format("{}/daqana/calibrations/TrkPanelMap.toml",getenv("SPACK_ENV"));
  toml::table tbl = toml::parse_file(fn);
    
  auto maps = tbl["TrkPanelMap"].as_array();
  //  int n_run_ranges = maps->size();

  for (auto&& node : *maps) {
    toml::table& entry_table = *node.as_table();

    auto* range = entry_table["run_range"].as_array();

    int min_run = range->at(0).as_integer()->get();
    int max_run = range->at(1).as_integer()->get();

    if ((RunNumber >= min_run) and (RunNumber <= max_run)) {
      // found run range - panel data - array of records

      auto panel_data_array = entry_table["panel_data"].as_array();
      for (auto&& item : *panel_data_array) {
        toml::table& p = *item.as_table();
        
        int mnid   = p["mnid" ].value_or(-1);
        Data_t* r = &_data[mnid];
        r->mnid = mnid;
        r->dtc_id = p["dtc_id"].value_or(-1);
        r->link   = p["link"  ].value_or(-1);
        r->plane  = p["plane" ].value_or(-1);
        r->ppid   = p["ppid"  ].value_or(-1);
        r->panel  = p["panel" ].value_or(-1);
        r->zface  = p["zface" ].value_or(-1);

        _panel_data_by_mnid[mnid] = r;
        _panel_data_by_online[r->dtc_id][r->link] = r;
        _panel_data_by_offline[r->plane][r->panel] = r;
      }
    }
  }
}

//-----------------------------------------------------------------------------
TrkPanelMap_t::~TrkPanelMap_t() {
}

//-----------------------------------------------------------------------------
TrkPanelMap_t* TrkPanelMap_t::Instance(int RunNumber) {
  if  (fgInstance == nullptr) {
    fgInstance = new TrkPanelMap_t(RunNumber);
  }
  return fgInstance;
}
