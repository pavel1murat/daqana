///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_mod_TrkPanelMap_t_hh__
#define __daqana_mod_TrkPanelMap_t_hh__

class TrkPanelMap_t {
public:
  struct Data_t {
    int  mnid;                        // 101='MN101' etc
    int  dtc_id;
    int  link;
    int  plane;                       // geo index of the plane
    int  ppid;                        // production plane   ID  (plane_21 ... )
    int  panel;                       // "geo" panel index
    int  zface;                       // z-ordered face (so far, random, could've calculated, *TODO*)
  };

  static TrkPanelMap_t* fgInstance;

  Data_t   _data[500];

  int      fRunNumber;
  Data_t* _panel_data_by_mnid   [500];    // [mnid]
  Data_t* _panel_data_by_online [36][6];  // [dtc_id][link]
  Data_t* _panel_data_by_offline[36][6];  // [geo_plane][geo_panel]

private:
  TrkPanelMap_t(int RunNumber);
  ~TrkPanelMap_t();

public:
  static TrkPanelMap_t* Instance(int RunNumber);

  Data_t* panel_data_by_mnid   (int mnid)              { return _panel_data_by_mnid   [mnid];          };
  Data_t* panel_data_by_online (int dtc_id, int link ) { return _panel_data_by_online [dtc_id][link ]; };
  Data_t* panel_data_by_offline(int plane , int panel) { return _panel_data_by_offline[plane ][panel]; };

};

#endif
