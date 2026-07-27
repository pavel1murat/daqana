///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_RunDb_hh__
#define __daqana_obj_RunDb_hh__

class RunDb {
public:
                                        // run types
  enum {
    e_BEAM            = 1,
    e_COSMICS         = 2,
    e_INTERNAL_PULSER = 3,
    e_PULSE_INJECTION = 4,
    e_MONTE_CARLO     = 5,
  };
  
  struct Data_t {
    int   run_number;                   // 
    int   run_type;                     // 1:beam 2:cosmics 3:internal_pulser 4:pulse_injection 5:noise 100:MC
    int   event_mode;                   // on/off spill
    int   cfo_rate;                     // event rate, Hz
    int   trigger_rate;                 // event rate to disk, Hz
    float ew_length;                    // event window length, sec
    int   ref_channel;                  // for pulse injection runs
    int   pulsed_channel[12];           // assume one per 8
    int   plane_flag[36];               // 
  };

  static RunDb* fgInstance;

private:
  RunDb();
  ~RunDb();

public:
  static RunDb* Instance();

  int GetRunInfo(int RunNumber, RunDb::Data_t* RunInfo);

};

#endif
