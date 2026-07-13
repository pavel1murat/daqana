///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_RunDb_hh__
#define __daqana_obj_RunDb_hh__

class RunDb {
public:
  struct Data_t {
    int  run_number;                    // 
    int  run_type;                      //
    int  ref_channel;                   // for pulse injection runs
    int  pulsed_channel[12];
    int  plane_flag[36];
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
