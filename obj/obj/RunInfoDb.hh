///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_RunInfoDb_hh__
#define __daqana_obj_RunInfoDb_hh__

class RunInfoDb {
public:
  struct Data_t {
    int  run_number;                    // 
    int  run_type;                      //
    int  ref_channel;                   // for pulse injection runs
    int  pulsed_channels[12];
    int  plane_flag[36];
  };

  static RunInfoDb* fgInstance;

private:
  RunInfoDb();
  ~RunInfoDb();

public:
  static RunInfoDb* Instance();

  int GetRunInfo(int RunNumber, RunInfoDb::Data_t* RunInfo);

};

#endif
