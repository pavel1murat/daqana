#ifndef __daqana_obj_fit_functions_hh__
#define __daqana_obj_fit_functions_hh__

// from fit_dt_01_msh.cc

extern int fit_dt01_msh    (int Slot, int FirstRun, int Panel1=0, int Panel2=6, int FirstChannel=0,int LastChannel=95, int PrintLevel=1);

extern int fit_dt01_cosmics(int RunNumber, int Panel1=0, int Panel2=12, int FirstChannel=0, int LastChannel=95, int PrintLevel=0);

#endif
