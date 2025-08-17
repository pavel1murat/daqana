///////////////////////////////////////////////////////////////////////////////
// for now, assume the time cluster in one station, update later
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_obj_DaqTimeCluster_hh__
#define __daqana_obj_DaqTimeCluster_hh__

#include "TObject.h"

class DaqTimeCluster : public TObject {
public:
  int     nsh ;
  int     nch ;
  float   t0  ;
  float   tmin;
  float   tmax;
  float   edep_max;

  float   y0;
  float   dydz;
  float   chi2yz;                       // per DOF

  int     nplanes;
  int     nfaces ;
  int     npanels;   // n panels with hits

  int     _nhf  [4];
  float   _timef[4];

  int     _nhp  [2];
  float   _timep[2];

  int     _mnid      [12];
  int     _nh_panel  [12];
  float   _time_panel[12];
  float   _edep_panel[12];

  DaqTimeCluster();

  virtual ~DaqTimeCluster();

  int   mnid       (int i) { return _mnid [i]; }
  int   nhf        (int i) { return _nhf  [i]; }
  float timef      (int i) { return _timef[i]; }
  int   nhp        (int i) { return _nhp  [i]; }
  float timep      (int i) { return _timep[i]; }
  int   nh_panel   (int i) { return _nh_panel  [i]; }
  float time_panel (int i) { return _time_panel[i]; }
  float edep_panel (int i) { return _edep_panel[i]; }

  ClassDef(DaqTimeCluster,1);
};

#endif
