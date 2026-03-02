///////////////////////////////////////////////////////////////////////////////
// change time cluster format to accomodata 18 stations
// the record is fairly sparse, but the number of time clusters is not too large
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

  int     _nhf  [18][4];
  float   _timef[18][4];

  int     _nhp  [18][2];
  float   _timep[18][2];

  int     _mnid      [18][12];
  int     _nh_panel  [18][12];
  float   _time_panel[18][12];
  float   _edep_panel[18][12];
  int     max_nh_panel;                 // max number of hits in one of the panels

  DaqTimeCluster();

  virtual ~DaqTimeCluster();

  int   mnid       (int stn, int ip) { return _mnid      [stn][ip]; }
  int   nhf        (int stn, int ip) { return _nhf       [stn][ip]; }
  float timef      (int stn, int ip) { return _timef     [stn][ip]; }
  int   nhp        (int stn, int ip) { return _nhp       [stn][ip]; }
  float timep      (int stn, int ip) { return _timep     [stn][ip]; }
  int   nh_panel   (int stn, int ip) { return _nh_panel  [stn][ip]; }
  float time_panel (int stn, int ip) { return _time_panel[stn][ip]; }
  float edep_panel (int stn, int ip) { return _edep_panel[stn][ip]; }

  ClassDef(DaqTimeCluster,1);
};

#endif
