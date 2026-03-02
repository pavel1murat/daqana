//

#include "daqana/obj/DaqTimeCluster.hh"

ClassImp(DaqTimeCluster)

//-----------------------------------------------------------------------------
DaqTimeCluster::DaqTimeCluster() : TObject(), nsh(-1), nch(-1), t0(-1.), tmin(-1), tmax(-1) {

  npanels = 0;
  nfaces  = 0;
  nplanes = 0;

  tmin     =  1.e12;
  tmax     = -1.e12;
  edep_max = -1.e12;
  chi2yz   = -1.e12;
  y0       = -1.e12;
  dydz     = -1.e12;

  for (int stn=0; stn<18; stn++) {
    for (int i=0; i<2; i++) {
      _nhp  [stn][i] = 0;
      _timep[stn][i] = 0.;
    }
    for (int i=0; i<4; i++) {
      _nhf  [stn][i] = 0;
      _timef[stn][i] = 0.;
    }

    for (int i=0; i<12; i++) {
      _nh_panel  [stn][i] = 0;
      _time_panel[stn][i] = 0.;
      _edep_panel[stn][i] = 0.;
      _mnid      [stn][i] = -1;
    }
  }
  max_nh_panel = 0;
}

//-----------------------------------------------------------------------------
DaqTimeCluster::~DaqTimeCluster() {
}
