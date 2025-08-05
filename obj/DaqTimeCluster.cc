//

#include "daqana/obj/DaqTimeCluster.hh"

ClassImp(DaqTimeCluster)

//-----------------------------------------------------------------------------
DaqTimeCluster::DaqTimeCluster() : TObject(), nsh(-1), nch(-1), t0(-1.), tmin(-1), tmax(-1) {

  npanels = 0;
  nfaces  = 0;
  nplanes = 0;

  tmin    =  1.e12;
  tmax    = -1.e12;
  
  for (int i=0; i<2; i++) {
    _nhp  [i] = 0;
    _timep[i] = 0.;
  }
  for (int i=0; i<4; i++) {
    _nhf  [i] = 0;
    _timef[i] = 0.;
  }

  for (int i=0; i<12; i++) {
    _nh_panel  [i] = 0;
    _time_panel[i] = 0.;
    _mnid      [i] = -1;
  }
}

//-----------------------------------------------------------------------------
DaqTimeCluster::~DaqTimeCluster() {
}
