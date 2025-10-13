#ifndef ComboHitData_t_hh
#define ComboHitData_t_hh

// format of the printout
struct ComboHitData_t {
                                        // read from the combohit printout
  int    i;
  int    nsh;
  int    sid;                           // parsing the sid gives the offline indices
  int    flags;
  double x;
  double y;
  double z;
  double phi;
  double time;
  double tcorr;
  double edep;
  double drtime;
  double prtime;
  double tres;
  double wdist;
  double wres;
  int    simid;
  double p;
  double pz;
  int    pdg;
  int    pdgm;
  int    genid;
                                        // so far, transients, undefined
  int    drs;                           // drift sign
  double r;                             // drift radius
  int    flag;                          //

  int    IsGood() { return 1; }
};

#endif
