#ifndef __daqana_ana_fit_dt01_hh__
#define __daqana_ana_fit_dt01_hh__
///////////////////////////////////////////////////////////////////////////////
// fit dt01 histograms 
// 
// TrkPreampStraw table sign convention: stored are deltaT01 = T(cal)-T(hv) ,
// so better plot histogrrams for those differences, to avoid inverting the sign
// when producing the constants.
//
// deltaT01 residuals do not account for the signal propagation time
// for 1m straw, that could be of the order of 4 ns, better to correct for that
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

class fit_dt01 {
public:
  
  struct fit_result_t {
    int    ich;                           // defined on input
    double p[3];
    double e[3];
    double chi2dof;
    float  n0;
    float  ntot;
    float  ineff;
  };

  int fRunNumber;

  fit_result_t fFr[216][96]; // fit results 

  int fit_histogram(TH1F* Hist, fit_result_t* Frr);
  int fit_msh      (int FirstRun, int Panel1=0, int Panel2=12, int PrintLevel=1);
  int fit_pin      (const char* Fn, int* RunNumber, int SlotLo, int SlotHi, int Panel, int PrintLevel=1);
  int fit_cosmics  (int  RunNumber, int Slot, int Panel, int Channel, const char* Fn = nullptr, int PrintLevel=1);

  // need run number to determine whether the calibrations have to be corrected for the signal proparation
  int write_TrkPreampStraw(int RunNumber, const char* Fn = "TrkPreampStraw.txt.new");
};

#endif
