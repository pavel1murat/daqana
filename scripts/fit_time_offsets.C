///////////////////////////////////////////////////////////////////////////////
// fit histograms produced by Stntuple/mod/TrkFragmentAna_module.cc
//-----------------------------------------------------------------------------

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

//-----------------------------------------------------------------------------
// assume the hist files are stored in the default directory
//-----------------------------------------------------------------------------
int init(int* Index) {
  int adc_index_1[96] = {
    91, 85, 79, 73, 67, 61, 55, 49,
    43, 37, 31, 25, 19, 13,  7,  1,
    90, 84, 78, 72, 66, 60, 54, 48,
      
    42, 36, 30, 24, 18, 12,  6,  0,
    93, 87, 81, 75, 69, 63, 57, 51,
    45, 39, 33, 27, 21, 15,  9,  3,
      
    94, 88, 82, 76, 70, 64, 58, 52,
    46, 40, 34, 28, 22, 16, 10,  4,
    95, 89, 83, 77, 71, 65, 59, 53,
      
    44, 38, 32, 26, 20, 14,  8,  2, 
    92, 86, 80, 74, 68, 62, 56, 50,
    47, 41, 35, 29, 23, 17, 11,  5
  };


  for (int i=0; i<96; i++) {
    Index[i] = adc_index_1[i];
  }

  // _referenceChannel[0] = 91;
  // _referenceChannel[1] = 94;

  return 0;
}


int fit_time_offsets(int RunNumber) {

  char fn[200];

  const char* hist_dir = "/srv/mu2e/data/projects/tracker/vst/hist";
  sprintf(fn,"%s/trkvst.annex.trk_fragment_ana.%06i.hist",hist_dir,RunNumber);

  TFile* f = TFile::Open(fn);

  int adc_index[96];

  init(adc_index);

  TH1D* hpy[96];

  struct fit_result_t {
    double p[3];
    double e[3];
    double chi2dof;
  };

  fit_result_t fr[96]; 

  TH2F* h0 = (TH2F*) f->Get("//TrkFragmentAna/trk/frag_0/dt0rc_vs_adc_0");

  for (int i=0; i<48; i++) {
    printf("-- channel %02i\n",i);
    int bin  = i+1;
    int ich  = adc_index[i];

    hpy[ich] = h0->ProjectionY(Form("hpy_%02i",i),bin,bin);
//-----------------------------------------------------------------------------
// fit and retrieve the fit results
//-----------------------------------------------------------------------------
    if (i == 0) {
      fr[i].chi2dof = -1;

      for (int ip=0; ip<3; ip++) {
        fr[i].p[ip] = 0.;
        fr[i].e[ip] = 0.;
      }          
      continue;
    }
    TFitResultPtr tfr = hpy[ich]->Fit("gaus","sq","");

    if ((! tfr->IsValid()) or tfr->IsEmpty()) {
      printf("# FIT ERROR: index: %2i channel: %2i \n",i,ich);
      continue;
    }

    fr[i].chi2dof = tfr->Chi2()/tfr->Ndf();
    double sf     = sqrt(fr[ich].chi2dof);

    for (int ip=0; ip<3; ip++) {
      fr[i].p[ip] = tfr->Parameter(ip);
      fr[i].e[ip] = tfr->Error(ip)*sf;
    }          
  }

//-----------------------------------------------------------------------------
// the second FPGA
//-----------------------------------------------------------------------------
  TH2F* h1 = (TH2F*) f->Get("//TrkFragmentAna/trk/frag_0/dt0rc_vs_adc_1");

  for (int i=48; i<96; i++) {
    int bin  = i+1;
    int ich  = adc_index[i];
    printf("-- channel %02i\n",i);

    hpy[ich] = h1->ProjectionY(Form("hpy_%02i",i),bin,bin);
//-----------------------------------------------------------------------------
// fit and retrieve the fit results
//-----------------------------------------------------------------------------
    if (i == 48) {
      fr[i].chi2dof = -1;

      for (int ip=0; ip<3; ip++) {
        fr[i].p[ip] = 0.;
        fr[i].e[ip] = 0.;
      }          
      continue;
    }

    TFitResultPtr tfr = hpy[ich]->Fit("gaus","sq","");

    if ((! tfr->IsValid()) or tfr->IsEmpty()) {
      printf("# FIT ERROR: index: %2i channel: %2i \n",i,ich);
      continue;
    }

    fr[i].chi2dof = tfr->Chi2()/tfr->Ndf();
    double sf     = sqrt(fr[i].chi2dof);

    for (int ip=0; ip<3; ip++) {
      fr[i].p[ip] = tfr->Parameter(ip);
      fr[i].e[ip] = tfr->Error(ip)*sf;
    }          
  }
//-----------------------------------------------------------------------------
// print fir resutls
//-----------------------------------------------------------------------------
  printf(" i ich       mean     emean     sig       esig       chi2dof\n");
  printf("------------------------------------------------------------\n");
  for (int i=0; i<96; i++) {
    int ich = adc_index[i];
    printf(" %2i %2i",i,ich);
    printf(" %12.5f %12.5f  %12.5f %12.5f  %12.5f\n",
           fr[i].p[1],fr[i].e[1], fr[i].p[2], fr[i].e[2], fr[i].chi2dof);
  }


  return 0;
}
