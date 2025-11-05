///////////////////////////////////////////////////////////////////////////////
// fit drho histograms to produce T0 channel-to-channel offsets
//-----------------------------------------------------------------------------

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

//-----------------------------------------------------------------------------
// example DsID: vst04s0s10r000
//-----------------------------------------------------------------------------
int fit_drho_vs_rdrift(const char* DsID, int RunNumber,int Panel) {

  char fn[200];

  // const char* hist_dir = "/srv/mu2e/data/projects/tracker/vst/hist";
  // const char* hist_dir = "./hist";
  const char* hist_dir = "/projects/mu2e/data/projects/vst/hist";
  sprintf(fn,"%s/hst.murat.%s.daqana.%06i_000001.root",hist_dir,DsID,RunNumber);

  printf("histfile:%s\n",fn);

  TFile* f = TFile::Open(fn);

  struct fit_result_t {
    double p[3];
    double e[3];
    double chi2dof;
  };

  fit_result_t fr[200];
  TH1D* hpy[200];

  TH2F* h0 = (TH2F*) f->Get(Form("//digis/psh_%i/drho_vs_r",Panel));

  int nx = h0->GetNbinsX();


  nx = 70; // force 3.5 mm
  
  for (int i=0; i<nx; i++) {
    // printf("-- straw %02i\n",i);
    int bin  = i+1;

    hpy[i] = h0->ProjectionY(Form("hpy_%02i",i),bin,bin);
//-----------------------------------------------------------------------------
// initialize and skip empty channels
//-----------------------------------------------------------------------------
    fr[i].chi2dof = -1;

    for (int ip=0; ip<3; ip++) {
      fr[i].p[ip] = 0.;
      fr[i].e[ip] = 0.;
    }        
    int nent = hpy[i]->GetEntries();
    if (nent<500) continue;
//-----------------------------------------------------------------------------
// fit and retrieve the fit results
//-----------------------------------------------------------------------------
    TFitResultPtr tfr = hpy[i]->Fit("gaus","sqw","");

    if ((! tfr->IsValid()) or tfr->IsEmpty()) {
      printf("# FIT ERROR: straw: %2i\n",i);
      continue;
    }
//-----------------------------------------------------------------------------
// if the fit was successfull, redo with +/ 600 um from the maximum
//-----------------------------------------------------------------------------
    double xmn =  tfr->Parameter(1);
    double err =  tfr->Error(1);

    printf("i nent xmn err : %3i %5i %10.4f %10.f\n",i,nent, xmn,err);
    
    tfr = hpy[i]->Fit("gaus","sq","",xmn-0.5,xmn+0.5);
    
    fr[i].chi2dof = tfr->Chi2()/tfr->Ndf();
    double sf     = sqrt(fr[i].chi2dof);

    for (int ip=0; ip<3; ip++) {
      fr[i].p[ip] = tfr->Parameter(ip);
      fr[i].e[ip] = tfr->Error(ip)*sf;
    }          
  }
//-----------------------------------------------------------------------------
// print fit resutlts with offsets in ns
// assume the drift velocity of 61 um/ns
//-----------------------------------------------------------------------------
  printf("panel ibin      rdrift       mean        emean          sig           esig       chi2dof\n");
  printf("#-----------------------------------------------------------------------------\n");
  double vdrift(61.);
  
  for (int i=0; i<nx; i++) {
    printf(" %4i  %4i %10.3f" ,Panel,i,h0->GetXaxis()->GetBinCenter(i+1));
    double dt   = fr[i].p[1]/vdrift*1e3;
    double edt  = fr[i].e[1]/vdrift*1e3;
    double sig  = fr[i].p[2]/vdrift*1e3;
    double esig = fr[i].e[1]/vdrift*1e3;
    
    printf(" %12.5f %12.5f  %12.5f %12.5f  %12.5f\n",
           //           fr[i].p[1],fr[i].e[1], fr[i].p[2],fr[i].e[2], fr[i].chi2dof);
           dt,edt,sig,esig, fr[i].chi2dof);
  }


  return 0;
}
