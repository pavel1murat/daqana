///////////////////////////////////////////////////////////////////////////////
// fit dt01 histograms produced by daqana/mod/MakeStationHist_module.cc
// 
// tested fitting the cosmic data (the second function)
//-----------------------------------------------------------------------------
#define __CLING__ 1
#include <iostream>
#include "daqana/mod/mod/TrkPanelMap_t.hh"
// #include "Offline/DataProducts/inc/StrawId.hh"


#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

struct fit_result_t {
  int    ich;                           // defined on input
  double p[3];
  double e[3];
  double chi2dof;
  float  n0;
  float  ntot;
  float  ineff;
};

//-----------------------------------------------------------------------------
// initialize the panel map for a given run
//-----------------------------------------------------------------------------
const TrkPanelMap_t* _panel_map[500];   // indexing - online (sparse): [mnid]

int init_trk_panel_map() {
  
  for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
    int mnid = tpm->mnid;
    _panel_map[mnid] = tpm;
    
    // int idtc   = tpm->dtc % 2;
    // int ilink  = tpm->link;
    // int ipanel = 6*idtc+ilink;
    
    // _edata.panel[ipanel].mnid = tpm->mnid;
  }
  return 0;
}

//-----------------------------------------------------------------------------
int fit_histogram(TH1F* Hist, fit_result_t* Frr) {
  int rc(0);
    
  Frr->chi2dof = -1;
        
  for (int ip=0; ip<3; ip++) {
    Frr->p[ip] = -999;
    Frr->e[ip] = -999;
  }          

  Hist->Draw();

  int nx   = Hist->GetNbinsX();
  int nent = Hist->GetEntries();
  if (nent < 100) {
    return -1;
  }
//-----------------------------------------------------------------------------
// find bin with max content and fit +/- 5 ns from it
//-----------------------------------------------------------------------------
  double nmax = -1.;
  float  tmax = -1e6;
  for (int ix=0; ix<nx; ix++) {
    double y = h->GetBinContent(ix+1);
    // std::cout << "ix:" << ix << " y:" << y << std::endl;
    if (y > nmax) {
      nmax = y;
      tmax = Hist->GetBinCenter(ix+1);
    }
  }

  
  // std::cout << "tmax:" << tmax << " nmax:" << nmax << std::endl;

  if (nmax > 0) {
    TFitResultPtr tfr = Hist->Fit("gaus","sq","",tmax-10,tmax+10);

    if ((! tfr->IsValid()) or tfr->IsEmpty()) {
      printf("# FIT ERROR: channel: %2i\n",Frr->ich);
      rc = -2;
      return rc;
    }
    
    Frr->chi2dof = tfr->Chi2()/tfr->Ndf();
    double sf     = sqrt(Frr->chi2dof);
          
    for (int ip=0; ip<3; ip++) {
      Frr->p[ip] = tfr->Parameter(ip);
      Frr->e[ip] = tfr->Error(ip)*sf;
    }
  }
  return rc;
}

//-----------------------------------------------------------------------------
// FirstChannel: first pulsed channel, they go by 8
// 'msh' : histograms produced by make_station_hist job (StationAna module)
// printout for calibration: straw_id, delay_hv, delay_cal, thr_hv, thr_cal, gain
// thr_hv, thr_cal, gain are not used by the reconstruction
// pin: pulse injection
//-----------------------------------------------------------------------------
int fit_dt01_msh(int FirstRun, int Panel1=0, int Panel2=12, int PrintLevel=1) {

  fit_result_t fr[2][6][96]; 

  const char* pnl_name[12] = { "MN261", "MN248", "MN224", "MN262", "MN273", "MN276",
                               "MN253", "MN101", "MN219", "MN213", "MN235", "MN247"
  };

  const char* hist_dir = "/data/tracker/vst/hist";

  char fn[200];

  for (int ir=0; ir<8; ir++) {
    sprintf(fn,"%s/hst.mu2e.vst00s000r000n000.make_station_hist.%06d_000001.root",hist_dir,FirstRun+ir);

    std::cout << "001: ir:" << ir << " fn=" << fn << std::endl;
    
    TFile* f = TFile::Open(fn);

    if (PrintLevel == 0) {
    }
    else if (PrintLevel == 1) {
      printf("name idtc ilink ic       mean        emean        sig        esig      chi2dof  n0   ntot   ineff\n");
      printf("-------------------------------------------------------------------------------------------------------\n");
    }

    //   for (int ipnl=0; ipnl<12; ipnl++) {
    for (int ipnl=Panel1; ipnl<Panel2; ipnl++) {
      int pcie_addr = ipnl / 6;
      int link      = ipnl % 6;
//-----------------------------------------------------------------------------
// first channel is the same as ir
//-----------------------------------------------------------------------------
      for (int ich=ir; ich<96; ich+=8) {
        fit_result_t* frr = &fr[pcie_addr][link][ich];
        frr->ich = ich;
 
        TH1F* h = (TH1F*) f->Get(Form("//StationAna/pnlset_00/%s/str_%02i/ch_%02i_dtchg",pnl_name[ipnl],ich,ich));

        int rc = fit_histogram(h,frr);
        if (rc < 0) continue;
        
        TH1F* h_nhitsg = (TH1F*) f->Get(Form("//StationAna/pnlset_00/%s/str_%02i/ch_%02i_nhitsg",pnl_name[ipnl],ich,ich));
        frr->n0   = h_nhitsg->GetBinContent(1);
        frr->ntot = h_nhitsg->GetEntries();
        frr->ineff = 0;
        if (frr->ntot > 0) frr->ineff = frr->n0/frr->ntot;
      }
    }
  }
//-----------------------------------------------------------------------------
// print fit resutls
//-----------------------------------------------------------------------------
  for (int pcie=0; pcie<2; pcie++) {
    for (int link=0; link<6; link++) {
      for (int ich=0; ich<96; ich++) {
        fit_result_t* frr = &fr[pcie][link][ich];
        if (PrintLevel == 1) {
          printf("%5s   %i    %i   %2i",pnl_name[6*pcie+link], pcie, link, ich);
          printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %5.0f %5.0f %10.4f\n",
                 frr->p[1],frr->e[1], frr->p[2], frr->e[2], frr->chi2dof, frr->n0, frr->ntot, frr->ineff);
        }
      }
    }
  }

  return 0;
}


//-----------------------------------------------------------------------------
// pin: pulse injection
//-----------------------------------------------------------------------------
int fit_dt01_pin(int Station, int Panel, int PrintLevel=1) {

  fit_result_t fr[216][96]; 

  const char* pnl_name[12] = { "MN261", "MN248", "MN224", "MN262", "MN273", "MN276",
                               "MN253", "MN101", "MN219", "MN213", "MN235", "MN247"
  };

  const char* hist_dir = "./pasha";

  char fn[200];

  int run1(122655), run2(122662);

  RunInfoDb::Data_t* rid;  

  for (int ir=0; ir<8; ir++) {
    
    sprintf(fn,"%s/pin_%06d_%06d.hist",hist_dir,run1,run2);

    std::cout << "001: ir:" << ir << " fn=" << fn << std::endl;
    
    TFile* f = TFile::Open(fn);

    if (PrintLevel == 0) {
    }
    else if (PrintLevel == 1) {
      printf("name idtc ilink ic       mean        emean        sig        esig      chi2dof  n0   ntot   ineff\n");
      printf("-------------------------------------------------------------------------------------------------------\n");
    }

    int pnl1(Station*12+Panel), pnl2(Station*12+Panel+1);
    if (Panel == -1) pnl2 = Station*12+11;
    
    //   for (int ipnl=0; ipnl<12; ipnl++) {
    for (int ipnl=pnl1; ipnl<pnl2; ipnl++) {
//-----------------------------------------------------------------------------
// first channel is the same as iri
//-----------------------------------------------------------------------------
      int nch = rid->n_pulsed_channels;
      for (int i=0; i<nch; i++) {
        int ich = rid->channel[i];
                                        // for each run process its channels
        fit_result_t* frr = &fr[ipnl][ich];
        frr->ich = ich;
 
        TH1F* h = (TH1F*) f->Get(Form("//pin/%06d/pnl_%03d/ch_%02d/dt01",run_number,panel,ich));

        int rc = fit_histogram(h,frr);
      }
    }
  }
//-----------------------------------------------------------------------------
// print fit resutls
//-----------------------------------------------------------------------------
  for (int pcie=0; pcie<2; pcie++) {
    for (int link=0; link<6; link++) {
      for (int ich=0; ich<96; ich++) {
        fit_result_t* frr = &fr[pcie][link][ich];
        if (PrintLevel == 1) {
          printf("%5s   %i    %i   %2i",pnl_name[6*pcie+link], pcie, link, ich);
          printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %5.0f %5.0f %10.4f\n",
                 frr->p[1],frr->e[1], frr->p[2], frr->e[2], frr->chi2dof, frr->n0, frr->ntot, frr->ineff);
        }
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// for cosmic run, process all channels
//-----------------------------------------------------------------------------
int fit_dt01_cosmics(int RunNumber, int Panel1=0, int Panel2=12, int FirstChannel=0, int LastChannel=95, int PrintLevel=0) {
                                        // initialize panel map

  init_trk_panel_map();
 
  char fn[200];

  const char* hist_dir = "/data/tracker/vst/hist";
  sprintf(fn,"%s/hst.mu2e.vst00s000r000n000.make_station_hist.%06d_000001.root",hist_dir,RunNumber);

  TFile* f = TFile::Open(fn);

  struct fit_result_t {
    double p[3];
    double e[3];
    double chi2dof;
  };

  fit_result_t fr[96]; 

  const char* pnl_name[12] = { "MN261", "MN248", "MN224", "MN262", "MN273", "MN276",
                               "MN253", "MN101", "MN219", "MN213", "MN235", "MN247"
  };

  const int mnid[12] = { 261, 248, 224, 262, 273, 276,
                         253, 101, 219, 213, 235, 247
  };

  if (PrintLevel == 0) {
    printf("ipnl  name idtc ilink ic       mean        emean        sig        esig      chi2dof\n");
    printf("------------------------------------------------------------------------------------\n");
  }
  else if (PrintLevel == 1) {
    printf("  sid  delay_hv    delay_cal    thr_hv     thr_cal    gain\n");
    printf("----------------------------------------------------------\n");
  }
  else if (PrintLevel == 2) {
    printf("  sid  dt      err_dr      sig_dt     err_sig_dt   chi2 \n");
    printf("--------------------------------------------------------\n");
  }

  for (int ipnl=Panel1; ipnl<Panel2; ipnl++) {
    
    const TrkPanelMap_t* tpm = _panel_map[mnid[ipnl]];
    int station = tpm->station;
    int plane   = tpm->plane;
    int panel   = tpm->panel;

    // printf("mnid:%i station:%i plane:%i panel:%i\n",mnid[ipnl],station,plane,panel);
    
    for (int ich=FirstChannel; ich<LastChannel+1; ich++) {
                                        // initialize to undefined (bad)
      fit_result_t* frr = &fr[ich];
      frr->ich = ich'

      TH1F* h = (TH1F*) f->Get(Form("//StationAna/pnlset_00/%s/str_%02i/ch_%02i_dtchg",pnl_name[ipnl],ich,ich));
//-----------------------------------------------------------------------------
// fit histogram
//-----------------------------------------------------------------------------
      int rc = fit_histogram(h,frr);
//-----------------------------------------------------------------------------
// print fit results, chid is the 'compact' channel index used in calib DB
//-----------------------------------------------------------------------------
      // mu2e::StrawId sid(plane,panel,ich);

      int chid = ich+96*(panel+plane*6);

      if (PrintLevel == 1) {
        printf("%5i,%8.3f,%8.3f,%8.3f,%8.3f,%10.3f\n",chid /*sid.asUint16()*/, frr->p[1], 0., 12., 12., 70000.);
      }
      else if (PrintLevel == 2) {
        printf("%5i,%8.3f,%8.3f,%8.3f,%8.3f,%10.3f\n",chid, frr->p[1], frr->e[1], frr->p[2], frr->e[2], frr->chi2dof);
      }
    }
  }

  return 0;
}
