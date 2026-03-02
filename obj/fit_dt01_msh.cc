///////////////////////////////////////////////////////////////////////////////
// fit histograms produced by daqana/mod/MakeStationHist_module.cc
// tested fitting the cosmic data (the second function)
//-----------------------------------------------------------------------------
#define __CLING__ 1
#include <iostream>
#include <format>
#include "daqana/obj/TrkPanelMap_t.hh"
// #include "Offline/DataProducts/inc/StrawId.hh"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "daqana/obj/fit_functions.hh"

//-----------------------------------------------------------------------------
// FirstChannel: first pulsed channel, they go by 8
// 'msh' : histograms produced by make_station_hist job (StationAna module)
// printout for calibration: straw_id, delay_hv, delay_cal, thr_hv, thr_cal, gain
// thr_hv, thr_cal, gain are not used by the reconstruction
// have two runs per calibration, masks 0xba and 0x55
// first channel in the group of 8 is always overlapping (if not masked out)
//-----------------------------------------------------------------------------
int fit_dt01_msh(int Slot, int FirstRun, int Panel1, int Panel2, int FirstChannel,int LastChannel, int PrintLevel) {

  struct fit_result_t {
    double p[3];
    double e[3];
    double chi2dof;
    float  n0;
    float  ntot;
    float  ineff;
  };

  int const nruns(2);
  fit_result_t fr[nruns][36][6][96]; 

  TrkPanelMap_t* tpm = TrkPanelMap_t::Instance(FirstRun);

  const char* hist_dir = "/projects/mu2e/data/projects/vst/hist";

  char fn[200];

  for (int ir=0; ir<nruns; ir++) {
    int run_number = FirstRun+ir;
    sprintf(fn,"%s/hst.mu2e.vst00s000r000n000.make_station_hist.%06d_000001.root",hist_dir,run_number);

    std::cout << std::format("001: run_number:{} fn:{}\n",run_number,fn);
    
    TFile* f = TFile::Open(fn);

    if (PrintLevel == 0) {
    }
    else if (PrintLevel == 1) {
      printf("name idtc ilink ic       mean        emean        sig        esig      chi2dof  n0   ntot   ineff\n");
      printf("-------------------------------------------------------------------------------------------------------\n");
    }
//-----------------------------------------------------------------------------
// plane and panel numbering - offline 
//-----------------------------------------------------------------------------
    for (int plane=2*Slot; plane<2*Slot+2; ++plane) {
      for (int panel=Panel1; panel<Panel2; panel++) {
        TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
//-----------------------------------------------------------------------------
// first channel is the same as ir
//-----------------------------------------------------------------------------
        for (int ich=FirstChannel; ich<LastChannel+1; ich++) {
                                        // fit 'good' histogram'
          char hname[128];
          sprintf(hname,"//StationAna/slot_%02d/MN%03d/str_%02i/ch_%02i_dtchg",Slot,tpmd->mnid,ich,ich);
          TH1F* h = (TH1F*) f->Get(hname);
//-----------------------------------------------------------------------------
// skip empty channels
//-----------------------------------------------------------------------------
          if (h->GetEntries() < 1000) continue;
          int nx = h->GetNbinsX();
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
              tmax = h->GetBinCenter(ix+1);
            }
          }

      // std::cout << "tmax:" << tmax << " nmax:" << nmax << std::endl;

          fit_result_t* frr = &fr[ir][plane][panel][ich];
      
          frr->chi2dof = -1;
        
          for (int ip=0; ip<3; ip++) {
            frr->p[ip] = -999;
            frr->e[ip] = -999;
          }          
          
          h->Draw();

          if (nmax > 0) {
            TFitResultPtr tfr = h->Fit("gaus","sq","",tmax-5,tmax+5);

            if ((! tfr->IsValid()) or tfr->IsEmpty()) {
              printf("# FIT ERROR: channel: %2i\n",ich);
              continue;
            }

            frr->chi2dof = tfr->Chi2()/tfr->Ndf();
            double sf     = sqrt(frr->chi2dof);

            for (int ip=0; ip<3; ip++) {
              frr->p[ip] = tfr->Parameter(ip);
              frr->e[ip] = tfr->Error(ip)*sf;
            }
          }
          
          TH1F* h_nhitsg = (TH1F*) f->Get(Form("//StationAna/slot_%02d/MN%03d/str_%02i/ch_%02i_nhitsg",Slot,tpmd->mnid,ich,ich));
          frr->n0   = h_nhitsg->GetBinContent(1);
          frr->ntot = h_nhitsg->GetEntries();
          frr->ineff = 0;
          if (frr->ntot > 0) frr->ineff = frr->n0/frr->ntot;
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// print fit results
//-----------------------------------------------------------------------------
  for (int plane=2*Slot; plane<2*Slot+2; plane++) {
    for (int panel=0; panel<6; panel++) {
      TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
      for (int ich=FirstChannel; ich<LastChannel+1; ich++) {
        fit_result_t* frr = &fr[0][plane][panel][ich];
        if (PrintLevel == 1) {
          int mnid = tpmd->mnid;
          printf(" %2i    %i   MN%03d  %2i",plane, panel, mnid, ich);
          printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %5.0f %5.0f %10.4f\n",
                 frr->p[1],frr->e[1], frr->p[2], frr->e[2], frr->chi2dof, frr->n0, frr->ntot, frr->ineff);
        }
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// for cosmic runs, process all channels
// this function needs an update
//-----------------------------------------------------------------------------
int fit_dt01_cosmics(int RunNumber, int Panel1, int Panel2, int FirstChannel, int LastChannel, int PrintLevel) {
                                        // initialize panel map

  //  init_trk_panel_map();
 
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
    
    TrkPanelMap_t* tpm       = TrkPanelMap_t::Instance(RunNumber);
    TrkPanelMap_t::Data_t* r = tpm->panel_data_by_mnid(mnid[ipnl]);
    int plane   = r->plane;
    int panel   = r->panel;
    // int station = plane / 2;

    // printf("mnid:%i station:%i plane:%i panel:%i\n",mnid[ipnl],station,plane,panel);
    
    for (int ich=FirstChannel; ich<LastChannel+1; ich++) {
                                        // initialize to undefined (bad)
      fit_result_t* frr = &fr[ich];

      frr->chi2dof = -1;
        
      for (int ip=0; ip<3; ip++) {
        frr->p[ip] = 0;
        frr->e[ip] = 0;
      }          

      // for (int ich=8; ich<9; ich++) {
      
      TH1F* h = (TH1F*) f->Get(Form("//StationAna/pnlset_00/%s/str_%02i/ch_%02i_dtchg",pnl_name[ipnl],ich,ich));
      int  nx = h->GetNbinsX();
      long int nent = h->GetEntries();
      // std::cout << "ich:" << ich << " nent:" << nent << std::endl;
      if (nent > 100) {
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
            tmax = h->GetBinCenter(ix+1);
          }
        }

      // std::cout << "tmax:" << tmax << " nmax:" << nmax << std::endl;

        h->Draw();

        if (nmax > 0) {
          TFitResultPtr tfr = h->Fit("gaus","sq","",tmax-10,tmax+10);

          if ((! tfr->IsValid()) or tfr->IsEmpty()) {
            // printf("# FIT ERROR: channel: %2i\n",ich);
            frr->chi2dof = -1.;
      
            for (int ip=0; ip<3; ip++) {
              frr->p[ip] = -999.;
              frr->e[ip] = -999.;
            }
            // continue;
          }

          frr->chi2dof = tfr->Chi2()/tfr->Ndf();
          double sf     = sqrt(frr->chi2dof);
      
          for (int ip=0; ip<3; ip++) {
            frr->p[ip] = tfr->Parameter(ip);
            frr->e[ip] = tfr->Error(ip)*sf;
          }
        }
      }
//-----------------------------------------------------------------------------
// print fit results, chid is the 'compact' channel index 
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
