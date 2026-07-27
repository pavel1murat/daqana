///////////////////////////////////////////////////////////////////////////////
// fit dt01 histograms produced by daqana/mod/MakeStationHist_module.cc
// 
// tested fitting the cosmic data (the second function)
// TrkPreampStraw sign convention: stored are delta(T01) = T(cal)-T(hv) ,
// so better plot histogrrams for those differences, to avoid inverting the sign
// when producing the constants.
// for channels with no / failed calibrations, store zeros
//-----------------------------------------------------------------------------

#include <iostream>
#include <fstream>

#include "daqana/obj/TrkPanelMap_t.hh"
#include "daqana/obj/RunDb.hh"
// #include "Offline/DataProducts/inc/StrawId.hh"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "daqana/obj/RunDb.hh"

#include "ana/fit_dt01.hh"

//-----------------------------------------------------------------------------
int fit_dt01::fit_histogram(TH1F* Hist, fit_result_t* Frr) {
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
    double y = Hist->GetBinContent(ix+1);
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
int fit_dt01::fit_msh(int FirstRun, int Panel1=0, int Panel2=12, int PrintLevel=1) {

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
// pin: pulse injection histograms
// Panel - offline numbering (0-215)
// RunNumber : array of 8 run numbers
//-----------------------------------------------------------------------------
int fit_dt01::fit_pin(const char* Fn, int* RunNumber, int SlotLo, int SlotHi, int Panel, int PrintLevel=1) {

  // const char* pnl_name[12] = { "MN261", "MN248", "MN224", "MN262", "MN273", "MN276",
  //                              "MN253", "MN101", "MN219", "MN213", "MN235", "MN247"
  // };

  const char* hist_dir = "/data/mu2e/mu2etrk/hist/plot_n001_pin";

  char fn[200];
  
  sprintf(fn,"%s/hst.mu2e.trk.plot_n001_pin.%06d_%06d.hist",hist_dir,RunNumber[0],RunNumber[7]);
  TFile* f = TFile::Open(Fn);

  RunDb::Data_t rdt;  

  TrkPanelMap_t* tpm = TrkPanelMap_t::Instance(RunNumber[0]);

  for (int ir=0; ir<8; ir++) {
    int rn = RunNumber[ir];
    RunDb::Instance()->GetRunInfo(rn,&rdt);

    int nch = 12; // rdt.n_pulsed_channels;
    
    //    std::cout << "001: ir:" << ir << " fn=" << fn << std::endl;
    
    for (int slot=SlotLo; slot<SlotHi+1; slot++) {

      for (int pnl=0; pnl<12; pnl++) {
//-----------------------------------------------------------------------------
// first channel is the same as iri
//-----------------------------------------------------------------------------
        int panel_index = 12*slot+pnl;
        for (int i=0; i<nch; i++) {
          int ich = rdt.pulsed_channel[i];
                                        // for each run process its channels
          fit_result_t* frr = &fFr[panel_index][ich];
          frr->ich = ich;
 
          TH1F* h = (TH1F*) f->Get(Form("//plot_n001_pin/%06d/pnl_%03d/ch_%02d/dt01",rn,panel_index,ich));

          /*int rc = */ fit_histogram(h,frr);
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// print fit results
//-----------------------------------------------------------------------------
  if (PrintLevel == 1) {
    printf("pnl216 slot plane panel MnID ich       mean        emean        sig        esig      chi2dof  n0   ntot   ineff\n");
    printf("-------------------------------------------------------------------------------------------------------\n");
  }

  for (int slot=SlotLo; slot<SlotHi+1; slot++) {
    for (int pnl12=0; pnl12<12; pnl12++) {
      int panel_index = 12*slot+pnl12;
      int plane = slot*2 + pnl12/6;
      int panel = pnl12 % 6;
      
      TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
      int mnid  = tpmd->mnid;
      
      for (int ich=0; ich<96; ich++) {
        fit_result_t* frr = &fFr[panel_index][ich];
        if (PrintLevel == 1) {
          printf("%4i   %4i %5i %5i MN%03d %3i",panel_index,slot,plane,panel,mnid,ich);
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
int fit_dt01::fit_cosmics(int RunNumber, int Slot, int Panel, int Channel, const char* Fn, int PrintLevel) {
  // int RunNumber, int Panel1=0, int Panel2=12, int FirstChannel=0, int LastChannel=95, int PrintLevel=0) {
                                        // initialize panel map

  //   init_trk_panel_map();
 
  // char fn[200];

  const char* hist_dir = "/data/mu2e/mu2etrk/hist/plot_occup";
  
  std::string fn;
  if (Fn == nullptr) fn = std::format("{}/hst.mu2e.vst00s000r000n007.plot_n007_occup.{:6d}_000001.root",hist_dir,RunNumber);
  else               fn = std::format("{}",Fn);

  TFile* f = TFile::Open(fn.data());
  RunDb::Data_t rdt;
  
  RunDb::Instance()->GetRunInfo(RunNumber,&rdt);
  
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

  TrkPanelMap_t* tpm = TrkPanelMap_t::Instance(RunNumber);

  int slot_lo(Slot), slot_hi(Slot);
  if (Slot == -1) {
    slot_lo = 0; slot_hi=17;
  }

  int panel_lo(Panel), panel_hi(Panel);
  if (Panel < 0) {
    panel_lo =  0;
    panel_hi = 11;
  }

  int ch_lo(Channel), ch_hi(Channel);
  if (Channel < 0) {
    ch_lo =  0;
    ch_hi = 95;
  }

  for (int slot=slot_lo; slot<slot_hi+1; slot++) {
    for (int pnl12=panel_lo; pnl12<panel_hi+1; pnl12++) {

      int panel = pnl12 % 6;
      int plane = 2*slot+(pnl12 / 6);
      
      // TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
                                                          
    // printf("mnid:%i station:%i plane:%i panel:%i\n",mnid[ipnl],station,plane,panel);
    
      for (int ich=ch_lo; ich<ch_hi+1; ich++) {
        // n007 : histograms in hist_00 have PH>50 cut
        std::string hist_name = std::format("//n007_occup/{:06d}/hist_00/slot_{:02d}/pnl_{:03d}/ch_{:02d}/dt01",
                                            RunNumber,slot,pnl12,ich);
        std::string clone_name = std::format("h_clone_{}",hist_name);
        TH1F* h_clone = (TH1F*) f->Get(hist_name.data())->Clone(clone_name.data());
        h_clone->Rebin(5);
//-----------------------------------------------------------------------------
// fit histogram
//-----------------------------------------------------------------------------
        int panel_index   = 12*slot+pnl12;
        fit_result_t* frr = &fFr[panel_index][ich];
        frr->ich          = ich;
        int rc            = fit_histogram(h_clone,frr);
        // delete h_clone;
      }
    }
  }
//-----------------------------------------------------------------------------
// print fit results
//-----------------------------------------------------------------------------
  if (PrintLevel == 1) {
    printf("pnl216 slot plane panel MnID ich       mean        emean        sig        esig      chi2dof  n0   ntot   ineff\n");
    printf("-------------------------------------------------------------------------------------------------------\n");
  }

  for (int slot=slot_lo; slot<slot_hi+1; slot++) {
    for (int pnl12=panel_lo; pnl12<panel_hi+1; pnl12++) {
      int panel_index = 12*slot+pnl12;
      int plane       = slot*2 + pnl12/6;
      int panel       = pnl12 % 6;
      
      TrkPanelMap_t::Data_t* tpmd = tpm->panel_data_by_offline(plane,panel);
      int mnid  = tpmd->mnid;
      
      for (int ich=ch_lo; ich<ch_hi+1; ich++) {
        fit_result_t* frr = &fFr[panel_index][ich];
        if (PrintLevel == 1) {
          printf("%4i   %4i %5i %5i MN%03d %3i",panel_index,slot,plane,panel,mnid,ich);
          printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %5.0f %5.0f %10.4f\n",
                 frr->p[1],frr->e[1], frr->p[2], frr->e[2], frr->chi2dof, frr->n0, frr->ntot, frr->ineff);
        }
      }
    }
  }

  return 0;
}


//-----------------------------------------------------------------------------
// to be with the pulse injection
//-----------------------------------------------------------------------------
int fit_dt01::write_TrkPreampStraw(int RunNumber, const char* Fn) {
  int rc(0);
  
  // as printed from Mu2e offline
  float straw_length[96] = {
    1175.755, 1171.692, 1167.580, 1163.421, 1159.213, 1154.956, 1150.649, 1146.292, 1141.884, 1137.425,
    1132.914, 1128.350, 1123.733, 1119.061, 1114.336, 1109.554, 1104.717, 1099.823, 1094.872, 1089.862,
    1084.793, 1079.664, 1074.474, 1069.222, 1063.908, 1058.530, 1053.088, 1047.580, 1042.005, 1036.363,
    1030.652, 1024.871, 1019.019, 1013.095, 1007.097, 1001.024,  994.874,  988.647,  982.341,  975.954,
     969.484,  962.931,  956.292,  949.565,  942.749,  935.841,  928.840,  921.743,  914.549,  907.255,
     899.858,  892.356,  884.746,  877.025,  869.191,  861.241,  853.171,  844.977,  836.656,  828.205,
     819.619,  810.894,  802.025,  793.008,  783.837,  774.507,  765.012,  755.347,  745.503,  735.475,
     725.254,  714.833,  704.202,  693.351,  682.271,  670.950,  659.375,  647.533,  635.408,  622.985,
     610.246,  597.169,  583.733,  569.911,  555.675,  540.993,  525.827,  510.133,  493.862,  476.955,
     459.341,  440.935,  421.634,  401.308,  379.794,  356.877 
  };

  float v(250.); // mm/ns
  
  std::ofstream file(Fn);   // creates/truncates file

  RunDb::Data_t     fRunInfo;
  RunDb::Instance()->GetRunInfo(RunNumber,&fRunInfo);

  for (int pnl=0; pnl<216; pnl++) {
    for (int ich=0; ich<96; ich++) {
      fit_result_t* frr = &fFr[pnl][ich];
      double dt01(0);
      if (frr->chi2dof > 0) dt01 = frr->p[1];
      if (fRunInfo.run_type == RunDb::e_PULSE_INJECTION) {
        // calibration based on pulse injection: correct deltaT for the signal propagation
        dt01 += straw_length[ich]/v;
      }
      int icch = 96*pnl+ich;
      file << std::format("{:5},  {:11.4f},   0,   12.0,   0,  12.0, 1528000.0\n",icch,dt01);
    }
  }
  return rc;
}
