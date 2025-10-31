///////////////////////////////////////////////////////////////////////////////
// PM: compater TRK timing calibrations :// 1. (TrkStrawPreamp+TrkDelayPanel)
// 1. TrkStrawPreamp: CAL-HV
// 2. (TrkStrawPreamp+TrkDelayPanel)
//
// .L compare_calibrations.C
// auto cc = new CompareCalibrations("v3_992","v3_992",
//                                   "v001/daqana/calibrations/TrkPreampStraw_v3.txt",
//                                   "v001/trkvst/stationfcl/calibrations/run-107992/TrkPreampStraw.txt");
// cc->h_dcal[0]->Draw()
//-----------------------------------------------------------------------------
#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TString.h"
#include "TFolder.h"

#include "daqana/scripts/read_trk_preamp_straw.C"
#include "daqana/obj/digis.hh"

class CompareCalibrations {
public:
  TString fFn1;
  TString fFn2;
  TString fName;
  TString fTitle;

  TFolder* fFolder;

  TDirectory* fDir;

  TH1F* h_dt_cal_hv_p [12];
  TH1F* h_dt_cal_hv_r [12];
  TH1F* h_dt_cal_hv_pr[12];
  TH1F* h_dcal        [12];

  TH1F* h_ctot_calp[12];
  TH1F* h_ctot_hvp [12];
  TH1F* h_ctot_calr[12];
  TH1F* h_ctot_hvr [12];

  std::vector<TrkPreampStrawData_t> fCalib1;
  std::vector<TrkPreampStrawData_t> fCalib2;

//-----------------------------------------------------------------------------
// name should be short
//-----------------------------------------------------------------------------
  CompareCalibrations(const char* Name, const char* Title, const char* Fn1, const char* Fn2, int Panel = -1) {
    double panel_offset_richie[12] = {0.0,   0.0,  0.0, 0.0,   0.0, 0.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
    double panel_offset_pasha [12] = {0.0, -8.00, -9.0, 0.0,  13.0, 3.0,  97.0,  96.0,  95.0, 105.0,  95.0,  95.0};

    TH1::AddDirectory(0);

    fFolder = gROOT->GetRootFolder()->AddFolder(Name,Name);

    fName  = Name;
    fTitle = Title;
    fFn1   = Fn1;
    fFn2   = Fn2;

    printf(">> comparing calibrations from (1):%s and (2):%s\n",Fn1,Fn2);

    read_trk_preamp_straw(Fn1, fCalib1);
    read_trk_preamp_straw(Fn2, fCalib2);

    for (int i=0; i<12; i++) {
//-----------------------------------------------------------------------------
// CAL-HV
//-----------------------------------------------------------------------------
      h_dt_cal_hv_p [i] = new TH1F(Form("h_%s_dt_cal_hvp_%02i" ,fName.Data(),i), Form("%s: h_dt_cal_hv_p_%02i",fTitle.Data(),i),500,-25,25);
      h_dt_cal_hv_r [i] = new TH1F(Form("h_%s_dt_cal_hvr_%02i" ,fName.Data(),i), Form("%s: h_dt_cal_hv_r_%02i",fTitle.Data(),i),500,-25,25);
      h_dt_cal_hv_pr[i] = new TH1F(Form("h_%s_dt_cal_hvpr_%02i",fName.Data(),i), Form("%s: #DeltaT (1)-(2): panel %02i",fName.Data(),i),200,-10,10);

      h_ctot_calp   [i] = new TH1F(Form("h_%s_ctot_calp_%02i",fName.Data(),i), Form("%s: h_ctot_calp_%02i",fTitle.Data(),i),750,-150,150);
      h_ctot_hvp    [i] = new TH1F(Form("h_%s_ctot_hvp_%02i" ,fName.Data(),i), Form("%s: h_ctot_hvp_%02i" ,fTitle.Data(),i),750,-150,150);
      h_ctot_calr   [i] = new TH1F(Form("h_%s_ctot_calr_%02i",fName.Data(),i), Form("%s: h_ctot_calr_%02i",fTitle.Data(),i),750,-150,150);
      h_ctot_hvr    [i] = new TH1F(Form("h_%s_ctot_hvr_%02i" ,fName.Data(),i), Form("%s: h_ctot_hvr_%02i" ,fTitle.Data(),i),750,-150,150);

      h_dcal        [i] = new TH1F(Form("h_%s_dcal_%02i"     ,fName.Data(),i), Form("%s: h_dcal_%02i"     ,fTitle.Data(),i),200,-10,10);
//-----------------------------------------------------------------------------
// add histograms to the folder
//-----------------------------------------------------------------------------
      fFolder->Add(h_dt_cal_hv_p [i]);
      fFolder->Add(h_dt_cal_hv_r [i]);
      fFolder->Add(h_dt_cal_hv_pr[i]);
      fFolder->Add(h_ctot_calp   [i]);
      fFolder->Add(h_ctot_hvp    [i]);
      fFolder->Add(h_ctot_calr   [i]);
      fFolder->Add(h_ctot_hvr    [i]);
      fFolder->Add(h_dcal        [i]);
    }

    // printf("channel_offsets[5]: %10.3f\n",channel_offsets[5]);
    int ipmin =  0;
    int ipmax = 12;
    if (Panel >= 0) {
      ipmin=Panel;
      ipmax=Panel+1;
    }
    for (int ip=ipmin; ip<ipmax; ++ip) {
      int offset = 96*ip;

      for (int i=offset; i<offset+96; ++i) {
        int   ich  = fCalib1[offset+i].ich;

        float cal1 = fCalib1[i].dtcal;
        float hv1  = fCalib1[i].dthv;
        float cal2 = fCalib2[i].dtcal;
        float hv2  = fCalib2[i].dthv;

        //        printf("ich:%i calp:%10.3f hvp:%10.3f calr:%10.3f hvr:%10.3f\n",ich, cal1, hv1,cal2, hv2);
      
        // CAL-HV for different sets
        float dtp = cal1-hv1;
        float dtr = cal2-hv2;
      
        h_dt_cal_hv_p[ip]->Fill(dtp);
        h_dt_cal_hv_r[ip]->Fill(dtr);
        h_dt_cal_hv_pr[ip]->Fill(dtp-dtr);

        float ctot_calp = cal1+panel_offset_pasha [ip];
        float ctot_hvp  = hv1 +panel_offset_pasha [ip];
        float ctot_calr = cal2+panel_offset_richie[ip];
        float ctot_hvr  = hv2 +panel_offset_richie[ip];

        h_ctot_calp[ip]->Fill(ctot_calp);
        h_ctot_hvp [ip]->Fill(ctot_hvp );
        h_ctot_calr[ip]->Fill(ctot_calr);
        h_ctot_hvr [ip]->Fill(ctot_hvr );

        h_dcal[ip]->Fill(ctot_calp-ctot_calr);
      }
    }
  }


//-----------------------------------------------------------------------------
  void SaveHist(const char* Filename) {
    // save histograms booked by all the modules into a file with the given name
    // Mode = 1: save folders
    
    TFile* f = new TFile(Filename,"recreate");
    
    digis::SaveFolder(fFolder,f);
    
    f->Close();
    delete f;
  }

};
