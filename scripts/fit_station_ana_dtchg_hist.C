///////////////////////////////////////////////////////////////////////////////
// StationAna produced delta_T(hv-cal) histograms
// those histograms make the second part of input for TrkPreampStraw tables
// need a "reference" channel which is read out and not masked
//
// StationAna dtchg histogtams plot T(CAL)-T(HV), take that into account when updating TrkPreampStraw
// example:

/*
.L v001/daqana/scripts/fit_station_ana_dtchg_hist.C
auto fsa = new FitStationAna(121103)
fsa->fit_dtchg_hist(10)
fsa->fit_dtchg_hist(11)
fsa->fit_dtchg_hist(12)
fsa->write_TrkPreampStraw_table("TrkPreampStraw_v00.txt")
 */

///////////////////////////////////////////////////////////////////////////////
#include "format"
#include "iostream"

#include "TFile.h"
#include "TH1F.h"
#include "daqana/obj/TrkPanelMap_t.hh"


//-----------------------------------------------------------------------------
// step 1: for each straw fit 'dtchg' distributions with a Gaussian,
//         check consistency across different panels
//-----------------------------------------------------------------------------

class FitStationAna {
public:
  
  struct FitResult_t {
    double p[3];
    double e[3];
    double chi2dof;
    float  n0;
    float  ntot;
    float  ineff;
  };

  struct PanelHist_t {
    TH1F*  h_dt;
    TH1F*  h_chi2;
  };

  int            fRunNumber;
  FitResult_t    fFitResults[18][12][96];
  PanelHist_t    fHist[18][12];
  TrkPanelMap_t* fTpm;

  FitStationAna(int RunNumber) {
    fRunNumber = RunNumber;
    fTpm       = TrkPanelMap_t::Instance(RunNumber);
  }
  
  int fit_dtchg_hist  (int Slot, int Panel = -1, int PrintLevel=0);
  int plot_fit_results(int Slot, int Panel = -1);

  int write_TrkPreampStraw_table(const char* Fn = "TrkPreampStraw.txt");
};

//-----------------------------------------------------------------------------
int FitStationAna::fit_dtchg_hist(int Slot, int Panel, int PrintLevel) {
  int rc(0);
  
  std::string hist_dir("/data/mu2e/mu2etrk/hist");
  std::string dsid    ("vst00s000r000n000");
  
  std::string fn = std::format("{}/station_ana/hst.mu2e.{}.make_station_hist.{}_000001.root",
                               hist_dir,dsid,fRunNumber);
  
  TFile* f = TFile::Open(fn.data());

  int ip1(Panel), ip2(Panel+1);
  if (Panel < 0) { ip1 = 0; ip2 = 12; }
  
  for (int ip=ip1; ip<ip2; ip++) {
                                        // "offline coordinates
    int plane = 2*Slot + ip/6;
    int panel = ip % 6;

    TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
    int mnid = tpmd->mnid;

    std::cout << std::format("plane:{} panel:{} mnid:{}\n",plane,panel,mnid);

    std::string panel_path = std::format("//StationAna/slot_{:02d}/MN{:03d}",Slot,mnid);
    
    // for (int ich=0; ich<10; ich++) {
    for (int ich=0; ich<96; ich++) {
      std::string hist_path = std::format("{}/str_{:02d}/ch_{:02d}_dtchg",panel_path,ich,ich);
//-----------------------------------------------------------------------------
// initialize fit results
//-----------------------------------------------------------------------------
      FitResult_t* fr = &fFitResults[Slot][ip][ich];
      
      fr->chi2dof = -1;
        
      for (int ip=0; ip<3; ip++) {
        fr->p[ip] = -999;
        fr->e[ip] = -999;
      }

      // std::cout << std::format("hist_path:{}\n",hist_path);
      
      TH1F* h  = (TH1F*) f->Get(hist_path.data());
      int   nx = h->GetNbinsX();

      // 1. draw
      // h->Draw();

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
//-----------------------------------------------------------------------------
// perform fit
//-----------------------------------------------------------------------------
        if (nmax > 0) {
          TFitResultPtr tfr = h->Fit("gaus","sq","",tmax-10,tmax+10);

          if ((! tfr->IsValid()) or tfr->IsEmpty()) {
            // printf("# FIT ERROR: channel: %2i\n",ich);
            fr->chi2dof = -1.;
            
            for (int ip=0; ip<3; ip++) {
              fr->p[ip] = -999.;
              fr->e[ip] = -999.;
            }
            // continue;
          }
//-----------------------------------------------------------------------------
// renormalize errors returned by the fit to make effective chi2/dof = 1
//-----------------------------------------------------------------------------
          fr->chi2dof = tfr->Chi2()/tfr->Ndf();
          double sf    = sqrt(fr->chi2dof);
          
          for (int ip=0; ip<3; ip++) {
            fr->p[ip] = tfr->Parameter(ip);
            fr->e[ip] = tfr->Error(ip)*sf;
          }
        }
        
        // TH1F* h_nhitsg = (TH1F*) f->Get(Form("//StationAna/pnlset_00/%s/str_%02i/ch_%02i_nhitsg",pnl_name[ipnl],ich,ich));
        // frr->n0   = h_nhitsg->GetBinContent(1);
        // frr->ntot = h_nhitsg->GetEntries();
        // frr->ineff = 0;
        // if (frr->ntot > 0) frr->ineff = frr->n0/frr->ntot;
      }
    }
  }

//-----------------------------------------------------------------------------
// print fit resutls
//-----------------------------------------------------------------------------
  for (int ip=ip1; ip<ip2; ip++) {
    int plane = 2*Slot+ip/6;
    int panel = ip%6;

    TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
    int mnid = tpmd->mnid;
    
    for (int ich=0; ich<96; ich++) {
      FitResult_t* fr = &fFitResults[Slot][ip][ich];
      if (PrintLevel == 1) {
        std::cout << std::format("MN{:03d} {:2d} {:2d} {:2d}",mnid,plane,panel,ich);
        printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %5.0f %5.0f %10.4f\n",
               fr->p[1],fr->e[1], fr->p[2], fr->e[2], fr->chi2dof, fr->n0, fr->ntot, fr->ineff);
      }
    }
  }

  return rc;
}

//-----------------------------------------------------------------------------
int FitStationAna::plot_fit_results(int Slot, int Panel) {
  int rc(0);

  int ip1(Panel), ip2(Panel+1);
  
  std::string name = std::format("c_slot_{:02}_{:02}",Slot,Panel);
  
  if (Panel < 0) {
    ip1 = 0; ip2 = 12;
    name = std::format("c_slot_{:02}_all",Slot);
  }
  
  TCanvas* c = new TCanvas(name.data(),"c",1000,1000);
  c->Divide(1,2);

  bool first_hist(true);

  for (int ip=ip1; ip<ip2; ip++) {
    
    PanelHist_t* ph = &fHist[Slot][ip];
  
    std::string hist_name, title;

    hist_name  = std::format("h_dt_slot_{:02d}_panel_{:02d}",Slot,ip);
    title      = std::format("run {:06d} {}",fRunNumber,hist_name);    // for now
    ph->h_dt   = new TH1F(hist_name.data(),title.data(),100,0,100);
    ph->h_dt->SetMarkerStyle(20);
    ph->h_dt->SetMarkerSize(0.6);

    hist_name  = std::format("h_chi2_slot_{:02d}_panel_{:02d}",Slot,ip);
    title      = std::format("run {:06d} {}",fRunNumber,hist_name);    // for now
    ph->h_chi2 = new TH1F(hist_name.data(),title.data(),100,0,100);
    ph->h_chi2->SetMarkerStyle(20);
    ph->h_chi2->SetMarkerSize(0.6);

    for (int ich=0; ich<96; ich++) {
      FitResult_t* fr = &fFitResults[Slot][ip][ich];
      ph->h_chi2->SetBinContent(ich+1,fr->chi2dof);
      ph->h_chi2->SetBinError  (ich+1,0);

      if (fr->chi2dof > 0) {
        ph->h_dt->SetBinContent(ich+1,fr->p[1]);
        ph->h_dt->SetBinError  (ich+1,fr->e[1]);
      }
      
    }

    if (first_hist) {
      c->cd(1);
      ph->h_dt->Draw("pe");
      c->cd(2);
      ph->h_chi2->Draw("pe");
      first_hist = false;
    }
    else {
      c->cd(1);
      ph->h_dt->Draw("pe,same");
      c->cd(2);
      ph->h_chi2->Draw("pe,same");
    }
  }

  c->Modified();
  c->Update();
  return rc;
}

//-----------------------------------------------------------------------------
// create "trivial" TrkPreampStraw table using only TCAL-THV offsets,
// no channel to channel offsets yet
//-----------------------------------------------------------------------------
int FitStationAna::write_TrkPreampStraw_table(const char* Fn) {
  int rc(0);
  
  std::ofstream file(Fn);

  file << std::format("TABLE TrkPreampStraw   120000-200000\n");
  
  // offline numbering
  for (int plane=0; plane<36; plane++) {
    int ist = plane / 2;
    
    for (int panel=0; panel<6; panel++) {
      
      file << std::format("#--------------------------\n");
      file << std::format("# plane:{} panel:{}\n",plane, panel);
      file << std::format("#--------------------------\n");
      
      int ipn = panel + 6*(plane%2);
      for (int straw=0; straw<96; straw++) {
        int loc = straw+96*(panel+6*plane);
        float dt_cal = 0;
        // what is the right sign ? looks like minus - the correction is ADDED to the raw time
        float dt_hv  = -fFitResults[ist][ipn][straw].p[1];
        if (fFitResults[ist][ipn][straw].chi2dof < 0) dt_hv = 0.; // not to confuse things too much...
        
        file << std::format("{:5d}, {:10.3f}, {:10.3f}, 12.0, 12.0, 1528000.0\n",
                            loc,dt_cal,dt_hv);
      }
    }
  }
  
  file.close();
  
  return rc;
}
