///////////////////////////////////////////////////////////////////////////////
// plot differences between the average hit times reconstructed in different planes
// in a timecluster.
// this gives a fairly good approximation for the time offsets between teh different
// planes which need to be calibrated out
// take plane 0 as a reference plane and define all offsets with rhespect to it
// then, within each plane, can look at the differences between the panels
// and the average
// calibration code ADDS the TrkDelayPanel timing offset to the hit time,
// so the sign of the dt05[ip][ip-1] has to be inversed in the calibration table
/////////////////////////////////////////////////////////////////////////////////
/*
root [0] gSystem->Load("v001/.spack-env/view/lib/libdaqana_ana.so");
root [1] auto x = new plot_n002_tc(124155,"results/2026-08-22-15-43.make_n0041.mu2e-dl-01.fnal.gov.2403029/nts.mu2e.trk.vst00s000r000n004.124155_000001.root","004")
tree: 0xd8ea4f0
(plot_n002_tc *) 0x3041c90
root [2] x->Loop()
nentries:47331
*/
#include "ana/plot_n002_tc.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_n002_tc"

//-----------------------------------------------------------------------------
plot_n002_tc::plot_n002_tc(int RunNumber, const char* Fn, const char* Label) :
  TNamed(Form("run_%06d_n002_tc",RunNumber),Form("run_%06d_%s_tc",RunNumber,Label)),
  fChain(0) {

  std::string dir = std::format("/data/mu2e/mu2etrk/datasets/vst00s000r000{}",Label); 
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000{}.{:06d}_000001.root",dir,Label,RunNumber);
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
    if (!f || !f->IsOpen()) {
      f = new TFile(fn.data());
    }
  }

  fRunNumber = RunNumber;
  fEvent     = nullptr;
//-----------------------------------------------------------------------------
// pulsed channels
//-----------------------------------------------------------------------------
  // RunData_t rd;
  // rd.run_number        = 122629;
  //  rd.ref_channel       = 13;

  //  fRefChannel          = 21;
  const char* name = GetName();

  fTopFolder = (TFolder*) gROOT->GetRootFolder()->FindObject(name);
  
  if (fTopFolder == nullptr) {
    fTopFolder = gROOT->GetRootFolder()->AddFolder(name,name);
  }

  std::string rns = std::to_string(fRunNumber);
  
  fRunFolder = fTopFolder->AddFolder(rns.data(),rns.data());
//-----------------------------------------------------------------------------
// allow histograms in different folders to have the same name
//-----------------------------------------------------------------------------
  TH1::AddDirectory(0);

  fTpm       = TrkPanelMap_t::Instance(RunNumber);
  fBook      = new Booking(fRunFolder);

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  BookHistograms(fRunFolder);
}

//-----------------------------------------------------------------------------
plot_n002_tc::~plot_n002_tc() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n002_tc::BookTimeClusterHistograms(TimeClusterHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d}",fRunNumber,Index->slot);
  std::string name, title;

  name  = "t0";
  title = std::format("{} : t0",prefix);
  fBook->HBook1F(Hist->h_t0,name.data(),title.data(),1000,1,200,Folder);   // in us...

  return 0;
}



//-----------------------------------------------------------------------------
int plot_n002_tc::BookHistograms(TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  Index_t index;

  fHist = new Hist_t;

  int book_tc_histset[10];
  int n_tc_histsets(10);

  for (int i=0; i<n_tc_histsets; i++) { book_tc_histset[i] = 0; }

  book_tc_histset[0] = 1;

  for (int i=0; i<n_tc_histsets; i++) {
    if (book_tc_histset[i] == 0) continue;
    std::string folder_name = std::format("tc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    fHist->tc[i] = new TimeClusterHist_t();
    BookTimeClusterHistograms(fHist->tc[i],&index,fol);
  }

  for (int i=1; i<36; i++) {
    for (int j=i-1; j>=0; j--) {
      name  = std::format("dt05_{:02}_{:02}",i,j);
      title = std::format("{} : T{:02} - T{:02}",prefix,i,j);
      fBook->HBook1F(fHist->h_dt05[i][j],name.data(),title.data(),1000,-2500,2500,Folder);
    }
  }

  for (int i=0; i<216; i++) {
    name  = std::format("panel_dt_{:03}",i);
    title = std::format("{} : T{:03} - T(tc_best)",prefix,i);
    fBook->HBook1F(fHist->h_panel_dt[i],name.data(),title.data(),200,-100,100,Folder);
  }

  name  = std::format("pdt_all");
  title = std::format("panel dt, all",prefix);
  fBook->HBook1F(fHist->h_pdt_216,name.data(),title.data(),200,-10,10,Folder);

  name  = std::format("dt05_36");
  title = std::format("plane dT = T(i)-T(i-1) vs plane",prefix);
  fBook->HBook2F(fHist->h_dt05_36,name.data(),title.data(),200,-200,200,36,0,36,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n002_tc::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n002_tc::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
int plot_n002_tc::FillTimeClusterHistograms(TimeClusterHist_t* Hist, DaqTimeCluster* Tc) {
  
  Hist->h_t0->Fill(Tc->t0);
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_tc::FillHistograms() {
  // filling histograms: plot time differences between

  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;

  //  Index_t index;
  
  for (int i=0; i<fEvent->ntc; i++) {
    DaqTimeCluster*  tc = (DaqTimeCluster* ) fEvent->tc->UncheckedAt(i);

    FillTimeClusterHistograms(fHist->tc[0],tc);
  }; 

  for (int i=0; i<216; i++) {
    if (fPanelNh[i] > 0) {
      fHist->h_panel_dt[i]->Fill(fPanelDt[i]);
    }
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
Long64_t plot_n002_tc::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

//-----------------------------------------------------------------------------
void plot_n002_tc::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;

  fNEvents  = nev;
  fMaxEvent = -1;

  // int current_event = -1;
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
         
    if (fEvent->evn > fMaxEvent) {
      fMaxEvent = fEvent->evn;
    }
//-----------------------------------------------------------------------------
//  calculate ave
//-----------------------------------------------------------------------------
    // and reinitializa
    for (int i=0; i<36; i++) {
      t05[i] = 0.;
      n05[i] = 0;
    }
    
    // if conly "close" hits are stored, the number of hits in the list is less than the total
    // number of hits reconstructed in the event
    int nsh = fEvent->sh->GetEntriesFast();
    for (int i=0; i<nsh; i++) {
      DaqStrawHit*  sh = fEvent->Sh(i);
      if (sh->edep > 0.0005) {
        int ip = sh->plane();
        t05[ip] += sh->time;
        n05[ip] += 1;
      }
    }
    // average times
    for (int i=0; i<36; i++) {
      t05[i] = t05[i]/(n05[i]+1.e-12);
    }
    
    // and calculate their residuals
    for (int i=0; i<35; i++) {
      for (int j=i+1; j<36; j++) {
        if ((n05[i] > 1) and (n05[j] > 1)) {
          dt05[j][i] = t05[j]-t05[i];
          fHist->h_dt05[j][i]->Fill(dt05[j][i]);
        }
      }
    }
    // station 7 (planes 14 and 15 missing), need 16-13
    for (int i=1; i<36; i++) {
      if ((i < 14) or (i > 16)) {
        fHist->h_dt05_36->Fill(dt05[i][i-1],i);
      }
      else if (i == 16) {
        fHist->h_dt05_36->Fill(dt05[16][i-3],16);
      }
    }

    // identify time cluster with the largest number of hits
    int             nsh_best(0);
    DaqTimeCluster* tc_best(nullptr);
    
    for (int itc=0; itc<fEvent->ntc; itc++) {
      auto tc = fEvent->Tc(itc);
      if ((tc->nsh > nsh_best) and (tc->nplanes > 3)) {
        tc_best  = tc;
        nsh_best = tc->nsh;
      }
    }

    if (tc_best) {
      // for each panel plot the time difference between the combohit time and the TC t0
      int nch = fEvent->ch->GetEntriesFast();
      
      for (int i=0; i<216; i++) {
        fPanelNh[i] = 0;
      }
      
      for (int i=0; i<nch; i++) {
        DaqComboHit*  ch = fEvent->Ch(i);
        if (ch->edep > 0.0005) {
          int pln        = ch->plane();
          int pnl        = ch->panel();
          int pnl216     = pln*6+pnl;
          fPanelDt[pnl216]  = ch->time-tc_best->t0;
          fPanelNh[pnl216] += 1;
        }
      }
      // average times
      for (int i=0; i<216; i++) {
        fPanelDt[i] = fPanelDt[i]/(fPanelNh[i]+1.e-12);
      }
    }

//-----------------------------------------------------------------------------
// prep done, now fill non-residual histograms
//-----------------------------------------------------------------------------
    FillHistograms();
  }
  
//-----------------------------------------------------------------------------
// post-loop
//-----------------------------------------------------------------------------
  // for (int ip=0; ip<216; ip++) {
  //   fHist.h_panel_dt_111_vs_evn[ip]->GetXaxis()->SetRangeUser(0,fMaxEvent+100);
  // }
}

//-----------------------------------------------------------------------------
int plot_n002_tc::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n002_tc::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
// one occupancy canvas per station Ip2>Ip1
//-----------------------------------------------------------------------------
int plot_n002_tc::FitHistogram(TH1F* Hist, fit_result_t* Fr, int Ip1, int Ip2, int NMin) {

  //  fit_result_t* fr = &fFr[Ip2][Ip1];

  Fr->chi2dof = -1;
      
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = 0;
    Fr->e[ip] = -1;
  }

  // TH1F* h = fHist->h_dt05[Ip2][Ip1];

  int nent = Hist->GetEntries();
      
  if (nent < NMin) {
    return -1;
  }

  // find max bin

  int nbins = Hist->GetNbinsX();
  int imax = -1;
  int qmax = -1;
  for (int i=0; i<nbins; i++) {
    float q = Hist->GetBinContent(i+1);
    if (q > qmax) {
      imax = i+1;
      qmax = q;
    }
  }

  float tmax = Hist->GetBinCenter(imax);
  
  // for some reason, "sq" is required for tfr be defined
  TFitResultPtr tfr = Hist->Fit("gaus","sq","",tmax-50,tmax+50);
  
  if ((! tfr->IsValid()) or tfr->IsEmpty()) {
    std::cout << std::format("# FIT ERROR: h_dt05[{}]:[{}]\n",Ip2,Ip1);
    return -1;
  }

  Fr->chi2dof = tfr->Chi2()/tfr->Ndf();
  double sf    = sqrt(Fr->chi2dof);
          
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = tfr->Parameter(ip);
    Fr->e[ip] = tfr->Error(ip)*sf;
  }
  return 0;
}


//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n002_tc::PrintHistograms(int Refit) {

  gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_n002_tc.pdf",fRunNumber);
  
  gStyle->SetStatW(0.30);   // wider (NDC)


  for (int ip1=0; ip1<35; ip1++) {
    for (int ip2=ip1+1; ip2<36; ip2++) {

      //      fit_result_t* fr = &fFr[ip2][ip1];

      if (Refit) FitHistogram(fHist->h_dt05[ip2][ip1],&fFr[ip2][ip1],ip2,ip1);
    }
  }


  // next: perform fits
  

  for (int ic=0; ic<3; ic++) {
    // ic : canvas index. Plot 12 distributions per canvas
    
    TCanvas c(Form("c_%02i",ic),Form("c_%02i",ic),1600,1800);
    c.Divide(3,4);
    for (int ip=0; ip<12; ip++) {
      
      // 12 plane plots per canvas
      int plane = ic*12+ip;
      if (plane == 0) continue;

      // by default, delta_t = dt05[plane]-dt05[plane-1]
      // dt05 : hits above 0.5 keV
      
      TH1F* h = fHist->h_dt05[plane][plane-1];
      h->GetXaxis()->SetRangeUser(-100,100);
      c.cd(ip+1);
      h->Draw();

      // float hmax = (int(fMaxEvent/1.e6)+1)*1e6;
      //      // normalization to the rate :

      //      float input_rate = 1.e4;   // 10 kHz
      //      float scale = input_rate/(fNEvents+1.e-12);
      //      h->Scale(scale);
      //      h->SetMaximum(hmax);
      //      gPad->SetLogy(kTRUE);

    // make statbox transparent
      gPad->Update();           // create stats box

      auto st = (TPaveStats*)h->FindObject("stats");
      if (st) {
        st->SetFillStyle(0);    // transparent
        st->SetBorderSize(1);   // optional
      }
      gPad->Modified();
      gPad->Update();
    }

    if (ic == 0) {
      c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
    }
    else if (ic == 2) {
      c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
    }
    else {
      c.Print(fn.data());     // or .pdf, .root, ...
    }
  }


  // finally, print the resulting table

  float corr[36];

  std::cout << std::format("# i       dt      corr[i]   fit_chi2\n");
  for (int ip=0; ip<36; ip++) {
    corr[ip] = -1;
    fit_result_t* fr = &fFr[ip][ip-1];
    float dt = 0;

    if (ip == 0) {
      corr[ip] = 0;
    }
    else if (ip > 0) {
      dt = fr->p[1];
      corr[ip] = corr[ip-1]+dt; // total correction wrt plane 0
    }
    
    // planes 14 and 15 - slot 7 - are missing 

    if ((ip == 14) or (ip == 15)) {
      dt       = -1;
      corr[ip] = -1;
    }
    else if (ip == 16) {

      // pick the closest
      fr = &fFr[16][13];
      dt = fr->p[1];
      corr[ip] = corr[13]+dt;
      
    }

    // inverse signs of the step and total offsets (to add to calibrations)
    std::cout  << std::format(" {:3d} {:10.4f} {:10.4f} {:10.4f} \n",ip,-dt,corr[ip],fr->chi2dof);
    
  }

  return 0;
}
//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n002_tc::PrintPanelDt() {

  // gROOT->SetBatch(kTRUE);   // no GUI windows

  // std::string fn = std::format("run_{:6d}_n002_tc.pdf",fRunNumber);
  
  gStyle->SetStatW(0.30);   // wider (NDC)


  fHist->h_dt_vs_panel = new TH1F("dt_vs_panel","dt vs panel",216,0,216);
  
  for (int ip=0; ip<216; ip++) {
    FitHistogram(fHist->h_panel_dt[ip],&fPanelFr[ip],ip,0);
    float y  = fPanelFr[ip].p[1];
    float ey = fPanelFr[ip].e[1];
    fHist->h_dt_vs_panel->SetBinContent(ip+1,y );
    fHist->h_dt_vs_panel->SetBinError  (ip+1,ey);
    if (fPanelFr[ip].chi2dof > 0) {
      fHist->h_pdt_216->Fill(y);
    }
  }

  fHist->h_dt_vs_panel->Draw();

  // finally, print the resulting table


  std::cout << std::format("# i       dt      err      fit_chi2\n");
  for (int ip=0; ip<216; ip++) {
    //    corr[ip] = -1;
    fit_result_t* fr = &fPanelFr[ip];

    // inverse signs of the step and total offsets (to add to calibrations)
    std::cout  << std::format(" {:3d} {:10.4f} {:10.4f} {:10.4f} \n",ip,fr->p[1],fr->e[1],fr->chi2dof);
   
  }

  return 0;
}
