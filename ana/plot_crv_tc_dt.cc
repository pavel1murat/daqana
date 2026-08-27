///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include <iostream>
#include <fstream>

#include "ana/plot_crv_tc_dt.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_crv_tc_dt"

//-----------------------------------------------------------------------------
plot_crv_tc_dt::plot_crv_tc_dt(int RunNumber, const char* Fn) :
  TNamed(Form("run_%06d_plot_crv_tc_dt",RunNumber),Form("run_%06d_plot_crv_tc_dt",RunNumber)),
  fChain(0) {

  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
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
  
  fCcm       = CrvChannelMap_t::Instance(RunNumber);
  
  fBook      = new Booking(fRunFolder);

  fFrRef     = nullptr;

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  
  fHist = new Hist_t;
  BookHistograms(fHist,fRunFolder);
}


//-----------------------------------------------------------------------------
plot_crv_tc_dt::plot_crv_tc_dt(int RunNumber, int SubrunNumber, const char* Label) :
  TNamed(Form("run_%06d_n002_tc",RunNumber),Form("run_%06d_%s_tc",RunNumber,Label)),
  fChain(0) {

  std::string dir = std::format("/data/mu2e/mu2etrk/datasets/vst00s000r000{}",Label); 
  
  TFile *f(nullptr);

  std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000{}.{:06d}_{:06d}.root",dir,Label,RunNumber,SubrunNumber);
  f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
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
  fFrRef     = nullptr;

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  fHist = new Hist_t;
  BookHistograms(fHist,fRunFolder);
}


//-----------------------------------------------------------------------------
plot_crv_tc_dt::~plot_crv_tc_dt() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookCrvcHistograms(CrvcHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBook->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookCrvpHistograms(CrvpHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} roc:{:02d} feb:{}",
                                   fRunNumber,Index->sel, Index->roc, Index->feb);
  std::string name, title;

  name  = "ph";
  title = std::format("{} : ph",prefix);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...

  name  = "npes";
  title = std::format("{} : npes",prefix);
  fBook->HBook1F(Hist->h_npes,name.data(),title.data(),100,0,500,Folder);   // in us...

  name  = "time";
  title = std::format("{} : time",prefix);
  fBook->HBook1F(Hist->h_time,name.data(),title.data(),100,0,1.e5,Folder);   // in us...

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBook->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  name  = "feb";
  title = std::format("{} : feb",prefix);
  fBook->HBook1F(Hist->h_feb,name.data(),title.data(),100,0,100,Folder);   // in us...

  name  = "ch";
  title = std::format("{} : ch",prefix);
  fBook->HBook1F(Hist->h_ch,name.data(),title.data(),2000,0,2000,Folder);   // in us...

  name  = "dt_vs_feb";
  title = std::format("{} : dt vs feb",prefix);
  fBook->HBook2F(Hist->h_dt_vs_feb,name.data(),title.data(),100,0,100,1000,-1000,1000,Folder);   // in us...

  return 0;
}



//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookFebHistograms(FebHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  // Index_t index;

  std::string prefix = std::format("run:{:06d} roc:{} feb:{:02d}",fRunNumber,Index->roc, Index->feb);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBook->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  // name  = "dt";
  // title = std::format("{} : T(pulse)-T(trk TC)",prefix);
  // fBook->HBook1F(Hist->h_dt,name.data(),title.data(),400,-1000,1000,Folder);

  name  = "ch_vs_dt";
  title = std::format("{} : channel vs [T(pulse)-T(trk TC)]",prefix);
  fBook->HBook2F(Hist->h_ch_vs_dt,name.data(),title.data(),400,-1000,1000,64,0,64,Folder);

  // name  = "feb_vs_ch";
  // title = std::format("{} : FEB vs CH",prefix);
  // fBook->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);
  
  return 0;
}

//-----------------------------------------------------------------------------
// by FEB: 32 scintillation bars per FEB 
//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookRocHistograms(RocHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  // Index_t index;

  std::string prefix = std::format("run:{:06d} roc:{:02d}",fRunNumber,Index->roc);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBook->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  name  = "feb_vs_dt";
  title = std::format("{} : FEB vs_dt",prefix);
  fBook->HBook2F(Hist->h_feb_vs_dt,name.data(),title.data(),1000,-1000,1000,30,0,30,Folder);

  // name  = "feb_vs_ch";
  // title = std::format("{} : FEB vs CH",prefix);
  // fBook->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);

//-----------------------------------------------------------------------------
// by FEB - 24 or 25 FEBs per ROC ?
//-----------------------------------------------------------------------------
  int n_feb_histsets(30);

  for (int i=0; i<n_feb_histsets; i++) {
    std::string folder_name = std::format("feb_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->feb[i] = new FebHist_t();
    Index->feb = i;
    BookFebHistograms(Hist->feb[i],Index,fol);
  }

  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CrvIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBook->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  name  = "feb_vs_ch";
  title = std::format("{} : FEB vs CH",prefix);
  fBook->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

  name  = "dt_vs_sbid";
  title = std::format("{} : dt vs SBID",prefix);
  fBook->HBook2F(Hist->h_dt_vs_sbid,name.data(),title.data(),1500,0,1500,200,0,1000,Folder);

  name  = "feb_vs_sbid_0";
  title = std::format("{} : dt vs SBID ROC=1",prefix);
  fBook->HBook2F(Hist->h_feb_vs_sbid[0],name.data(),title.data(),600,0,600,30,0,30,Folder);

  name  = "feb_vs_sbid_1";
  title = std::format("{} : dt vs SBID ROC=2",prefix);
  fBook->HBook2F(Hist->h_feb_vs_sbid[1],name.data(),title.data(),600,0,600,30,0,30,Folder);

  name  = "och";
  title = std::format("{} : offline channel ID",prefix);
  fBook->HBook1F(Hist->h_och,name.data(),title.data(),2400,0,2400,Folder);

  // good for now
//-----------------------------------------------------------------------------
// CRV coincidence clusters
//-----------------------------------------------------------------------------
  int book_crvc_histset[10];
  int n_crvc_histsets(10);

  for (int i=0; i<n_crvc_histsets; i++) { book_crvc_histset[i] = 0; }

  book_crvc_histset[0] = 1;

  for (int i=0; i<n_crvc_histsets; i++) {
    if (book_crvc_histset[i] == 0) continue;
    std::string folder_name = std::format("crvc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvc[i] = new CrvcHist_t();
    BookCrvcHistograms(Hist->crvc[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  int book_crvp_histset[10];
  int n_crvp_histsets(10);

  for (int i=0; i<n_crvp_histsets; i++) { book_crvp_histset[i] = 0; }

  book_crvp_histset[0] = 1;             // all
  book_crvp_histset[1] = 1;             // ntc=1, dt_530 < 30
  book_crvp_histset[2] = 1;             // ntc=1, dt_590 < 30

  for (int i=0; i<n_crvp_histsets; i++) {
    if (book_crvp_histset[i] == 0) continue;
    std::string folder_name = std::format("crvp_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvp[i] = new CrvpHist_t();
    index.sel = i;
    BookCrvpHistograms(Hist->crvp[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// by ROC, links 0 and 3 --> rocs #1 and #4
//-----------------------------------------------------------------------------
  int book_roc_histset[10];
  int n_roc_histsets(10);

  for (int i=0; i<n_roc_histsets; i++) { book_roc_histset[i] = 0; }

  book_roc_histset[1] = 1;
  book_roc_histset[2] = 1;
  book_roc_histset[4] = 1;

  for (int i=0; i<n_roc_histsets; i++) {
    if (book_roc_histset[i] == 0) continue;
    std::string folder_name = std::format("roc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->roc[i] = new RocHist_t();
    index.roc = i;
    BookRocHistograms(Hist->roc[i],&index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_crv_tc_dt::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_crv_tc_dt::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int plot_crv_tc_dt::FillCrvpHistograms(CrvpHist_t* Hist, DaqCrvRecoPulse* Crvp) {
  // filling histograms: plot time differences between
  Hist->h_ph->Fill(Crvp->ph);
  Hist->h_npes->Fill(Crvp->npes);
  Hist->h_time->Fill(Crvp->time);
  Hist->h_feb->Fill(Crvp->feb);
  Hist->h_ch->Fill(Crvp->ch);
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int plot_crv_tc_dt::FillHistograms() {
  // filling histograms: plot time differences between

  //  Index_t index;
  
  for (int i2=0; i2<fEvent->ncrvp; i2++) {
    DaqCrvRecoPulse*  crvp = fEvent->Crvp(i2);

    // this is a global histogram
    fHist->h_feb_vs_ch->Fill(crvp->ch,crvp->feb);
                                        // have two histograms - one per roc , to color them
                                        // ROCs 1 and 2 --> hists 0 and 1
    fHist->h_feb_vs_sbid[crvp->roc-1]->Fill(crvp->sbid,crvp->feb);
    
    fHist->h_sbid->Fill(crvp->sbid);
    fHist->roc[crvp->roc]->h_sbid->Fill(crvp->sbid);
    fHist->roc[crvp->roc]->feb[crvp->feb]->h_sbid->Fill(crvp->sbid);
    int och = crvp->OfflineChID();
    fHist->h_och->Fill(och);
//-----------------------------------------------------------------------------
// CRVP[0] : all pulses
//-----------------------------------------------------------------------------
    FillCrvpHistograms(fHist->crvp[0],crvp);
  }
//-----------------------------------------------------------------------------
// double-nested loops start here
//-----------------------------------------------------------------------------
  for (int i=0; i<fEvent->ntc; i++) {
    DaqTimeCluster*  tc = fEvent->Tc(i);

    for (int i2=0; i2<fEvent->ncrvc; i2++) {
      DaqCrvCoincidenceCluster*  crvc = fEvent->Crvc(i2);
      float dt = crvc->time-tc->t0;
      
      fHist->crvc[0]->h_dt->Fill(dt);
    }

    for (int i2=0; i2<fEvent->ncrvp; i2++) {
      DaqCrvRecoPulse*  crvp = fEvent->Crvp(i2);
      float dt       = crvp->time-tc->t0;
      //int feb = crvp->feb;

      CrvpHist_t* hr = fHist->crvp[0];
      
      hr->h_dt->Fill(dt);
      fHist->h_dt_vs_sbid->Fill(crvp->sbid,dt);

      RocHist_t* roc_hr = fHist->roc[crvp->roc];
      roc_hr->h_feb_vs_dt->Fill(dt,crvp->feb);

      FebHist_t* feb_hr = fHist->roc[crvp->roc]->feb[crvp->feb];
      // feb_hr->h_dt->Fill(dt);
      feb_hr->h_ch_vs_dt->Fill(dt,crvp->ch);

      if (fEvent->ntc == 1) {
        if      (fabs(dt - 530) < 30) {
//-----------------------------------------------------------------------------
// CRVP[1] : first peak
//-----------------------------------------------------------------------------
          FillCrvpHistograms(fHist->crvp[1],crvp);
        }
        else if (fabs(dt - 590) < 30) {
//-----------------------------------------------------------------------------
// CRVP[2] : second peak
//-----------------------------------------------------------------------------
          FillCrvpHistograms(fHist->crvp[2],crvp);
        }
      }
    }
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
Long64_t plot_crv_tc_dt::LoadTree(Long64_t entry) {
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
void plot_crv_tc_dt::Loop(int NEvents) {

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
int plot_crv_tc_dt::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_crv_tc_dt::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
// Ip1, Ip2 
//-----------------------------------------------------------------------------
int plot_crv_tc_dt::FitHistogram(TH1* Hist, fit_result_t* Fr, float XMin, float XMax, int NMin) {

  //  fit_result_t* fr = &fFr[Ip2][Ip1];

  Fr->chi2dof = -1;
      
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = 0;
    Fr->e[ip] = -1;
  }

  // TH1F* h = fHist->h_dt05[Ip2][Ip1];

  int nbins     = Hist->GetNbinsX();
  int integral  = Hist->Integral(1,nbins);
      
  if (integral < NMin) {
    return -1;
  }

  // find max bin

  int   imax = -1;
  float qmax = -1;
  for (int i=0; i<nbins; i++) {
    float q = Hist->GetBinContent(i+1);
    if (q > qmax) {
      imax = i+1;
      qmax = q;
    }
  }

  if (qmax < 3) return -2;

  // estimate integral of the expected gaussian

  int i=0;
  //  int imin(0), imax(0);

  float sum = qmax;
  
  while (1) {
    i += 1;
                                        // check to the right of the maximum
    int iplus  = imax+i;
    if (iplus <= nbins) {
      float qplus = Hist->GetBinContent(iplus);
      if (qplus/qmax > 0.2) {
        sum += qplus;
      }
      else {
        // done
        iplus = nbins+1;
      }
    }
                                        // check the left side
    int iminus = imax-i;
    if (iminus > 0) {
                                        // bins start from 1
      
      float qminus = Hist->GetBinContent(iminus);
      if (qminus/qmax > 0.2) {
        sum += qminus;
      }
      else {
        // done
        iminus = -1;
      }
    }
    if ((iplus > nbins) and (iminus < 0)) break;
  }

  if (sum < 100) return -3;

  float t0 = Hist->GetBinCenter(imax);
  
  // for some reason, "sq" is required for tfr be defined
  float tmin{t0-30}, tmax{t0+30};
  if (XMax > XMin) {
    tmin = XMin;
    tmax = XMax;
  }
  
  TFitResultPtr tfr = Hist->Fit("gaus","sq","",tmin,tmax);
  
  if ((! tfr->IsValid()) or tfr->IsEmpty()) {
    // assume tha all indices are in the name/title
    std::cout << std::format("# FIT ERROR: Hist->name:{} Hist->title:{}\n",Hist->GetName(),Hist->GetTitle());
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
// do that for all ROCs and all FEBs
//-----------------------------------------------------------------------------
/*
# roc  feb   fit(crvp.time-tc.t0)  chi2dof
   1    1           530.724         3.004
   1    2           530.498         2.150
*/
int plot_crv_tc_dt::FitFebTimeOffsets(float TMin, float TMax) {

  for (int i=1; i<3; i++) {
    TH2F* h2 = fHist->roc[i]->h_feb_vs_dt;
    
    // fit Y-slices, use Ralf's integers
    for (int j=1; j<25; j++) {
      std::string hpname = std::format("hpx_{:02d}",j);
      TH1D* hp = h2->ProjectionX(hpname.data(),j+1,j+1);
      fit_result_t* fr = &fFr[i][j];
      FitHistogram(hp,fr,TMin,TMax,100);
    }
  }
  
  // done fitting, print results

  std::cout << std::format("# roc  feb   fit(crvp.time-tc.t0) dT(i-0)    sigma_i       chi2dof\n");
  
// find the first converged fit

  fFrRef = nullptr;

  bool ref_ch_found(false);
  
  for (int i=1; i<3; i++) {
    for (int j=1; j<25; j++) {
      if (fFr[i][j].chi2dof > 0) {
        fFrRef = &fFr[i][j];
        std::cout << std::format("reference channel: roc:{:2} feb:{:2} dt0:{:8.3f}\n",i,j,fFrRef->p[1]);
        ref_ch_found = true;
        break;
      }
    }
    if (ref_ch_found) break;
  }

  for (int i=1; i<3; i++) {
    for (int j=1; j<25; j++) {
      fit_result_t* fr = &fFr[i][j] ;
      float dt(0);
      if (fr->chi2dof > 0) {
        dt = fr->p[1]-fFrRef->p[1];
      }
      std::cout << std::format(" {:2d} {:4d}      {:9.3f}      {:9.3f}   {:9.3f}    {:9.3f}\n",
                               i,j,fr->p[1],dt,fr->p[2],fr->chi2dof);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// corrections are aimed to align in time all FEBs with FEB[1][1]
// to be called AFTER FitFebTimeOffset - that defines fFrRef
//------------------------------------------------------------------------------
int plot_crv_tc_dt::PrintTimeCorrections() {
  
  std::ofstream os("CrvTime_corr.txt");

  float dt0 = fFrRef->p[1];  // roc=1 feb=1

  for (int i=0; i<2304; i++) {
    
    CrvChannelMap_t::Data_t* dat = fCcm->ch_data_by_offline(i);

    float dt{0};
    if (dat != nullptr) {
      if (fFr[dat->roc][dat->feb].chi2dof > 0) {
        // dont correct FEBs with no fit
        dt = fFr[dat->roc][dat->feb].p[1] - dt0; // should be initialized to zero
      }
    }

    os << std::format("{:5}   {:8.3f}\n",i,dt);
  }

  os.close();
  
  return 0;
}
