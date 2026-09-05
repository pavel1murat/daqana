///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include <iostream>
#include <fstream>

#include "ana/plot_calo_tc_dt.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_calo_tc_dt"

// #include "daqana/obj/DaqEvent.hh"

//-----------------------------------------------------------------------------
plot_calo_tc_dt::plot_calo_tc_dt(int RunNumber, const char* Fn) :
  TNamed(Form("run_%06d_plot_calo_tc_dt",RunNumber),Form("run_%06d_plot_calo_tc_dt",RunNumber)),
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
  
  fCalCm     = CalChannelMap_t::Instance();
  
  fBook      = new Booking(fRunFolder);

  fFrRef     = nullptr;

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  
  fHist = new Hist_t;
  BookHistograms(fHist,fRunFolder);
}


//-----------------------------------------------------------------------------
plot_calo_tc_dt::plot_calo_tc_dt(int RunNumber, int SubrunNumber, const char* Label) :
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

  fCalCm->Init(RunNumber);

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
plot_calo_tc_dt::~plot_calo_tc_dt() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_calo_tc_dt::BookCalcHistograms(CalcHist_t* Hist, CalIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  name  = "edep";
  title = std::format("{} : edep",prefix);
  fBook->HBook1F(Hist->h_edep,name.data(),title.data(),100,0,1000,Folder);   // in MeV

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBook->HBook1F(Hist->h_dt_tc,name.data(),title.data(),100,-100,100,Folder);   // in ns...

  name  = "dt_crvc";
  title = std::format("{} : dt CRVC",prefix);
  fBook->HBook1F(Hist->h_dt_crvc,name.data(),title.data(),100,-100,100,Folder);   // in ns...

  name  = "dt_crvc_vs_dt_tc";
  title = std::format("{} : dt CRVC vs dt TC",prefix);
  fBook->HBook2F(Hist->h_dt_crvc_vs_dt_tc,name.data(),title.data(),100,-100,100,100,-100,100,Folder);   // in ns...

  // name  = "feb";
  // name  = "npes";
  // title = std::format("{} : npes",prefix);
  // fBook->HBook1F(Hist->h_npes,name.data(),title.data(),100,0,500,Folder);   // in us...

  // name  = "time";
  // title = std::format("{} : time",prefix);
  // fBook->HBook1F(Hist->h_time,name.data(),title.data(),100,0,1.e5,Folder);   // in us...

  // title = std::format("{} : feb",prefix);
  // fBook->HBook1F(Hist->h_feb,name.data(),title.data(),100,0,100,Folder);   // in us...

  // name  = "ch";
  // title = std::format("{} : ch",prefix);
  // fBook->HBook1F(Hist->h_ch,name.data(),title.data(),2000,0,2000,Folder);   // in us...

  // name  = "dt_vs_feb";
  // title = std::format("{} : dt vs feb",prefix);
  // fBook->HBook2F(Hist->h_dt_vs_feb,name.data(),title.data(),100,0,100,1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
int plot_calo_tc_dt::BookCalhHistograms(CalhHist_t* Hist, CalIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  // name  = "ph";
  // title = std::format("{} : ph",prefix);
  // fBook->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...

  // name  = "npes";
  // title = std::format("{} : npes",prefix);
  // fBook->HBook1F(Hist->h_npes,name.data(),title.data(),100,0,500,Folder);   // in us...

  // name  = "time";
  // title = std::format("{} : time",prefix);
  // fBook->HBook1F(Hist->h_time,name.data(),title.data(),100,0,1.e5,Folder);   // in us...

  // name  = "dt";
  // title = std::format("{} : dt",prefix);
  // fBook->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  // name  = "feb";
  // title = std::format("{} : feb",prefix);
  // fBook->HBook1F(Hist->h_feb,name.data(),title.data(),100,0,100,Folder);   // in us...

  // name  = "ch";
  // title = std::format("{} : ch",prefix);
  // fBook->HBook1F(Hist->h_ch,name.data(),title.data(),2000,0,2000,Folder);   // in us...

  // name  = "dt_vs_feb";
  // title = std::format("{} : dt vs feb",prefix);
  // fBook->HBook2F(Hist->h_dt_vs_feb,name.data(),title.data(),100,0,100,1000,-1000,1000,Folder);   // in us...

  return 0;
}




//-----------------------------------------------------------------------------
// by FEB: 32 scintillation bars per FEB 
//-----------------------------------------------------------------------------
int plot_calo_tc_dt::BookDiskHistograms(DiskHist_t* Hist, CalIndex_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  // Index_t index;

  // std::string prefix = std::format("run:{:06d} roc:{:02d}",fRunNumber,Index->roc);
  // std::string name, title;

  // name  = "sbid";
  // title = std::format("{} : SBID",prefix);
  // fBook->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  // name  = "feb_vs_dt";
  // title = std::format("{} : FEB vs_dt",prefix);
  // fBook->HBook2F(Hist->h_feb_vs_dt,name.data(),title.data(),1000,-1000,1000,30,0,30,Folder);

  // name  = "feb_vs_ch";
  // title = std::format("{} : FEB vs CH",prefix);
  // fBook->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

// //-----------------------------------------------------------------------------
// // by FEB - 24 or 25 FEBs per ROC ?
// //-----------------------------------------------------------------------------
//   int n_feb_histsets(30);

//   for (int i=0; i<n_feb_histsets; i++) {
//     std::string folder_name = std::format("feb_{:02d}",i);
//     TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
//     if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
//     Hist->feb[i] = new FebHist_t();
//     Index->feb = i;
//     BookFebHistograms(Hist->feb[i],Index,fol);
//   }

  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_calo_tc_dt::BookTrkHistograms(TrkHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "nhits";
  title = std::format("{} : nhits",prefix);
  fBook->HBook1F(Hist->h_nhits,name.data(),title.data(),100,0,100,Folder);   // in MeV

  name  = "chi2d";
  title = std::format("{} : chi2/ndod",prefix);
  fBook->HBook1F(Hist->h_chi2d,name.data(),title.data(),100,0,20,Folder);   // in ns...

  name  = "t0";
  title = std::format("{} : t0",prefix);
  fBook->HBook1F(Hist->h_t0,name.data(),title.data(),500,0,2.5e6,Folder);   // in ns...

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBook->HBook1F(Hist->h_dt_tc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dt_calc";
  title = std::format("{} : dt CALC",prefix);
  fBook->HBook1F(Hist->h_dt_calc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dt_crvc";
  title = std::format("{} : dt CRVC",prefix);
  fBook->HBook1F(Hist->h_dt_crvc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dx_calc";
  title = std::format("{} : dx_calc = x(trk)-x(calc)",prefix);
  fBook->HBook1F(Hist->h_dx_calc,name.data(),title.data(),200,-1000,1000,Folder);   // in ns...

  name  = "dy_calc";
  title = std::format("{} : dy_calc = y(trk)-y(calc)",prefix);
  fBook->HBook1F(Hist->h_dy_calc,name.data(),title.data(),100,-2500,2500,Folder);   // in ns...

  name  = "dx_calc_vs_dxdz";
  title = std::format("{} : dx_calc vs dxdz",prefix);
  fBook->HBook2F(Hist->h_dx_calc_vs_dxdz,name.data(),title.data(),200,-1,1,200,-1000,1000,Folder);

  name  = "dy_calc_vs_dydz";
  title = std::format("{} : dy_calc vs dydz",prefix);
  fBook->HBook2F(Hist->h_dy_calc_vs_dydz,name.data(),title.data(),100,-2.5,2.5,100,-2500,2500,Folder);

  name  = "dxdz";
  title = std::format("{} : dxdz",prefix);
  fBook->HBook1F(Hist->h_dxdz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "dydz";
  title = std::format("{} : dydz",prefix);
  fBook->HBook1F(Hist->h_dydz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "xcrv";
  title = std::format("{} : xcrv",prefix);
  fBook->HBook1F(Hist->h_xcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "zcrv";
  title = std::format("{} : zcrv",prefix);
  fBook->HBook1F(Hist->h_zcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "dx_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBook->HBook1F(Hist->h_dx_crvc,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "dz_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBook->HBook1F(Hist->h_dz_crvc,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "dx_crvc_vs_dxdy";
  title = std::format("{} : dx_crvc vs dxdy",prefix);
  fBook->HBook2F(Hist->h_dx_crvc_vs_dxdy,name.data(),title.data(),100,-5,5,200,-10000,10000,Folder);

  name  = "dz_crvc_vs_dzdy";
  title = std::format("{} : dz_crvc vs dzdy",prefix);
  fBook->HBook2F(Hist->h_dz_crvc_vs_dzdy,name.data(),title.data(),100,-2,2,100,-500,500,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int plot_calo_tc_dt::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CalIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "dt_vs_sipmid";
  title = std::format("{} : dt vs crystal ID",prefix);
  fBook->HBook2F(Hist->h_dt_vs_sipmid,name.data(),title.data(),3000,0,3000,1000,-1000,1000,Folder);

  name  = "sipmid";
  title = std::format("{} : Sipm ID",prefix);
  fBook->HBook1F(Hist->h_sipmid,name.data(),title.data(),3000,0,3000,Folder);

  name  = "n2_vs_n1";
  title = std::format("{} : N1:N1 calh E>10",prefix);
  fBook->HBook2F(Hist->h_n2_vs_n1,name.data(),title.data(),20,0,20,20,0,20,Folder);

  name  = "ntrk_0";
  title = std::format("{} : ntrk[0]",prefix);
  fBook->HBook1F(Hist->h_ntrk[0],name.data(),title.data(),10,0,10,Folder);

  name  = "ntrk_1";
  title = std::format("{} : ntrk[1]",prefix);
  fBook->HBook1F(Hist->h_ntrk[1],name.data(),title.data(),10,0,10,Folder);
//-----------------------------------------------------------------------------
// calorimeter hits 
//-----------------------------------------------------------------------------
  int book_calh_histset[10];
  int n_calh_histsets(10);

  for (int i=0; i<n_calh_histsets; i++) { book_calh_histset[i] = 0; }

  book_calh_histset[0] = 1;             // all

  for (int i=0; i<n_calh_histsets; i++) {
    if (book_calh_histset[i] == 0) continue;
    std::string folder_name = std::format("calh_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calh[i] = new CalhHist_t();
    index.sel = i;
    BookCalhHistograms(Hist->calh[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CALC: calorimeter clusters 
//-----------------------------------------------------------------------------
  int book_calc_histset[100];
  int n_calc_histsets(100);

  for (int i=0; i<n_calc_histsets; i++) { book_calc_histset[i] = 0; }

  book_calc_histset[0] = 1;             // all
  book_calc_histset[1] = 1;             // all DISK 0
  book_calc_histset[2] = 1;             // all DISK 1
  book_calc_histset[3] = 1;             // |dt_tc| < 30, all 
  book_calc_histset[4] = 1;             // |dt_tc| < 30, DISK 0
  book_calc_histset[5] = 1;             // |dt_tc| < 30, DISK 1
  book_calc_histset[6] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, all
  book_calc_histset[7] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, DISK 0
  book_calc_histset[8] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, DISK 1

  int sel_disk[100];
  sel_disk[0] = -1;
  sel_disk[1] =  0;
  sel_disk[2] =  1;
  sel_disk[3] = -1;
  sel_disk[4] =  0;
  sel_disk[5] =  1;
  sel_disk[6] = -1;
  sel_disk[7] =  0;
  sel_disk[8] =  1;

  for (int i=0; i<n_calc_histsets; i++) {
    if (book_calc_histset[i] == 0) continue;
    std::string folder_name = std::format("calc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calc[i] = new CalcHist_t();
    index.sel  = i;
    index.disk = sel_disk[i];
    BookCalcHistograms(Hist->calc[i],&index,fol);
  }
//-----------------------------------------------------------------------------
// Trk: tracks 
//-----------------------------------------------------------------------------
  TrkIndex_t trk_index;

  // std::string prefix = std::format("run:{:06d}",fRunNumber);
  int book_trk_histset[kMaxTrkHistSets];

  for (int i=0; i<kMaxTrkHistSets; i++) { book_trk_histset[i] = 0; }

  book_trk_histset[0] = 1;             // all
  
  book_trk_histset[1] = 1;             // in-time tracks
  book_trk_histset[2] = 1;             // in-time tracks DISK0
  book_trk_histset[3] = 1;             // in-time tracks DISK1

  book_trk_histset[4] = 1;             // in-time CAL
  book_trk_histset[5] = 1;             // in-time CAL DISK0
  book_trk_histset[6] = 1;             // in-time CAL DISK1

  book_trk_histset[7] = 1;             // in-time trk_nhits>10 CAL
  book_trk_histset[8] = 1;             // in-time trk_nhits>10 CAL DISK0
  book_trk_histset[9] = 1;             // in-time trk_nhits>10 CAL DISK1

  book_trk_histset[10] = 1;             // in-time CRVC trk_nhits>10

  for (int i=0; i<kMaxTrkHistSets; i++) {
    if (book_trk_histset[i] == 0) continue;
    std::string folder_name = std::format("trk_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->trk[i]   = new TrkHist_t();
    trk_index.sel  = i;
    BookTrkHistograms(Hist->trk[i],&trk_index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_calo_tc_dt::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_calo_tc_dt::Init(TTree *tree) {
  
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
int plot_calo_tc_dt::FillCalcHistograms(CalcHist_t* Hist, DaqCaloCluster* Calc, calc_param_t* Cp) {
  // filling histograms: plot time differences between
  
  Hist->h_edep->Fill(Calc->edep);
  Hist->h_dt_tc->Fill(Cp->dtmin_tc);
  Hist->h_dt_crvc->Fill(Cp->dtmin_crvc);
  Hist->h_dt_crvc_vs_dt_tc->Fill(Cp->dtmin_tc,Cp->dtmin_crvc);
  
  // Hist->h_npes->Fill(Crvp->npes);
  // Hist->h_time->Fill(Crvp->time);
  // Hist->h_feb->Fill(Crvp->feb);
  // Hist->h_ch->Fill(Crvp->ch);
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int plot_calo_tc_dt::FillDiskHistograms(DiskHist_t* Hist, DaqCaloRecoDigi* Calrd) {
  // filling histograms: plot time differences between
  // Hist->h_ph->Fill(Crvp->ph);
  // Hist->h_npes->Fill(Crvp->npes);
  // Hist->h_time->Fill(Crvp->time);
  // Hist->h_feb->Fill(Crvp->feb);
  // Hist->h_ch->Fill(Crvp->ch);
  return 0;
}

//-----------------------------------------------------------------------------
int plot_calo_tc_dt::FillTrkHistograms(TrkHist_t* Hist, DaqTrack* Trk, trk_param_t* Tp) {
  // filling histograms: plot time differences between
  
  Hist->h_nhits->Fill(Trk->nhits);
  Hist->h_chi2d->Fill(Trk->chi2/Trk->ndof);
  Hist->h_t0->Fill(Trk->t0);
  Hist->h_dt_tc->Fill(Tp->dtmin_tc);

  Hist->h_dt_calc->Fill(Tp->dtmin_calc);
  Hist->h_dt_crvc->Fill(Tp->dtmin_crvc);

  Hist->h_dx_calc->Fill(Tp->dx_calc);
  Hist->h_dy_calc->Fill(Tp->dy_calc);

  float dxdz = Trk->nx/Trk->nz;
  float dydz = Trk->ny/Trk->nz;
  
  Hist->h_dx_calc_vs_dxdz->Fill(dxdz,Tp->dx_calc);
  Hist->h_dy_calc_vs_dydz->Fill(dydz,Tp->dy_calc);

  float dxdy = dxdz/dydz;
  float dzdy = 1./dydz;

  Hist->h_dxdz->Fill(dxdz);
  Hist->h_dydz->Fill(dydz);
  
  Hist->h_xcrv->Fill(Tp->xcrv);
  Hist->h_zcrv->Fill(Tp->zcrv);

  Hist->h_dx_crvc->Fill(Tp->dx_crvc);
  Hist->h_dz_crvc->Fill(Tp->dz_crvc);

  Hist->h_dx_crvc_vs_dxdy->Fill(dxdy,Tp->dx_crvc);
  Hist->h_dz_crvc_vs_dzdy->Fill(dzdy,Tp->dz_crvc);
  
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int plot_calo_tc_dt::FillHistograms() {
  // filling histograms: plot time differences between

  //  Index_t index;
  
//-----------------------------------------------------------------------------
// double-nested loops start here
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fEvent->ncalrd; i1++) {
    DaqCaloRecoDigi*  calrd = fEvent->Calrd(i1);

    // find closest time cluster
    float dtmin = 1.e6;
    for (int i2=0; i2<fEvent->ntc; i2++) {
      DaqTimeCluster*  tc = fEvent->Tc(i2);
      if (tc->nsh < 8) continue;

      float dt       = calrd->time-tc->t0;
      if (fabs(dt) < dtmin) {
        dtmin = dt;
      }
    }
    fHist->h_dt_vs_sipmid->Fill(calrd->sipmid,dtmin);
    fHist->h_sipmid->Fill(calrd->sipmid);
  }
  
  fHist->h_n2_vs_n1->Fill(fNCalh10[0],fNCalh10[1]);

//-----------------------------------------------------------------------------
// fill cluster histograms
//-----------------------------------------------------------------------------
  int n_good_tc = 0;
  
  for (int i=0; i<fEvent->ncalc; i++) {
    DaqCaloCluster* calc     = fEvent->Calc(i);
    calc_param_t*   calc_par = &fListOfCalcParam[i];
    FillCalcHistograms(fHist->calc[0],calc,calc_par);
    if (calc->disk == 0) FillCalcHistograms(fHist->calc[1],calc,calc_par);
    else                 FillCalcHistograms(fHist->calc[2],calc,calc_par);

    if (fabs(calc_par->dtmin_tc) < 30) {
      FillCalcHistograms(fHist->calc[3],calc,calc_par);
      if (calc->disk == 0) FillCalcHistograms(fHist->calc[4],calc,calc_par);
      else                 FillCalcHistograms(fHist->calc[5],calc,calc_par);

      if (fabs(calc_par->dtmin_crvc) < 30) {
        FillCalcHistograms(fHist->calc[6],calc,calc_par);
        if (calc->disk == 0) FillCalcHistograms(fHist->calc[7],calc,calc_par);
        else                 FillCalcHistograms(fHist->calc[8],calc,calc_par);

        if (calc->edep > 20) {
          n_good_tc += 1;
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// fill track histograms
//-----------------------------------------------------------------------------
  fHist->h_ntrk[0]->Fill(fEvent->ntrk);
  if (n_good_tc > 0) {
    fHist->h_ntrk[1]->Fill(fEvent->ntrk);
  }

  for (int i=0; i<fEvent->ntrk; i++) {
    
   
    DaqTrack* trk = fEvent->Trk(i);
    trk_param_t* tp = &fListOfTrkParam[i];
    FillTrkHistograms(fHist->trk[0],trk,tp);
    if (tp->intime) {
      
      FillTrkHistograms(fHist->trk[1],trk,tp);
      if (tp->calc->disk == 0) {
        std::cout << std::format("in-time disk0: run:srn:evt : {:6}:{:06}:{}\n",fEvent->run,fEvent->srn,fEvent->evn);
        FillTrkHistograms(fHist->trk[2],trk,tp);
      }
      else {
        std::cout << std::format("in-time disk1: run:srn:evt : {:6}:{:06}:{}\n",fEvent->run,fEvent->srn,fEvent->evn);
        FillTrkHistograms(fHist->trk[3],trk,tp);
      }
    }
    
    if (tp->intime_calc) {
      FillTrkHistograms(fHist->trk[4],trk,tp);
      if (tp->calc->disk == 0) {
        FillTrkHistograms(fHist->trk[5],trk,tp);
      }
      else if (tp->calc->disk == 1) {
        FillTrkHistograms(fHist->trk[6],trk,tp);
      }
    }

    if ((trk->nhits > 10) and (tp->intime_calc)) {
      FillTrkHistograms(fHist->trk[7],trk,tp);
      if (tp->calc->disk == 0) {
        FillTrkHistograms(fHist->trk[8],trk,tp);
      }
      else if (tp->calc->disk == 1) {
        FillTrkHistograms(fHist->trk[9],trk,tp);
      }
    }

    if ((trk->nhits > 10) and (tp->intime_crvc)) {
      FillTrkHistograms(fHist->trk[10],trk,tp);
    }    
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
Long64_t plot_calo_tc_dt::LoadTree(Long64_t entry) {
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
int plot_calo_tc_dt::CalculateMissingTrkParameters() {
//-----------------------------------------------------------------------------
// extra track parameters
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fEvent->ntrk; i1++) {
    DaqTrack*  trk = fEvent->Trk(i1);

    trk_param_t* tp = &fListOfTrkParam[i1];
//-----------------------------------------------------------------------------
// initialize parameter record
//-----------------------------------------------------------------------------
    tp->dtmin_tc    = 1.e6;
    tp->dtmin_crvc  = 1.e6;
    tp->dtmin_calc  = 1.e6;
    tp->tc          = nullptr;
    tp->calc        = nullptr;
    tp->crvc        = nullptr;
    tp->intime      = 0;
    tp->intime_calc = 0;
    tp->intime_crvc = 0;
    tp->dx_calc     = 1.e6;
    tp->dy_calc     = 1.e6;
    tp->xcrv        = 1.e6;
    tp->zcrv        = 1.e6;
    tp->dx_crvc     = 1.e6;
    tp->dz_crvc     = 1.e6;
//-----------------------------------------------------------------------------    
// find the closest time cluster - should always be there
//-----------------------------------------------------------------------------    
    for (int i2=0; i2<fEvent->ntc; i2++) {
      DaqTimeCluster* tc = fEvent->Tc(i2);
      float dt = trk->t0-tc->t0;
      if (fabs(dt) < fabs(tp->dtmin_tc)) {
        tp->dtmin_tc = dt;
        tp->tc       = tc;
      }
    }
//-----------------------------------------------------------------------------
// determine the closest CRV coincidence
// as the calibration used time clusters, look at the time cluster T0
//-----------------------------------------------------------------------------
    float crv_time_offset = 0.; // 21; // today
    
    for (int i2=0; i2<fEvent->ncrvc; i2++) {
      DaqCrvCoincidenceCluster* crvc = fEvent->Crvc(i2);
      float dt = tp->tc->t0-(crvc->time-crv_time_offset);
      if (fabs(dt) < fabs(tp->dtmin_crvc)) {
        tp->dtmin_crvc = dt;
        tp->crvc       = crvc;
      }
    }
    
    if (tp->crvc) {
      if (fabs(tp->dtmin_crvc) <  30) {
        tp->intime_crvc = 1;
      }
    
      float ycrv  = tp->crvc->y ; // -80.;//-4280; // for kicks
      float dy    = ycrv-trk->y0;
      
      tp->xcrv    = trk->x0+(trk->nx/trk->ny)*dy;
      tp->zcrv    = trk->z0+(trk->nz/trk->ny)*dy;
      tp->dx_crvc = tp->xcrv-tp->crvc->x;
      tp->dz_crvc = tp->zcrv-(tp->crvc->z-23000.);  // approx
    }
//------------------------------;-----------------------------------------------
// determine the closest calorimeter cluster
//-----------------------------------------------------------------------------
    for (int i2=0; i2<fEvent->ncalc; i2++) {
      DaqCaloCluster* calc = fEvent->Calc(i2);

      float dt = tp->tc->t0-calc->time;
      if (fabs(dt) < fabs(tp->dtmin_calc)) {
        tp->dtmin_calc = dt;
        tp->calc       = calc;
      }
    }
//-----------------------------------------------------------------------------
// extrapolate track to the closest cluster
//-----------------------------------------------------------------------------
    if (tp->calc) {
      double zc = 2383.; // 3560.; // guesswork for the disk0 position , not sure it is correct
      if (tp->calc->disk == 1) {
        zc = 3517;
      }
      float dz    = zc-trk->z0;
      float tx    = trk->x0+trk->nx/trk->nz*dz;
      float ty    = trk->y0+trk->nx/trk->nz*dz;
      tp->dx_calc = tx-tp->calc->x;
      tp->dy_calc = ty-tp->calc->y;
    }

    if (fabs(tp->dtmin_calc) <  30) {
      tp->intime_calc = 1;
    }
    
    if ((fabs(tp->dtmin_tc  ) <  30) and
        (fabs(tp->dtmin_calc) <  30) and
        (fabs(tp->dtmin_crvc) <  30) and
        (tp->calc->edep       >= 20) and
        (trk->nhits           >= 10)
        ) {
      tp->intime = 1;
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_calo_tc_dt::CalculateMissingParameters() {

  fNCalh10[0] = 0;
  fNCalh10[1] = 0;

  fListOfCalcParam.clear();
  if (fEvent->ncalc > 0) fListOfCalcParam.resize(fEvent->ncalc);
  
  fListOfTrkParam.clear();
  if (fEvent->ntrk > 0) fListOfTrkParam.resize(fEvent->ntrk);
  
  for (int i1=0; i1<fEvent->ncalh; i1++) {
    DaqCaloHit*  calh = fEvent->Calh(i1);

    if (calh->edep > 10) {
      int disk = calh->Disk();
      fNCalh10[disk] += 1;
    }
  }
 
//-----------------------------------------------------------------------------
// extra parameters of the calorimeter clusters
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fEvent->ncalc; i1++) {
    DaqCaloCluster*  calc = fEvent->Calc(i1);

    calc_param_t& cp = fListOfCalcParam[i1];
//-----------------------------------------------------------------------------
// determine the closest time cluster
//-----------------------------------------------------------------------------
    cp.dtmin_tc = 1.e6;
    cp.tc       = nullptr;
    // determine the closest time cluster
    for (int i2=0; i2<fEvent->ntc; i2++) {
      DaqTimeCluster* tc = fEvent->Tc(i2);
      float dt = calc->time-tc->t0;
      if (fabs(dt) < fabs(cp.dtmin_tc)) {
        cp.dtmin_tc = dt;
        cp.tc       = tc;
      }
    }
//-----------------------------------------------------------------------------
// determine the closest CRV coincidence
//-----------------------------------------------------------------------------
    cp.dtmin_crvc = 1.e6;
    cp.crvc       = nullptr;

    float crv_time_offset = 0; // 21. // today
    
    for (int i2=0; i2<fEvent->ncrvc; i2++) {
      DaqCrvCoincidenceCluster* crvc = fEvent->Crvc(i2);
      float dt = calc->time-(crvc->time-crv_time_offset);
      if (fabs(dt) < fabs(cp.dtmin_crvc)) {
        cp.dtmin_crvc = dt;
        cp.crvc       = crvc;
      }
    }
//------------------------------;-----------------------------------------------
// determine the closest track
//-----------------------------------------------------------------------------
    cp.dtmin_trk = 1.e6;
    cp.trk       = nullptr;

    for (int i2=0; i2<fEvent->ntrk; i2++) {
      DaqTrack* trk = fEvent->Trk(i2);
      
      float dt = calc->time-trk->t0;
      if (fabs(dt) < fabs(cp.dtmin_trk)) {
        cp.dtmin_trk = dt;
        cp.trk       = trk;
      }
    }
  }

  CalculateMissingTrkParameters();

  return 0;
}

//-----------------------------------------------------------------------------
void plot_calo_tc_dt::Loop(int NEvents) {

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
// calculate missing parameters
//-----------------------------------------------------------------------------
    CalculateMissingParameters();
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
int plot_calo_tc_dt::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_calo_tc_dt::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
// Ip1, Ip2 
//-----------------------------------------------------------------------------
int plot_calo_tc_dt::FitHistogram(TH1* Hist, fit_result_t* Fr, float XMin, float XMax, int MinSum) {

  //  fit_result_t* fr = &fFr[Ip2][Ip1];

  Fr->chi2dof = -1;
      
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = 0;
    Fr->e[ip] = -1;
  }

  // TH1F* h = fHist->h_dt05[Ip2][Ip1];

  int nbins     = Hist->GetNbinsX();
  int integral  = Hist->Integral(1,nbins);
      
  if (integral < MinSum) {
    std::cout << std::format("ERROR: integral:{} < {}. BAIL OUT\n",integral,MinSum);
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

  if (qmax < 3) {
    std::cout << std::format("ERROR: qmax:{} < 3. BAIL OUT\n",qmax);
    return -2;
  }

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

  if (sum < MinSum) {
    std::cout << std::format("ERROR: sum:{} < {}. BAIL OUT\n",sum,MinSum);
    return -3;
  }

  float t0 = Hist->GetBinCenter(imax);
  
  // for some reason, "sq" is required for tfr be defined
  float tmin{t0-50}, tmax{t0+50};
  if (XMax > XMin) {
    tmin = XMin;
    tmax = XMax;
  }
  
  TFitResultPtr tfr = Hist->Fit("gaus","sql","",tmin,tmax);
  
  if ((! tfr->IsValid()) or tfr->IsEmpty()) {
    // assume tha all indices are in the name/title
    std::cout << std::format("# FIT ERROR: Hist->name:{} Hist->title:{}\n",Hist->GetName(),Hist->GetTitle());
    return -4;
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
int plot_calo_tc_dt::FitTimeOffsets(float TMin, float TMax) {

  TH2F* h2 = fHist->h_dt_vs_sipmid;

  int nch = kMaxNChannels;
  for (int i=0; i<nch; i++) {
    
    std::string hpname = std::format("hpx_{:04d}",i);
    TH1D* hp = h2->ProjectionY(hpname.data(),i+1,i+1);
    fit_result_t* fr = &fFr[i];
    FitHistogram(hp,fr,TMin,TMax,10);
  }
  
  // done fitting, write results

  std::ofstream os("CaloTimeCalib_corr.txt");

  os << std::format("TABLE CalTimeCalib {}\n",fRunNumber);
  os << std::format("# corrections to  CalTimeCalib cid= XXX, to be subtracted \n");
  os << std::format("# sipmid     dT      sigma      chi2 \n");
  
  for (int i=0; i<nch; i++) {
    
    // to do no harm if no fit, dt should be initialized to zero 
    float dt{0}; 
    if (fFr[i].chi2dof > 0) {
      dt = fFr[i].p[1];
    }

    os << std::format("{:5}   {:8.3f}  {:8.3f}  {:8.3f}\n",i,dt,fFr[i].p[2],fFr[i].chi2dof);
  }
  os.close();
  
  return 0;
}
