///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include "ana/plot_crv_tc_dt.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_crv_tc_dt"

//-----------------------------------------------------------------------------
plot_crv_tc_dt::plot_crv_tc_dt(int RunNumber, const char* Fn, const char* Label) :
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
plot_crv_tc_dt::~plot_crv_tc_dt() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookCrvcHistograms(CrvcHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d}",fRunNumber,Index->slot);
  std::string name, title;

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBook->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookCrvpHistograms(CrvpHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d}",fRunNumber,Index->slot);
  std::string name, title;

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBook->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  return 0;
}



//-----------------------------------------------------------------------------
int plot_crv_tc_dt::BookHistograms(TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  Index_t index;

  fHist = new Hist_t;

  int book_crvc_histset[10];
  int n_crvc_histsets(10);

  for (int i=0; i<n_crvc_histsets; i++) { book_crvc_histset[i] = 0; }

  book_crvc_histset[0] = 1;

  for (int i=0; i<n_crvc_histsets; i++) {
    if (book_crvc_histset[i] == 0) continue;
    std::string folder_name = std::format("crvc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    fHist->crvc[i] = new CrvcHist_t();
    BookCrvcHistograms(fHist->crvc[i],&index,fol);
  }

  int book_crvp_histset[10];
  int n_crvp_histsets(10);

  for (int i=0; i<n_crvp_histsets; i++) { book_crvp_histset[i] = 0; }

  book_crvp_histset[0] = 1;

  for (int i=0; i<n_crvp_histsets; i++) {
    if (book_crvp_histset[i] == 0) continue;
    std::string folder_name = std::format("crvp_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    fHist->crvp[i] = new CrvpHist_t();
    BookCrvpHistograms(fHist->crvp[i],&index,fol);
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
int plot_crv_tc_dt::FillHistograms() {
  // filling histograms: plot time differences between

  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;

  Index_t index;
  
  for (int i=0; i<fEvent->ntc; i++) {
    DaqTimeCluster*  tc = (DaqTimeCluster* ) fEvent->Tc(i);

    for (int i2=0; i2<fEvent->ncrvc; i2++) {
      DaqCrvCoincidenceCluster*  crvc = (DaqCrvCoincidenceCluster* ) fEvent->Crvc(i2);
      float dt = crvc->time-tc->t0;
      
      fHist->crvc[0]->h_dt->Fill(dt);
    }
  }
  
  for (int i=0; i<fEvent->ntc; i++) {
    DaqTimeCluster*  tc = (DaqTimeCluster* ) fEvent->Tc(i);

    for (int i2=0; i2<fEvent->ncrvc; i2++) {
      DaqCrvRecoPulse*  crvp = (DaqCrvRecoPulse* ) fEvent->Crvp(i2);
      float dt = crvp->time-tc->t0;
      
      fHist->crvp[0]->h_dt->Fill(dt);
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

