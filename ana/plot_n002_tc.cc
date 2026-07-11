///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include "ana/plot_n002_tc.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_n002_tc"

//-----------------------------------------------------------------------------
plot_n002_tc::plot_n002_tc(int RunNumber, const char* Fn) :
  TNamed(Form("run_%06d_n002_tc",RunNumber),Form("run_%06d_n002_tc",RunNumber)),
  fChain(0) {

  std::string dir("/data/mu2e/mu2etrk/datasets/vst00s000r000n002"); 
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000n002.{:06d}_000001.root",dir,RunNumber);
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

  Index_t index;
  
  for (int i=0; i<fEvent->ntc; i++) {
    DaqTimeCluster*  tc = (DaqTimeCluster* ) fEvent->tc->UncheckedAt(i);

    FillTimeClusterHistograms(fHist->tc[0],tc);
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
    
    // the number of straw hits is the same as teh number of straw digis
    for (int i=0; i<fEvent->nshtot; i++) {
      DaqStrawHit*  sh = (DaqStrawHit* ) fEvent->sh->UncheckedAt(i);
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
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n002_tc::PrintHistograms(int ISet) {

 //  gROOT->SetBatch(kTRUE);   // no GUI windows

 //  std::string fn = std::format("run_{:6d}_set_{:02}_occup.pdf",fRunNumber,ISet);
  
 //  gStyle->SetStatW(0.30);   // wider (NDC)

 // for (int is=0; is<18; is++) {
 //    TCanvas c(Form("c_%02i",is),Form("c_%02i",is),1600,1800);
 //    c.Divide(3,4);
 //    for (int ip=0; ip<12; ip++) {
 //      TH1F* h = fHist[ISet]->slot[is]->panel[ip]->h_occup;
 //      c.cd(ip+1);
 //      float hmax = (int(fMaxEvent/1.e6)+1)*1e6;
 //      // normalization to the rate :

 //      float input_rate = 1.e4;   // 10 kHz
 //      float scale = input_rate/(fNEvents+1.e-12);
 //      h->Scale(scale);
 //      h->SetMaximum(hmax);
 //      gPad->SetLogy(kTRUE);
 //      h->Draw("hist");
 //      // make statbox transparent
 //      gPad->Update();           // create stats box

 //      auto st = (TPaveStats*)h->FindObject("stats");
 //      if (st) {
 //        st->SetFillStyle(0);    // transparent
 //        st->SetBorderSize(1);   // optional
 //      }
 //      gPad->Modified();
 //      gPad->Update();
 //    }

 //    if (is == 0) {
 //      c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
 //    }
 //    else if (is == 17) {
 //      c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
 //    }
 //    else {
 //      c.Print(fn.data());     // or .pdf, .root, ...
 //    }
 //  }

  return 0;
}
