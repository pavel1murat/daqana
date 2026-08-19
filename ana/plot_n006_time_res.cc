///////////////////////////////////////////////////////////////////////////////
//  plot timing residuals across the subsytems
/////////////////////////////////////////////////////////////////////////////////
#include "ana/plot_n006_time_res.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_n006_time_res"

//-----------------------------------------------------------------------------
plot_n006_time_res::plot_n006_time_res(int RunNumber, const char* Fn, const char* Label) :
  
  TNamed(Form("%s_occup",Label),Form("%s_occup",Label)), fChain(0) {

  std::string dir = std::format("/data/mu2e/mu2etrk/datasets/kpp00s000r000{}",Label); 
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/nts.mu2e.cosmics.kpp00s000r000{}.{:06d}_000001.root",dir,Label,RunNumber);
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
  RunDb::Instance()->GetRunInfo(RunNumber,&fRunInfo);

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
plot_n006_time_res::~plot_n006_time_res() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//-----------------------------------------------------------------------------
int plot_n006_time_res::BookSelHistograms(SelHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{} sel:{}",fRunNumber,Index->sel);
  std::string name, title;

  title = std::format("{}: ecal",prefix);
  fBook->HBook1F(Hist->h_ecal,"ecal",title.data(),100,0,10000,Folder);

  title = std::format("{}: edisk[0]",prefix);
  fBook->HBook1F(Hist->h_edisk[0],"edisk_0",title.data(),100,0,10000,Folder);

  title = std::format("{}: edisk[1]",prefix);
  fBook->HBook1F(Hist->h_edisk[1],"edisk_1",title.data(),100,0,10000,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n006_time_res::BookHistograms(TFolder* Folder) {

  std::string prefix = std::format("run {}",fRunNumber);
  std::string name, title;

  Index_t index;

  int book_sel_histset[10];
  int n_sel_histsets(10);

  for (int i=0; i<n_sel_histsets; i++) { book_sel_histset[i] = 0; }

  book_sel_histset[0] = 1;
  book_sel_histset[1] = 1;

  fHist = new Hist_t();

  title = std::format("{}: n close calo hits",prefix);
  fBook->HBook1F(fHist->h_n_close_calh,"n_close_calh",title.data(),20,0,20,Folder);

  title = std::format("{}: n close CRV CC",prefix);
  fBook->HBook1F(fHist->h_n_close_crvc,"n_close_crvc",title.data(),20,0,20,Folder);

  title = std::format("{}: dt_calh_tc",prefix);
  fBook->HBook1F(fHist->h_dt_calh_tc,"dt_calh_tc",title.data(),800,-2000,2000,Folder);

  title = std::format("{}: dt_crvc_tc[0]",prefix);
  fBook->HBook1F(fHist->h_dt_crvc_tc[0],"dt_crvc_tc[0]",title.data(),800,-2000,2000,Folder);

  title = std::format("{}: dt_crvc_tc[1]",prefix);
  fBook->HBook1F(fHist->h_dt_crvc_tc[1],"dt_crvc_tc[1]",title.data(),800,-2000,2000,Folder);

  for (int i=0; i<n_sel_histsets; i++) {
    if (book_sel_histset[i] == 0) continue;
    std::string folder_name = std::format("sel_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    fHist->sel[i]  = new SelHist_t;
    index.sel      =  i;
    index.slot     = -1;
    index.panel    = -1;
    index.mnid     = -1;
    index.ch       = -1;
    BookSelHistograms(fHist->sel[i],&index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n006_time_res::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n006_time_res::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
int plot_n006_time_res::FillSelHistograms(SelHist_t* Hist) {
  // filling histograms: plot time differences between

  Hist->h_ecal->Fill(fEvent->ecal);
  Hist->h_edisk[0]->Fill(fEvent->edisk[0]);
  Hist->h_edisk[1]->Fill(fEvent->edisk[1]);

  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n006_time_res::FillHistograms() {
  // filling histograms: plot time differences between


  fHist->h_n_close_calh->Fill(fNCloseCalh);
  fHist->h_n_close_crvc->Fill(fNCloseCrvc);
//-----------------------------------------------------------------------------
// sel[0]: all events
//-----------------------------------------------------------------------------
  FillSelHistograms(fHist->sel[0]);

//-----------------------------------------------------------------------------
// sel[1]: CAL and CRV timed in
//-----------------------------------------------------------------------------
  if ((fNCloseCalh > 0) and (fNCloseCrvc > 0)) {
    FillSelHistograms(fHist->sel[1]);
  }

  return 0;
}


//-----------------------------------------------------------------------------
Long64_t plot_n006_time_res::LoadTree(Long64_t entry) {
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
// load a given event
//-----------------------------------------------------------------------------
int plot_n006_time_res::LoadEvent(int Run, int Subrun, int Event) {
  int found = 0;
  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  for (int jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ((fEvent->run == Run) && (fEvent->srn == Subrun) and (fEvent->evn == Event)) {
      found = 1;
      break;
    }
  }

  return found;
}

//-----------------------------------------------------------------------------
void plot_n006_time_res::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;

  fNEvents  = nev;
  fMaxEvent = -1;

  //  int current_event = -1;
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
         
    if (fEvent->evn > fMaxEvent) {
      fMaxEvent = fEvent->evn;
    }

    fNCloseCrvc = 0;
    fNCloseCalh = 0;
    
    for (int itc=0; itc<fEvent->ntc; itc++) {
      DaqTimeCluster* tc = fEvent->Tc(itc);
      // select good TC
      if ((tc->nch < 8) or (tc->ngh < 2)) continue;
//-----------------------------------------------------------------------------
// prep done, now fill non-residual histograms
//-----------------------------------------------------------------------------
      for (int icalh=0; icalh<fEvent->ncalh; icalh++) {
        DaqCaloHit* calh = fEvent->Calh(icalh);
        float dt = calh->time - tc->t0;
        fHist->h_dt_calh_tc->Fill(dt);
        if (fabs(dt) < 100) {
          fNCloseCalh++;
        }
      }

      for (int icrvc=0; icrvc<fEvent->ncrvc; icrvc++) {
        DaqCrvCoincidenceCluster* crvc = fEvent->Crvc(icrvc);
        float dt = crvc->time-tc->t0;
        fHist->h_dt_crvc_tc[0]->Fill(dt);
        if (fNCloseCalh > 0) {
          fHist->h_dt_crvc_tc[1]->Fill(dt);
        }
        if (fabs(dt) < 100) {
          fNCloseCrvc++;
        }
      }
    }
    
    FillHistograms();
//-----------------------------------------------------------------------------
// post-loop
//-----------------------------------------------------------------------------
  // for (int ip=0; ip<216; ip++) {
  //   fHist.h_panel_dt_111_vs_evn[ip]->GetXaxis()->SetRangeUser(0,fMaxEvent+100);
  // }

    //   PrintNoisyChannels();
  }
}

//-----------------------------------------------------------------------------
int plot_n006_time_res::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n006_time_res::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n006_time_res::PrintHistograms(int ISet) {

  gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_set_{:02}_occup.pdf",fRunNumber,ISet);
  
  gStyle->SetStatW(0.30);   // wider (NDC)

  //  float hmax = (int(fMaxEvent/1.e4)+1)*1.e2;
  
  // for (int is=0; is<18; is++) {
  //   TCanvas c(Form("c_%02i",is),Form("c_%02i",is),1600,1800);
  //   c.Divide(3,4);
  //   for (int ip=0; ip<12; ip++) {
  //     TH1F* h = fHist[ISet]->slot[is]->panel[ip]->h_occup;
  //     c.cd(ip+1);
     
  //     // normalization to the rate in Hz:

  //     float total_time = fMaxEvent*fRunInfo.ew_length;    // assume no prescale
  //     float rejf       = fRunInfo.cfo_rate/fRunInfo.trigger_rate;
  //     float eff_time   = total_time/rejf;
  //     h->Scale(1./eff_time);
  //     h->SetMaximum(hmax);
  //     h->SetMinimum(0.5);
  //     gPad->SetLogy(kTRUE);
  //     h->Draw("hist");
  //     // make statbox transparent
  //     gPad->Update();           // create statbox

  //     auto st = (TPaveStats*)h->FindObject("stats");
  //     if (st) {
  //       st->SetFillStyle(0);    // transparent
  //       st->SetBorderSize(1);   // optional
  //     }
  //     gPad->Modified();
  //     gPad->Update();
  //   }

  //   if (is == 0) {
  //     c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
  //   }
  //   else if (is == 17) {
  //     c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
  //   }
  //   else {
  //     c.Print(fn.data());     // or .pdf, .root, ...
  //   }
  // }

  return 0;
}


