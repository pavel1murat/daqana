///////////////////////////////////////////////////////////////////////////////
// 
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include "ana/plot_n001_occup.hh"

#include "TPaveStats.h"
#include "TStyle.h"

#include "TRACE/tracemf.h"
#define TRACE_NAME "plot_n001_occup"

//-----------------------------------------------------------------------------
plot_n001_occup::plot_n001_occup(int RunNumber, const char* Fn, const char* Label) :
  
  TNamed(Form("%s_occup",Label),Form("%s_occup",Label)), fChain(0) {

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
plot_n001_occup::~plot_n001_occup() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n001_occup::BookChannelHistograms(ChannelHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run {:6d} slot:plane:panel:{:02d}:{:2d}:{:02d} MNID:{} channel:{:02d}",
                                   fRunNumber,Index->slot,Index->plane,Index->panel,Index->mnid,Index->ch);
  std::string name, title;

  name  = "ph";
  title = std::format("{} : pulse height",prefix);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),200,0,200,Folder);

  name  = "bl";
  title = std::format("{} : baseline",prefix);
  fBook->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,500,Folder);

  name  = "fs";
  title = std::format("{} : first sample",prefix);
  fBook->HBook1F(Hist->h_fs,name.data(),title.data(),50,0,50,Folder);

  name  = "tdc0";
  title = std::format("{} : tdc0, us",prefix);
  fBook->HBook1F(Hist->h_tdc0,name.data(),title.data(),1500,0,150,Folder);

  name  = "dt01";
  title = std::format("{} : dt01=tdc(CAL)-tdc(HV), ns",prefix);
  fBook->HBook1F(Hist->h_dt01,name.data(),title.data(),1000,-50,50,Folder);

  name  = "edep";
  title = std::format("{} : edep, keV",prefix);
  fBook->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n001_occup::BookPanelHistograms(PanelHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d} plane:{:02d} panel:{:03d} mnid:MN{:03d}",
                                   fRunNumber,Index->slot,Index->plane,Index->panel,Index->mnid);
  std::string name, title;

  name  = "occup";
  title = std::format("{} : occupancy",prefix);
  fBook->HBook1F(Hist->h_occup,name.data(),title.data(),100,0,100,Folder);

  name  = "edep";
  title = std::format("{} : edep, keV",prefix);
  fBook->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  for (int i=0; i<96; i++) {
    std::string folder_name = std::format("ch_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->ch[i] = new ChannelHist_t;
    Index->ch   = i;
    BookChannelHistograms(Hist->ch[i],Index,fol);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int plot_n001_occup::BookSlotHistograms(SlotHist_t* Hist, Index_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d}",fRunNumber,Index->slot);
  std::string name, title;

  name  = "edep";
  title = std::format("{} : edep",prefix);
  fBook->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  for (int i=0; i<12; i++) {
    std::string folder_name = std::format("pnl_{:03d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->panel[i] = new PanelHist_t;
    
    Index->plane = 2*Index->slot + i/6;
    Index->panel = i % 6;
    Index->pnl12 = i;
    TLOG(TLVL_DEBUG+1) << std::format("Index->plane:{:02d} Index->panel:{}",Index->plane,Index->panel);

    TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(Index->plane,Index->panel);
    Index->mnid    = tpmd->mnid;
    
    BookPanelHistograms(Hist->panel[i],Index,fol);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int plot_n001_occup::BookSelHistograms(Hist_t* Hist, Index_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  std::string prefix = std::format("run:{:06d} selection:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "occup_2d";
  title = std::format("{} : occupancy, panel vs ch",prefix);
  fBook->HBook2F(Hist->h_occup_2d,name.data(),title.data(),100,0,100,216,0,216,Folder);

  name  = "occup";
  title = std::format("{} : per panel occupancy",prefix);
  fBook->HBook1F(Hist->h_occup,name.data(),title.data(),216,0,216,Folder);

  Index_t index;

  for (int i=1; i<36; i++) {
    for (int j=i-1; j>=0; j--) {
      name  = std::format("dt05_{:02}_{:02}",i,j);
      title = std::format("{} : T{:02} - T{:02}",prefix,i,j);
      fBook->HBook1F(Hist->h_dt05[i][j],name.data(),title.data(),1000,-2500,2500,Folder);
    }
  }
  
  for (int i=0; i<18; i++) {
    std::string folder_name = std::format("slot_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->slot[i] = new SlotHist_t;
    index.slot  =  i;
    index.panel = -1;
    index.mnid  = -1;
    index.ch    = -1;
    BookSlotHistograms(Hist->slot[i],&index,fol);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n001_occup::BookHistograms(TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  Index_t index;

  int book_histset[10];
  int n_histsets(10);

  for (int i=0; i<n_histsets; i++) { book_histset[i] = 0; }

  book_histset[0] = 1;
  book_histset[1] = 1;

  for (int i=0; i<n_histsets; i++) {
    if (book_histset[i] == 0) continue;
    std::string folder_name = std::format("hist_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    fHist[i]    = new Hist_t;
    index.sel   =  i;
    index.slot  = -1;
    index.panel = -1;
    index.mnid  = -1;
    index.ch    = -1;
    BookSelHistograms(fHist[i],&index,fol);
  }
  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n001_occup::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n001_occup::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
int plot_n001_occup::FillChannelHistograms(ChannelHist_t* Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh) {
  Hist->h_ph->Fill(Sd->ph);
  Hist->h_bl->Fill(Sd->bl);
  Hist->h_fs->Fill(Sd->fs);

  float tdc0_ns = Sd->tdc0*5./256;      // CAL
  float tdc1_ns = Sd->tdc1*5./256;      // HV
  float dt01    = tdc0_ns-tdc1_ns;
  float tdc0_us = tdc0_ns/1000.;
  
  Hist->h_tdc0->Fill(tdc0_us);
  Hist->h_dt01->Fill(dt01);

  if (Sh) {
    Hist->h_edep->Fill(Sh->edep);
  }
                                        // no edep for the moment
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n001_occup::FillPanelHistograms(PanelHist_t* Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh) {
  Hist->h_occup->Fill(Index->ch);
  if (Sh) {
    Hist->h_edep->Fill(Sh->edep);
  }
  return 0;
}

//-----------------------------------------------------------------------------
// there are ntuples w/o straw hits (n007)
//-----------------------------------------------------------------------------
int plot_n001_occup::FillSlotHistograms(SlotHist_t* Hist, Index_t* Index, DaqStrawDigi* Sd, DaqStrawHit* Sh) {

  if (Sh) {
    Hist->h_edep->Fill(Sh->edep);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n001_occup::FillHistograms() {
  // filling histograms: plot time differences between

  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;

  Index_t index;

  int nsd = fEvent->sd->GetEntriesFast();
  for (int i=0; i<nsd; i++) {
    DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
    
    DaqStrawHit*  sh(nullptr);
    if (fEvent->nshtot > 0) sh = (DaqStrawHit* ) fEvent->sh->UncheckedAt(i);
    // std::cout << std::format("i:{:5d} sd_mnid[i]:{:03d}\n",i,sd->mnid);

    index.plane = sd->plane();
    index.panel = sd->panel();
    index.ch    = sd->straw();
    index.slot  = index.plane/2;
    index.pnl12 = 6*(index.plane%2)+index.panel;
    index.mnid  = sd->mnid;

    int offline_panel = index.slot*12+index.pnl12;
    fHist[0]->h_occup->Fill(offline_panel);
    fHist[0]->h_occup_2d->Fill(index.ch,offline_panel);

    SlotHist_t* slot_h = fHist[0]->slot[index.slot];
    FillSlotHistograms(slot_h,&index,sd,sh);

    PanelHist_t* panel_h = slot_h->panel[index.pnl12];
    FillPanelHistograms(panel_h,&index,sd,sh);

    ChannelHist_t* ch_h = panel_h->ch[index.ch];
    FillChannelHistograms(ch_h,&index,sd,sh);

    // there is one-to-one correspondence between hits and digis

    if (sh and (sh->edep > 0.0005)) {
      fHist[1]->h_occup->Fill(offline_panel);
      fHist[1]->h_occup_2d->Fill(index.ch,offline_panel);

      slot_h = fHist[1]->slot[index.slot];
      FillSlotHistograms(slot_h,&index,sd,sh);
     
      panel_h = slot_h->panel[index.pnl12];
      FillPanelHistograms(panel_h,&index,sd,sh);
     
      ch_h = panel_h->ch[index.ch];
      FillChannelHistograms(ch_h,&index,sd,sh);
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
Long64_t plot_n001_occup::LoadTree(Long64_t entry) {
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
void plot_n001_occup::Loop(int NEvents) {

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
//-----------------------------------------------------------------------------
//  calculate ave
//-----------------------------------------------------------------------------
    // and reinitialize
    for (int i=0; i<36; i++) {
      t05[i] = 0.;
      n05[i] = 0;
    }
    if (fEvent->nshtot > 0) {
      // if straw hits are present at all, the number of straw hits is the same
      // as the number of straw digis
      for (int i=0; i<fEvent->nsdtot; i++) {
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
            fHist[0]->h_dt05[j][i]->Fill(dt05[j][i]);
          }
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

  PrintNoisyChannels();
}

//-----------------------------------------------------------------------------
int plot_n001_occup::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n001_occup::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n001_occup::PrintHistograms(int ISet) {

  gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_set_{:02}_occup.pdf",fRunNumber,ISet);
  
  gStyle->SetStatW(0.30);   // wider (NDC)

  float hmax = (int(fMaxEvent/1.e4)+1)*1.e2;
  
  for (int is=0; is<18; is++) {
    TCanvas c(Form("c_%02i",is),Form("c_%02i",is),1600,1800);
    c.Divide(3,4);
    for (int ip=0; ip<12; ip++) {
      TH1F* h = fHist[ISet]->slot[is]->panel[ip]->h_occup;
      c.cd(ip+1);
     
      // normalization to the rate in Hz:

      float total_time = fMaxEvent*fRunInfo.ew_length;    // assume no prescale
      float rejf       = fRunInfo.cfo_rate/fRunInfo.trigger_rate;
      float eff_time   = total_time/rejf;
      h->Scale(1./eff_time);
      h->SetMaximum(hmax);
      h->SetMinimum(0.5);
      gPad->SetLogy(kTRUE);
      h->Draw("hist");
      // make statbox transparent
      gPad->Update();           // create statbox

      auto st = (TPaveStats*)h->FindObject("stats");
      if (st) {
        st->SetFillStyle(0);    // transparent
        st->SetBorderSize(1);   // optional
      }
      gPad->Modified();
      gPad->Update();
    }

    if (is == 0) {
      c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
    }
    else if (is == 17) {
      c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
    }
    else {
      c.Print(fn.data());     // or .pdf, .root, ...
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// print channels which contain more than Percentage of the total, default - 1%
//-----------------------------------------------------------------------------
int plot_n001_occup::PrintNoisyChannels(float Percentage) {

  int   nbx  = fHist[0]->h_occup_2d->GetNbinsX();
  int   nby  = fHist[0]->h_occup_2d->GetNbinsY();
  float qtot = fHist[0]->h_occup_2d->GetEntries();

  std::cout << std::format("qtot:{}\n",qtot);
  
  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      float nxy = fHist[0]->h_occup_2d->GetBinContent(ix+1,iy+1);
      if (nxy > qtot*Percentage) {
        // ix is the channel number
        int plane = iy / 6;
        int panel = iy % 6;
        TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
        
        std::cout << std::format("ix:{:3} iy:{:3} plane:{:02d} panel:{} channel:{:02d} dtc_id:{:02d} link:{} MN{:03d} nxy:{:8} percent:{:8.4f}\n",
                                 ix,iy,plane,panel,ix,tpmd->dtc_id,tpmd->link,tpmd->mnid,nxy,float(nxy)/qtot);
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n001_occup::PrintDt05Histograms() {

  gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_n002_tc.pdf",fRunNumber);
  
  gStyle->SetStatW(0.30);   // wider (NDC)

  for (int ic=0; ic<3; ic++) {
    TCanvas c(Form("c_%02i",ic),Form("c_%02i",ic),1600,1800);
    c.Divide(4,3);
    for (int ip=0; ip<12; ip++) {
      int plane = ic*12+ip;
      if (plane == 0) continue;
      TH1F* h = fHist[0]->h_dt05[plane][plane-1];
      h->GetXaxis()->SetRangeUser(-100,100);
      c.cd(ip+1);
      // float hmax = (int(fMaxEvent/1.e6)+1)*1e6;
      //      // normalization to the rate :

      //      float input_rate = 1.e4;   // 10 kHz
      //      float scale = input_rate/(fNEvents+1.e-12);
      //      h->Scale(scale);
      //      h->SetMaximum(hmax);
      //      gPad->SetLogy(kTRUE);
      h->Fit("gaus","","",-100,100);
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

  return 0;
}


//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int plot_n001_occup::PlotOccupMap(int ISet, int Print) {

  // gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_set_{:02}_occup_map.pdf",fRunNumber,ISet);
  
  gStyle->SetStatW(0.30);   // wider (NDC)

  TCanvas* c = new TCanvas(Form("c_occup_map"),Form("c_occup_map"),1600,900);

  TH2F* h2 = new TH2F("h_occup_map","occup map: slot vs panel",12,0,12,18,0,18);

  for (int is=0; is<18; is++) {
    for (int ip=0; ip<12; ip++) {
      TH1F* h = fHist[ISet]->slot[is]->panel[ip]->h_occup;
      int nent = h->GetEntries();
      h2->SetBinContent(ip+1,is+1,nent);
    }
  }
  h2->Draw();
                                        // make statbox transparent
  gPad->Update();                       // create statbox

  auto sbox = (TPaveStats*) h2->FindObject("stats");
  if (sbox) {
    sbox->SetFillStyle(0);    // transparent
    sbox->SetBorderSize(1);   // optional
  }
  gPad->Modified();
  gPad->Update();

  if (Print != 0) {
    // if (is == 0) {
    //   c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
    // }
    // else if (is == 17) {
    //   c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
    // }
    // else {
    c->Print(fn.data());     // or .pdf, .root, ...
    //    }
  }

  return 0;
}

