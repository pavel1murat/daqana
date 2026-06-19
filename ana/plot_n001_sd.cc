///////////////////////////////////////////////////////////////////////////////
// run on digi_ntuple , process data of pulse injection runs
// for each straw, fill 4 histograms:
// - tcal : T(CAL),
// - bl   : baseline from the waveform fit
// - ph   : pulse height from the waveform fit
// - dt   : T(CAL)-T(HV)
// also for each panel plot "occupancy" histogram to see which channels
// have been pulsed in a given run
// for a given run, the histogtams are stored in the folder with the run number
// reference channel can be different for different panels and it also depends
// on the run number
// for each of 7-8 runs each panel has its own reference channel and then they will
// linked together
// 
// |------------+------------------+-------------|
// | run number | pulsed channels  | ref channel |
// |------------+------------------+-------------|
// |     121170 | 0, 4,  8, 12     |          12 |
// |     121171 | 1, 4,  9, 12     |          12 |
// |     121172 | 1, 5,  9, 13 ... |           9 |
// |     121173 | 2, 5, 10, 13     |          13 |
// |     121174 | missing          |             |
// |     121175 | 2, 6, 10, 14     |             |
// |     121176 | 3, 6, 11, 14     |             |
// |     121177 | 3, 7, 11, 15 ... |             |
// |------------+------------------+-------------|
// 
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include "ana/plot_n001_sd.hh"



//-----------------------------------------------------------------------------
plot_n001_sd::plot_n001_sd(const char* Name, int RunNumber, const char* Fn = 0) : TNamed(Name,Name), fChain(0) {

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
  RunData_t rd;
  rd.run_number        = 122629;
  //  rd.ref_channel       = 13;

  fRefChannel          = 21;

  fTopFolder = (TFolder*) gROOT->GetRootFolder()->FindObject(Name);
  
  if (fTopFolder == nullptr) {
    fTopFolder = gROOT->GetRootFolder()->AddFolder(Name,Name);
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
  BookHistograms(&fHist,fRunFolder);
}

//-----------------------------------------------------------------------------
plot_n001_sd::~plot_n001_sd() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n001_sd::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  name  = "h_ch";
  title = std::format("run:{} : pulsed channel",fRunNumber);
  fBook->HBook1F(Hist->h_ch,name.data(),title.data(),96,0,96,Folder);

  name  = "h_ph";
  title = std::format("run:{} : pulse height (ALL)",fRunNumber);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),200,0,200,Folder);

  name  = "h_tdc0";
  title = std::format("run:{} : TDC0 (all), us",fRunNumber);
  fBook->HBook1F(Hist->h_tdc0,name.data(),title.data(),1000,0,100,Folder);

  name  = "h_plane";
  title = std::format("run:{} : plane number",fRunNumber);
  fBook->HBook1F(Hist->h_plane,name.data(),title.data(),36,0,36,Folder);

  for (int i=0; i<6; i++) {
    name  = std::format("h_dt20_{}_panel_{}",GetName(),i);
    title = std::format("run:{} : plane dt20 panel:{}",fRunNumber,i);
    fBook->HBook1F(Hist->h_dt20[i],name.data(),title.data(),36,0,36,Folder);
  }

  name  = std::format("h_panel_dt");
  title = std::format("run:{} {} : panel vs dt=t_i-t(21)",fRunNumber,prefix);
  fBook->HBook2F(Hist->h_panel_dt,name.data(),title.data(),500,-25,25,216,0,216,Folder);
 
  name  = std::format("h_panel_dt_111");
  title = std::format("run:{} {} : panel vs dt_111",fRunNumber,prefix);
  fBook->HBook2F(Hist->h_panel_dt_111,name.data(),title.data(),500,-25,25,216,0,216,Folder);

  for (int ip=0; ip<216; ip++) {
    name  = std::format("h_panel_dt_111_vs_evn_{:03d}",ip);
    title = std::format("run:{} panel {:03d} dt_111 vs evn",fRunNumber,ip);
    fBook->HProf(Hist->h_panel_dt_111_vs_evn[ip],name.data(),title.data(),10000,0,1e6,-50,50,Folder);
  }
 
  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n001_sd::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n001_sd::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
void plot_n001_sd::fill_histograms() {
  // filling histograms: plot time differences between

  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;
  
  for (int i=0; i<fEvent->nsdtot; i++) {
    DaqStrawDigi* sdi = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
    // std::cout << std::format("i:{:5d} sd_mnid[i]:{:03d}\n",i,sd->mnid);

    fHist.h_ph->Fill(sdi->ph);

    if (sdi->ph < 60) continue;
    
    fHist.h_ch->Fill(sdi->straw());
    fHist.h_plane->Fill(sdi->plane());
    
    int pln = sdi->plane();
    int pnl = sdi->panel();
    int ich = sdi->straw();

    float ti = sdi->tdc0*5./256*1.e-3;
    fHist.h_tdc0->Fill(ti);

    int ip = 6*pln+pnl;                 // index of the panel
    DaqStrawDigi* sdr = fSdr[ip];
    
    if ((sdr != nullptr) and (sdi != sdr)) {
      float tr = sdr->tdc0*(5./256.)*1.e-3;
      float dt = ti-tr;
      fHist.h_panel_dt->Fill(dt,ip);
    }

    if (sdi == sdr) {
                                        // sdi: first hit in the panel in the reference channel
                                        // compare it to the one in panel 111
      DaqStrawDigi* sdr_111 = fSdr[111];
      if ((sdr_111 != nullptr) and (sdi != sdr_111)) {
        float tr_111 = sdr_111->tdc0*(5./256.)*1.e-3;
        float dt_111 = ti-tr_111;
        fHist.h_panel_dt_111->Fill(dt_111,ip);

        // std::cout << std::format("evn:{:6d} tr_111:{:10.5f} tr_112:{:10.5f}\n",fEvent->evn,tr_111,ti);

        float dt = tr_111-ti;
        if (dt < 0) dt += 20;
        fHist.h_panel_dt_111_vs_evn[ip]->Fill(fEvent->evn,dt);
      }
    }

  }
}


//-----------------------------------------------------------------------------
Long64_t plot_n001_sd::LoadTree(Long64_t entry) {
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
void plot_n001_sd::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;

  fMaxEvent = -1;
  
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (fEvent->evn > fMaxEvent) {
      fMaxEvent = fEvent->evn;
    }

    // std::cout << std::format("------- run:{:6d} srn:{:6d} evn:{:8d}\n",
    //                          fEvent->run,fEvent->srn,fEvent->evn);

                                        // fSdr[ip] first ch=21 digi in a panel=ip in this event
    for (int ip=0; ip<216; ip++) {
      fSdr[ip] = nullptr;
    }

    for (int i=0; i<fEvent->nsdtot; i++) {
      DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
      // std::cout << std::format("i:{:5d} sd_mnid[i]:{:03d}\n",i,sd->mnid);
      
      // TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_mnid(sd->mnid);
      // if (tpmd == nullptr) {
      //   std::cout << std::format("event:{:08d} ERROR: for digi #{} : wrong MNID:{}, SKIP\n",
      //                            fEvent->evn,i,sd->mnid);
      //   continue;
      // }

      int pln = sd->plane();
      int pnl = sd->panel();
      int ich = sd->straw();

      int ip  = 6*pln+pnl;

      if ((ich == fRefChannel) and (sd->ph > 60)) { 
        // first good hit, channel=21 for each plane/panel
        // std::cout << std::format("event:{} plane:{:2d} panel:{} index:{}\n",
        //                          fEvent->evn,pln,pnl,i);
        if (fSdr[ip] == nullptr) fSdr[ip] = sd;
      }
    }
//-----------------------------------------------------------------------------
// prep done, now fill histograms
//-----------------------------------------------------------------------------
    fill_histograms();
  }
//-----------------------------------------------------------------------------
// post-loop
//-----------------------------------------------------------------------------
  for (int ip=0; ip<216; ip++) {
    fHist.h_panel_dt_111_vs_evn[ip]->GetXaxis()->SetRangeUser(0,fMaxEvent+100);
  }
  
}

//-----------------------------------------------------------------------------
int plot_n001_sd::ResetHistograms() {
  fHist.h_ch->Reset();
  fHist.h_ph->Reset();
  fHist.h_plane->Reset();
  fHist.h_panel_dt->Reset();
  fHist.h_panel_dt_111->Reset();
  for (int ip=0; ip<216; ip++) {
    fHist.h_panel_dt_111_vs_evn[ip]->Reset();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n001_sd::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}
