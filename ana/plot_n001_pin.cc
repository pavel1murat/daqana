///////////////////////////////////////////////////////////////////////////////
// run on digi_ntuple , process data of pulse injection runs
// plot histograms:
// - sd.tdc0                  : for all (just to see the time distribution)
// - sd.straw()               : to see pulsed channels, assumed to be the same for all panels
// - sd.ph                    : pulse height, to see where to cut
// - plane                    : to see which stations were included
// - h_panel_dt               : 2D, panel number vs deltaT, deltaT = T_i-T(first hit in a reference channel)
//                              reference channel is run-dependent, here comes in a run DB
//                              so far, ran only on data of a single run (122629)
// - h_panel_dt_111           : 2D, panel number vs deltaT, deltaT = T_i-T(first hit in a reference channel of panel 111)
// - h_panel_dt_111_vs_evn[i] : deltaT (above) vs event number

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
// |     122655 | 5, 13, 21...     |          12 |
// |     122656 | 6, 14, 22        |          12 |
// |     122657 | 7, 15, 23        |           9 |
// |     122658 | 0,  8, 16,       |          13 |
// |     122659 | 1,  9,           |             |
// |     122660 | 2, 10, 18, ...   |             |
// |     122661 | 3, 11, 19, ...   |             |
// |     122662 | 4, 12, 20, ...   |             |
// |------------+------------------+-------------|
//
// - make ntuples
// - run daqana/scripts/make_pin_dt_hist.C to produce a histogram file
// - run ... on a histogram file to produce calibrations
//   - if fit didn't converge (or a channel has no data) set deltaT to 0
//   - need to figure out the sign of correction
/////////////////////////////////////////////////////////////////////////////////
/*
  .L v001/daqana/scripts/plot_n002_hist_001.C
  //
  x->SaveHist("pulse_injection_120807_120808.hist");
*/
#include "ana/plot_n001_pin.hh"

//-----------------------------------------------------------------------------
plot_n001_pin::plot_n001_pin(const char* Name, int RunNumber, const char* Fn) : TNamed(Name,Name), fChain(0) {
  std::string dsid("vst00s000r000n001");
  std::string dir ("/data/mu2e/mu2etrk/datasets");
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/{}/nts.mu2e.trk.{}.{:06d}_000001.root",dir,dsid,dsid,RunNumber);
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

  RunInfoDb::Instance()->GetRunInfo(RunNumber,&fRunInfo);

  fRefChannel      = fRunInfo.ref_channel;    // ok for 122655

  fBook      = new Booking(fRunFolder);

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  printf("tree: %p\n",(void*) tree);

  Init(tree);
  BookHistograms(&fHist,fRunFolder);
}

//-----------------------------------------------------------------------------
plot_n001_pin::~plot_n001_pin() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//-----------------------------------------------------------------------------
int plot_n001_pin::BookChannelHistograms(ChannelHist_t* Hist, TrkPanelMap_t::Data_t* Tpmd, int Ich, TFolder* Folder) {
  int rc(0);

  std::string prefix = std::format("run:{} MN{:03d} plane:{:02d} panel:{} ich:{:02d}",fRunNumber,Tpmd->mnid,Tpmd->plane,Tpmd->panel,Ich);
  std::string name, title;

  name  = "tdc0";
  title = std::format("{}: TDC0",prefix);
  fBook->HBook1F(Hist->h_tdc0,name.data(),title.data(),96,0,96,Folder);

  name  = "bl";
  title = std::format("{}: baseline",prefix);
  fBook->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,500,Folder);

  name  = "ph";
  title = std::format("{}: pulse height",prefix);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),200,0,200,Folder);

  name  = "dt01";
  title = std::format("{}: dt01 =T(cal)-T(hv)",prefix);
  fBook->HBook1F(Hist->h_dt01,name.data(),title.data(),400,-20,20,Folder);
  
  return rc;
}


//-----------------------------------------------------------------------------
// Ip : offline index: 6*SidPlane + SidPanel
// 'Sid'='used in offline sid'
//-----------------------------------------------------------------------------
int plot_n001_pin::BookPanelHistograms(PanelHist_t* Hist, TrkPanelMap_t::Data_t* Tpmd, TFolder* Folder) {
  int rc(0);

  std::string prefix = std::format("");
  std::string name, title;

  name  = "occup";
  title = std::format("run:{} : MN{:03d} plane:{:02d} panel:{}",fRunNumber,Tpmd->mnid,Tpmd->plane,Tpmd->panel);
  fBook->HBook1F(Hist->h_ch,name.data(),title.data(),96,0,96,Folder);

  for (int i=0; i<96; i++) {
    std::string folder_name = std::format("ch_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());

    BookChannelHistograms(&Hist->channel[i],Tpmd,i,fol);
  }
  
  return rc;
}

//-----------------------------------------------------------------------------
int plot_n001_pin::BookHistograms(Hist_t* Hist, TFolder* Folder) {
  int rc(0);
                                        // use offline panel/plane numbering
  for (int i=0; i<216; i++) {
    int plane = i / 6;
    int panel = i % 6;
                                        // use only powered on panels
    if (fRunInfo.plane_flag[plane] == 0) continue;

    TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
        
    std::string folder_name = std::format("pnl_{:03d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookPanelHistograms(&Hist->panel[i],tpmd,fol);
  }

  return rc;
}

//-----------------------------------------------------------------------------
Int_t plot_n001_pin::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n001_pin::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);
}

//-----------------------------------------------------------------------------
int plot_n001_pin::FillChannelHistograms(ChannelHist_t* Hist, DaqStrawDigi* Sd) {
  int rc(0);

  Hist->h_tdc0->Fill(Sd->tdc0);
  Hist->h_bl->Fill(Sd->bl);
  Hist->h_ph->Fill(Sd->ph);
  float dt01_ns = (Sd->tdc0-Sd->tdc1)*5./256.;
  Hist->h_dt01->Fill(dt01_ns);
  
  return rc;
}

//-----------------------------------------------------------------------------
int plot_n001_pin::FillPanelHistograms(PanelHist_t* Hist, DaqStrawDigi* Sd) {
  int rc(0);

  int ich = Sd->straw();
  Hist->h_ch->Fill(ich);

  FillChannelHistograms(&Hist->channel[ich],Sd);
  
  return rc;
}

//-----------------------------------------------------------------------------
int plot_n001_pin::FillHistograms() {
  int rc(0);
  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;
  
  for (int i=0; i<fEvent->nsdtot; i++) {
    DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);

    int ch    = sd->straw();
    int plane = sd->plane();
    int panel = sd->panel();

    int ip    = plane*6+panel;

    FillPanelHistograms(&fHist.panel[ip],sd);
  }

  return rc;
}


//-----------------------------------------------------------------------------
Long64_t plot_n001_pin::LoadTree(Long64_t entry) {
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
void plot_n001_pin::Loop(int NEvents) {

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

    //                                     // fSdr[ip] first ch=21 digi in a panel=ip in this event
    // for (int ip=0; ip<216; ip++) {
    //   fSdr[ip] = nullptr;
    // }

//     for (int i=0; i<fEvent->nsdtot; i++) {
//       DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
//       // std::cout << std::format("i:{:5d} sd_mnid[i]:{:03d}\n",i,sd->mnid);
//       
//       // TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_mnid(sd->mnid);
//       // if (tpmd == nullptr) {
//       //   std::cout << std::format("event:{:08d} ERROR: for digi #{} : wrong MNID:{}, SKIP\n",
//       //                            fEvent->evn,i,sd->mnid);
//       //   continue;
//       // }
// 
//       int pln = sd->plane();
//       int pnl = sd->panel();
//       int ich = sd->straw();
// 
//       int ip  = 6*pln+pnl;
// 
//       if ((ich == fRefChannel) and (sd->ph > 60)) { 
//         // first good hit, channel=21 for each plane/panel
//         // std::cout << std::format("event:{} plane:{:2d} panel:{} index:{}\n",
//         //                          fEvent->evn,pln,pnl,i);
//         if (fSdr[ip] == nullptr) fSdr[ip] = sd;
//       }
//     }
//-----------------------------------------------------------------------------
// prep done, now fill histograms
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
int plot_n001_pin::ResetHistograms() {
  // fHist.h_ch->Reset();
  // fHist.h_ph->Reset();
  // fHist.h_plane->Reset();
  // fHist.h_panel_dt->Reset();
  // fHist.h_panel_dt_111->Reset();
  // for (int ip=0; ip<216; ip++) {
  //   fHist.h_panel_dt_111_vs_evn[ip]->Reset();
  // }
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n001_pin::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

