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
#include "ana/plot_n001_noise.hh"


//-----------------------------------------------------------------------------
plot_n001_noise::plot_n001_noise(const char* Name, int RunNumber, const char* Fn) : TNamed(Name,Name), fChain(0) {

  std::string dir("/data/mu2e/mu2etrk/datasets/vst00s000r000n001"); 
  
  TFile *f(nullptr);
  if (Fn != nullptr) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(Fn);
    if (!f || !f->IsOpen()) {
      f = new TFile(Fn);
    }
  }
  else {
    std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000n001.{:06d}_000001.root",dir,RunNumber);
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
plot_n001_noise::~plot_n001_noise() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n001_noise::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  // name  = "h_ch";
  // title = std::format("run:{} : pulsed channel",fRunNumber);
  // fBook->HBook1F(Hist->h_ch,name.data(),title.data(),96,0,96,Folder);

  // name  = "h_ph";
  // title = std::format("run:{} : pulse height (ALL)",fRunNumber);
  // fBook->HBook1F(Hist->h_ph,name.data(),title.data(),200,0,200,Folder);

  // name  = "h_bl";
  // title = std::format("run:{} : baseline (ALL)",fRunNumber);
  // fBook->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,500,Folder);

  // name  = "h_tdc0";
  // title = std::format("run:{} : TDC0 (all), us",fRunNumber);
  // fBook->HBook1F(Hist->h_tdc0,name.data(),title.data(),1000,0,100,Folder);

  // name  = "h_plane";
  // title = std::format("run:{} : plane number",fRunNumber);
  // fBook->HBook1F(Hist->h_plane,name.data(),title.data(),36,0,36,Folder);

  // for (int i=0; i<6; i++) {
  //   name  = std::format("h_dt20_{}_panel_{}",GetName(),i);
  //   title = std::format("run:{} : plane dt20 panel:{}",fRunNumber,i);
  //   fBook->HBook1F(Hist->h_dt20[i],name.data(),title.data(),36,0,36,Folder);
  // }

  // name  = std::format("h_panel_dt");
  // title = std::format("run:{} {} : panel vs dt=t_i-t(21)",fRunNumber,prefix);
  // fBook->HBook2F(Hist->h_panel_dt,name.data(),title.data(),500,-25,25,216,0,216,Folder);
 
  // name  = std::format("h_panel_dt_111");
  // title = std::format("run:{} {} : panel vs dt_111",fRunNumber,prefix);
  // fBook->HBook2F(Hist->h_panel_dt_111,name.data(),title.data(),500,-25,25,216,0,216,Folder);

  // for (int ip=0; ip<216; ip++) {
  //   name  = std::format("h_panel_dt_111_vs_evn_{:03d}",ip);
  //   title = std::format("run:{} panel {:03d} dt_111 vs evn",fRunNumber,ip);
  //   fBook->HProf(Hist->h_panel_dt_111_vs_evn[ip],name.data(),title.data(),10000,0,1e6,-50,50,Folder);
  // }
 
  return 0;
}

//-----------------------------------------------------------------------------
Int_t plot_n001_noise::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n001_noise::Init(TTree *tree) {
  
  // #include "daqana/scripts/daqana_nt_init.C"

  fChain = tree;
  fCurrent = -1;
  // fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt",&fEvent);
  // fChain->SetBranchAddress("evt.sd",&fSd);

  for (int i=0; i<500; i++) {
    fDat[i] = nullptr;
  }
}

//-----------------------------------------------------------------------------
void plot_n001_noise::fill_histograms() {
  // filling histograms: plot time differences between

  // DaqStrawDigi* sdr = (DaqStrawDigi*) fEvent->sd->UncheckedAt(fHitIndex[fRefPlane][0]);
  // float tr = sdr->tdc0*(5./256.)*1.e-3;

  

}


//-----------------------------------------------------------------------------
Long64_t plot_n001_noise::LoadTree(Long64_t entry) {
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
void plot_n001_noise::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << std::format("nentries:{}\n",nentries);

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;

  //
  std::vector<noise_t>* pdat;
  noise_t* dr;                          // data record

  int first_event(-1);
  int index;
  
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    if (first_event == -1) {
      first_event = (fEvent->evn/1000)*1000;
      index = 0;
    }
    
    for (int i=0; i<fEvent->nsdtot; i++) {
      DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
   
      if (fDat[sd->mnid] == nullptr) {
        fDat[sd->mnid] = new std::vector<noise_t>();
        pdat          = fDat[sd->mnid];
        fListOfDat.push_back(pdat);
        
        pdat->push_back(noise_t(first_event));
      }
      else {
        pdat          = fDat[sd->mnid];
      }
                                        // last element
      dr = &pdat->back();
      // one data point - per 1000 events
      /*
      current_index = (fEvent->evn-first_event)/1000;
      
      if (current_index != index) {
                                        // and add a new record
        pdat->push_back(noise_t(current_index));
        dr = &pdat->back();
      }
                                        // accumulating
      dr->evmax    = fEvent->evn;
      dr->nevents += 1;
      dr->sumx    += sd->bl;
      */
      // pnos->sumx2   += sd->bl*sd->bl;
    }
  }
//-----------------------------------------------------------------------------
// second loop
//-----------------------------------------------------------------------------
/*
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    if (first_event == -1) {
      first_event = (fEvent->evn/1000)*1000;
      index = 0;
    }
    
    for (int i=0; i<fEvent->nsdtot; i++) {
      DaqStrawDigi* sd = (DaqStrawDigi*) fEvent->sd->UncheckedAt(i);
   
      if (fDat[sd->mnid] == nullptr) {
        fDat[sd->mnid] = new std::vector<noise_t>();
        pdat          = fDat[sd->mnid];
        fListOfDat.push_back(pdat);
        
        pdat->push_back(noise_t(first_event));
      }
      else {
        pdat          = fDat[sd->mnid];
      }
                                        // last element
      dr = &pdat->back();
      // one data point - per 1000 events
      if (fEvent->evn-dr->evmin >= 1000) {
                                        // and add a new record
        first_event = (fEvent->evn/1000)*1000;
        pdat->push_back(noise_t(first_event));
        dr = &pdat->back();
      }
                                        // accumulating
      dr->evmax    = fEvent->evn;
      dr->nevents += 1;
      dr->sumx    += sd->bl;
      // pnos->sumx2   += sd->bl*sd->bl;
    }
  }
*/
//-----------------------------------------------------------------------------
// post-loop
//-----------------------------------------------------------------------------
  int npanels = fListOfDat.size();
  for (int i=0; i<npanels; i++) {
    std::vector<noise_t>* pdat = fListOfDat[i];
    int nbins = pdat->size();

    fHist.h_noise[i] = new TH1D(Form("h_%02i",i),Form("h_%02i",i),nbins,0,nbins);

    for (int ib=0; ib<nbins; ib++) {
      noise_t* dr = &pdat->at(ib);
      double sigm = dr->sigm();
      fHist.h_noise[i]->SetBinContent(ib+1,sigm);
      fHist.h_noise[i]->SetBinError(ib+1,0);
    }
  }
  
}

//-----------------------------------------------------------------------------
int plot_n001_noise::ResetHistograms() {
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n001_noise::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}
