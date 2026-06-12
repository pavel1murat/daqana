#ifndef plot_n002_hist_h
#define plot_n002_hist_h
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
#include <format>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.
// #include "daqana/obj/DaqEvent.hh"
// #include "TObject.h"
// #include "daqana/obj/DaqStrawDigi.hh"
// #include "daqana/obj/DaqStrawHit.hh"
// #include "daqana/obj/DaqComboHit.hh"
// #include "daqana/obj/DaqTimeCluster.hh"
// #include "daqana/obj/DaqTrack.hh"
// #include "daqana/obj/DaqTrkStrawHit.hh"
// #include "daqana/obj/DaqSegment.hh"
#include "daqana/obj/TrkPanelMap_t.hh"

#include "booking.C"

class plot_n002_hist {
public :

  enum {
    kNStations         = 18,
    kNPanelsPerStation = 12,
  };
//-----------------------------------------------------------------------------
// fit results - data structures
//-----------------------------------------------------------------------------
  struct ChannelFitResult_t {
    double p[3];                        // gaussian fit
    double e[3];
    double chi2dof;
  };

  struct FitResult_t {
    int                 ref_channel;
    ChannelFitResult_t  cfr[96];
  };

  FitResult_t  fFitResult[kNStations][kNPanelsPerStation];

  static FitResult_t  fgFitResult[kNStations][kNPanelsPerStation];
//-----------------------------------------------------------------------------
// data structures
//-----------------------------------------------------------------------------
  struct PanelData_t {
    int    ref_channel;
    int    ref_tdc;                     // in units of 5./256 ns , approx 20 ps
  };

  struct TrackerData_t {
    PanelData_t   panel[kNStations][kNPanelsPerStation];
  } fTrackerData;
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct StrDigiHist_t {
    TH1F* h_tcal;
    TH1F* h_bl;
    TH1F* h_ph;
    TH1F* h_dt01;                       // T(CAL,i)-T(HV,i)
    TH1F* h_dtref;                      // T(CAL,i)-T(CAL,ref_channel)
  };
  
  struct PanelHist_t {
    StrDigiHist_t  ch[96];
    TH1F* h_occup;
  };
  
  struct StationHist_t {
    PanelHist_t  panel[kNPanelsPerStation];
  };
  
  struct Hist_t {
    StationHist_t  station[kNStations];
  } fHist;
//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
  TFolder*       fTopFolder;
  TFolder*       fRunFolder;
  
  Booking*       fBook;

  TrkPanelMap_t* fTpm;
  
  int            fRunNumber;
  int            fStation1;
  int            fStation2;
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

                                        // will the same directory work ?
#include "daqana_nt_format.hh"
                                        // for independent runs, the name should eb the same..
                                        // make it different to process the same run with different refence channels
  
  plot_n002_hist(const char* Name, int RunNumber, int RefChannel, int Station1, int Station2);
  
  plot_n002_hist(const char* Fn = "nts.mu2e.trk.vst00s000r000n001.120718_000001.root");

  virtual ~plot_n002_hist();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void             Loop          (int NEvents = -1);

  
  int              BookStrDigiHistograms(StrDigiHist_t* Hist, int Station, int Panel, int Channel, TFolder* Folder);
  int              BookPanelHistograms  (PanelHist_t*   Hist, int Station, int Panel, TFolder* Folder);
  int              BookStationHistograms(StationHist_t* Hist, int Station, TFolder* Folder);
  int              BookHistograms       (Hist_t* Hist, TFolder* Folder);

  int              CloneFitResults      ();

  int              FillStrDigiHistograms(StrDigiHist_t* Hist, int Ind, PanelData_t* Pd);
  int              FillPanelHistograms  (PanelHist_t*   Hist, int Ind);

  int              GetRefChannel(int Station, int Panel);

  int              PrintFitResults(int Station = -1, int Panel = -1);
  int              FitHistograms(int Station1, int Station2, int PrintLevel = 0);
  int              ResetHistograms();
  int              SaveHistograms(const char* Filename);
  int              AddFitResults (plot_n002_hist* Ph);
};

plot_n002_hist::FitResult_t  plot_n002_hist::fgFitResult[kNStations][kNPanelsPerStation];

//-----------------------------------------------------------------------------
plot_n002_hist::plot_n002_hist(const char* Name, int RunNumber, int RefChannel, int Station1, int Station2) : fChain(0) {

  std::string dir("/data/mu2e/mu2etrk/datasets/vst00s000r000n001"); 
  std::string fn = std::format("{}/nts.mu2e.trk.vst00s000r000n001.{:06d}_000001.root",dir,RunNumber);
  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
  }

  // plane 20: MN 122, 080, 133, 150, 199, 278
  // plane 21: MN 112, 072, 202, 193, 200, 270
  // plane 22: MN 053, 052, 169, 035, 039, 043
  // plane 23: MN 151, 254, 225, 229, 275, 274
  // plane 24: MN 132, 145, 082, 207, 139, 155
  // plane 25: MN 055, 069, 064, 062, 063, 067

  fRunNumber = RunNumber;
  fStation1  = Station1;
  fStation2  = Station2;

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

  Init(tree);
  BookHistograms(&fHist,fRunFolder);

  for (int is=Station1; is<=Station2; is++) {
    for (int ip=0; ip<kNPanelsPerStation; ip++) {
      PanelData_t* pdt = &fTrackerData.panel[is][ip];
//-----------------------------------------------------------------------------
// try to process most of the channels, deal with exceptions later
// how do I identify exceptions ?
// for example run 121170 panel 1 (MN199) has channel 12 masked off
//-----------------------------------------------------------------------------
      pdt->ref_channel = RefChannel; // ref_channels[is][ip];
    }
  }
//-----------------------------------------------------------------------------
// init fit result structure
//-----------------------------------------------------------------------------
  for (int is=Station1; is<=Station2; is++) {
    for (int ip=0; ip<kNPanelsPerStation; ip++) {
      FitResult_t*  fr = &fFitResult   [is][ip];
      PanelData_t* pdt = &fTrackerData.panel[is][ip];
      fr->ref_channel      = pdt->ref_channel;
      for (int ich=0; ich<96; ich++) {
        ChannelFitResult_t* cfr = &fr->cfr[ich];
        cfr->chi2dof = -1;
        for (int i=0; i<3; i++) {
          cfr->p[i] = -999.0;           // to be seen easily when printed
          cfr->e[i] = -999.0;
        }
      }
    }
  }
}


//-----------------------------------------------------------------------------
plot_n002_hist::plot_n002_hist(const char* Fn) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  std::string fn = std::format("/data/mu2e/mu2etrk/datasets/vst00s000r000n001/{}",Fn);
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(fn.data());
  }

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  Init(tree);

  BookHistograms(&fHist,fRunFolder);
}


//-----------------------------------------------------------------------------
plot_n002_hist::~plot_n002_hist() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
int plot_n002_hist::BookStrDigiHistograms(StrDigiHist_t* Hist, int Station, int Panel, int Channel, TFolder* Folder) {
  
  std::string prefix = std::format("s{:02d}_p{:02d}_{:02d}",Station,Panel,Channel);
  std::string name, title;
  
  name  = std::format("h_tcal_{}",prefix);
  title = std::format("run:{} {} : tcal",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_tcal,name.data(),title.data(),1000,0,100000,Folder);
  
  name  = std::format("h_bl_{}",prefix);
  title = std::format("run:{} {} : bl",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,1000,Folder);
  
  name  = std::format("h_ph_{}",prefix);
  title = std::format("run:{} {} : ph",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_ph,name.data(),title.data(),500,0,1000,Folder);
  
  name  = std::format("h_dt01_{}",prefix);
  title = std::format("run:{} {} : t_cal[i]-t_hv[i]",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_dt01,name.data(),title.data(),500,-25,25,Folder);
  
  name  = std::format("h_dtref_{}",prefix);
  title = std::format("run:{} {} : t_cal[i]-t_cal[ref]",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_dtref,name.data(),title.data(),500,-25,25,Folder);
  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookPanelHistograms(PanelHist_t* Hist, int Station, int Panel, TFolder* Folder) {

  std::string prefix = std::format("s{:02d}_p{:02d}",Station,Panel);
  std::string name, title;
  
  name  = std::format("h_{}_occup",prefix);
  title = std::format("run:{} {} : occupancy",fRunNumber,prefix);
  fBook->HBook1F(Hist->h_occup,name.data(),title.data(),96,0,96,Folder);
  
  for (int i=0; i<96; i++) {
    std::string folder_name = std::format("ch_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookStrDigiHistograms(&Hist->ch[i],Station,Panel,i,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookStationHistograms(StationHist_t* Hist, int Station, TFolder* Folder) {

  for (int pnl=0; pnl<kNPanelsPerStation; pnl++) {
    std::string folder_name = std::format("pnl_{:02d}",pnl);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookPanelHistograms(&Hist->panel[pnl],Station,pnl,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  for (int i=fStation1; i<=fStation2; i++) {
    std::string folder_name = std::format("stn_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    BookStationHistograms(&Hist->station[i],i,fol);
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// 'Ind' : straw digi index in an array
// Pd has the reference time
//-----------------------------------------------------------------------------
int plot_n002_hist::FillStrDigiHistograms(StrDigiHist_t* Hist, int Ind, PanelData_t* Pd) {
  Hist->h_ph->Fill(sd_ph[Ind]);
  Hist->h_bl->Fill(sd_bl[Ind]);
  Hist->h_tcal->Fill(sd_tdc0[Ind]);

  float dt01  = (sd_tdc0[Ind]-sd_tdc1[Ind])*5./256.; // follow convention of mu2e::StrawHit : dt = CAL-HV
  float dtref = (sd_tdc0[Ind]-Pd->ref_tdc)*5./256.; // follow convention of mu2e::StrawHit : dt = CAL-HV

  Hist->h_dt01->Fill(dt01);
  Hist->h_dtref->Fill(dtref);

  return 0;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::FillPanelHistograms(PanelHist_t* Hist, int Ind) {

  int straw = sd_sid[Ind] & 0x7f;
  Hist->h_occup->Fill(straw);

  return 0;
}


//-----------------------------------------------------------------------------
Int_t plot_n002_hist::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
void plot_n002_hist::Init(TTree *tree) {
#include "daqana/scripts/daqana_nt_init.C"
}

//-----------------------------------------------------------------------------
Long64_t plot_n002_hist::LoadTree(Long64_t entry) {
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
void plot_n002_hist::Loop(int NEvents) {

  ResetHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  int nev = NEvents;
  if (NEvents <= 0) nev = nentries;
   
  for (int jentry=0; jentry<nev; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // std::cout << std::format("----------- event:{:8} nsdtot:{:6}\n",evn,nsdtot);
//-----------------------------------------------------------------------------
// event initialization
//-----------------------------------------------------------------------------
    for (int is=fStation1; is<=fStation2; is++) {
      for (int ip=0; ip<kNPanelsPerStation; ip++) {
        PanelData_t* pdt = &fTrackerData.panel[is][ip];
        pdt->ref_tdc = -1;
      }
    }
//-----------------------------------------------------------------------------
// for each panel determine the time in its reference channel
//-----------------------------------------------------------------------------
    for (int i=0; i<nsdtot; i++) {
      int straw =  sd_sid[i]        & 0x7f;
      int panel = (sd_sid[i] >>  7) & 0x7;
      int plane = (sd_sid[i] >> 10) & 0x3f;
      int slot  = plane / 2;

      int ip    = (plane % 2)*6 +panel;
//-----------------------------------------------------------------------------
// reference time for panel in a given event: first time in a reference channel
//-----------------------------------------------------------------------------
      PanelData_t* pdt = &fTrackerData.panel[slot][ip];
      if ((straw == pdt->ref_channel) and (pdt->ref_tdc < 0)) {
        pdt->ref_tdc = sd_tdc0[i];
      }
    }
    
    for (int i=0; i<nsdtot; i++) {
      TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_mnid(sd_mnid[i]);
      if (tpmd == nullptr) {
        std::cout << std::format("event:{:08d} ERROR: for digi #{} : wrong MNID:{}, SKIP\n",evn,i,sd_mnid[i]);
        continue;
      }
      int ipln = tpmd->plane;
      int ist  = ipln / 2;
      int ipnl = tpmd->panel + 6*(ipln % 2); // "offline" index, runs from 0 to 11

      if ((ist < fStation1) or (ist > fStation2)) continue;
      
      int ich  = sd_sid[i] & 0x7f;
        
      PanelHist_t*   panel_hist = &fHist.station[ist].panel[ipnl];
      StrDigiHist_t* sd_hist    = &fHist.station[ist].panel[ipnl].ch[ich];
        
      // std::cout << std::format("----------- ist:{} ipnl:{} ich:{}\n",ist,ipnl,ich);
      if (sd_ph[i] > 50.) {
        FillPanelHistograms  (panel_hist,i);
        auto pdt = &fTrackerData.panel[ist][ipnl];
        FillStrDigiHistograms(sd_hist   ,i, pdt);
        // if (ipnl == 11) {
        //                                 // something is wrong with the panel 11, print few first events
        //   if (evn < 10) {
        //     std::cout << std::format("evn:{} mnid:{} ref_channel:{} ref_time0:{:10.3f} time0:{}\n",
        //                              evn,sd_mnid[i],pdt->ref_channel,pdt->ref_tdc*5./256,sd_tdc0[i]*5./256);
        //   }
        // }
      }
    }
  }
}

//-----------------------------------------------------------------------------
int plot_n002_hist::FitHistograms(int Station1, int Station2, int PrintLevel) {
  int rc(0);
  
  for (int is=Station1; is<=Station2; is++) {
    for (int ip=0; ip<kNPanelsPerStation; ip++) {
      for (int ich=0; ich<96; ich++) {
        TH1F* h = fHist.station[is].panel[ip].ch[ich].h_dtref;
        int nx = h->GetNbinsX();
        int nentries = h->GetEntries();
        // skip empty and low stat channels
        if (nentries < 100) continue;
//-----------------------------------------------------------------------------
// apply some straw-man intelligence when fitting - look for the highest bin
// and define the fit window aroung it
//-----------------------------------------------------------------------------
        double nmax = -1.;
        float  tmax = -1.e6;
        for (int ix=0; ix<nx; ix++) {
          double y = h->GetBinContent(ix+1);
          // std::cout << "ix:" << ix << " y:" << y << std::endl;
          if (y > nmax) {
            nmax = y;
            tmax = h->GetBinCenter(ix+1);
          }
        }

        // std::cout << "tmax:" << tmax << " nmax:" << nmax << std::endl;
//-----------------------------------------------------------------------------
// perform fit
//-----------------------------------------------------------------------------
        if (nmax > 0) {
          TFitResultPtr tfr = h->Fit("gaus","sq","",tmax-15,tmax+15);
          
          ChannelFitResult_t* cfr = &fFitResult[is][ip].cfr[ich];
          if ((! tfr->IsValid()) or tfr->IsEmpty()) {
            // printf("# FIT ERROR: channel: %2i\n",ich);
            cfr->chi2dof = -1.;
            
            // for (int ip=0; ip<3; ip++) {
            //   fr->p[ip] = -999.;
            //   fr->e[ip] = -999.;
            // }
            // continue;
          }
          else {
//-----------------------------------------------------------------------------
// renormalize errors returned by the fit to make effective chi2/dof = 1
//-----------------------------------------------------------------------------
            cfr->chi2dof = tfr->Chi2()/tfr->Ndf();
            double sf    = sqrt(cfr->chi2dof);
            
            for (int ip=0; ip<3; ip++) {
              cfr->p[ip] = tfr->Parameter(ip);
              cfr->e[ip] = tfr->Error(ip)*sf;
            }
          }
        }
      }
    }
  }

  if (PrintLevel != 0) PrintFitResults();
  return rc;
}

//-----------------------------------------------------------------------------
// print fit results 
//-----------------------------------------------------------------------------
int plot_n002_hist::PrintFitResults(int Station, int Panel) {
  int rc(0);

  int is1(Station), is2(Station+1);
  if (Station == -1) { is1=fStation1; is2 = fStation2; }
  
  int ip1(Panel), ip2(Panel+1);
  if (Panel == -1) { ip1=0; ip2 = kNPanelsPerStation; }
  
  for (int is=is1; is<is2; is++) {
    for (int ip=ip1; ip<ip2; ip++) {
      
      int plane = 2*is+ip/6;
      int panel = ip%6;
      
      TrkPanelMap_t::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
      int mnid = tpmd->mnid;
      
      for (int ich=0; ich<96; ich++) {
        ChannelFitResult_t* cfr = &fFitResult[is][ip].cfr[ich];
        std::cout << std::format("MN{:03d} {:2d} {:2d} {:2d}",mnid,plane,panel,ich);
        printf(" %11.4f %11.4f %11.4f %11.4f %11.4f\n",
               cfr->p[1],cfr->e[1], cfr->p[2], cfr->e[2], cfr->chi2dof);
      }
    }
  }
  
  return rc;
}

//-----------------------------------------------------------------------------
int plot_n002_hist::ResetHistograms() {
  for (int is=fStation1; is<=fStation2; is++) {
    for (int ip=0; ip<kNPanelsPerStation; ip++) {
      fHist.station[is].panel[ip].h_occup->Reset();
      for (int ich=0; ich<96; ich++) {
        fHist.station[is].panel[ip].ch[ich].h_tcal->Reset();
        fHist.station[is].panel[ip].ch[ich].h_bl->Reset();
        fHist.station[is].panel[ip].ch[ich].h_ph->Reset();
        fHist.station[is].panel[ip].ch[ich].h_dt01->Reset();
        fHist.station[is].panel[ip].ch[ich].h_dtref->Reset();
      }
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// assume several similar jobs
//-----------------------------------------------------------------------------
int plot_n002_hist::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");
  fBook->SaveFolder(fTopFolder,f);
  f->Close();
  delete f;

  return 0;
}

#endif


//-----------------------------------------------------------------------------
int make_calib_hist(int Run1, int Run2, int Station1, int Station2) {
  int k = 0;
  
  plot_n002_hist* p[100]; // not more than 8 calibration runs in a series

  int i = 0;
  for (int run=Run1; run<Run2+1; run++) {
    std::cout << std::format("FIX REFERENCE CHANNELS\n");
    return -1;
    int ref_channel = 11;
    p[i] = new plot_n002_hist("n002",run,ref_channel,Station1,Station2);
    p[i]->Loop();
    i++;
  }

  std::string fn = std::format("pulse_injection_{:6d}_{:6d}.hist",Run1,Run2);
  p[0]->SaveHistograms(fn.data());
  
  return 0;
}


//-----------------------------------------------------------------------------
// clone fit results into a global structure which will be used
// to produce the final calibraton table
//-----------------------------------------------------------------------------
int plot_n002_hist::CloneFitResults () {
  int rc(0);

  for (int is=fStation1; is<=fStation2; is++) {
    for (int ip=0; ip<kNPanelsPerStation; ip++) {
      fgFitResult[is][ip].ref_channel = fFitResult[is][ip].ref_channel;
      for (int ich=0; ich<96; ich++) {
        ChannelFitResult_t* cfr1 = &fFitResult[is][ip].cfr[ich];
        ChannelFitResult_t* cfr2 = &fgFitResult[is][ip].cfr[ich];
        cfr2->chi2dof = cfr1->chi2dof;
        for (int i=0; i<3; i++) {
          cfr2->p[i] = cfr1->p[i];
          cfr2->e[i] = cfr1->e[i];
        }
      }
    }
  }
  return rc;
}

//-----------------------------------------------------------------------------
// clone fit results into a global structure which will be used
// to produce the final calibraton table
// fit results stored in Ph to this
//-----------------------------------------------------------------------------
int plot_n002_hist::AddFitResults (plot_n002_hist* Ph, int Station, int Panel) {
  int rc(0);
  
  int is1(Station), is2(Station+1);
  if (Station == -1) { is1=fStation1; is2 = fStation2; }
  
  int ip1(Panel), ip2(Panel+1);
  if (Panel == -1) { ip1=0; ip2 = kNPanelsPerStation; }

  for (int is=is1; is<=is2; is++) {
    for (int ip=ip1; ip<ip2; ip++) {
      float ref_dt_offset(0);
      int ref_channel_ph = Ph->fFitResult[is][ip].ref_channel;
      int ref_channel    = fFitResult[is][ip].ref_channel;
      if (ref_channel != ref_channel_ph) {
//-----------------------------------------------------------------------------
// figure out offset to be added
// cfr stores offset of the Ph reference channel 
//-----------------------------------------------------------------------------
        ChannelFitResult_t* cfr = &fFitResult[is][ip].cfr[ref_channel_ph];
        float ref_dt_offset     = cfr->p[1]; // dt(ch9 -ch12) - to be added to all added fit means
        
        std::cout << std::format("is:{} ip:{} ref_channel_ph:{} ref_dt_offset:{}\n",
                                 is,ip,ref_channel_ph,ref_dt_offset);
      }
      
      for (int ich=0; ich<96; ich++) {
        ChannelFitResult_t* cfr = &fFitResult[is][ip].cfr[ich];
//-----------------------------------------------------------------------------
// skip channels with calibrations already defined and Ph channels w/o fits
//-----------------------------------------------------------------------------
        if (cfr->chi2dof > 0)                               continue;
        
        ChannelFitResult_t* cfr_ph = &Ph->fFitResult[is][ip].cfr[ich];
        if (cfr_ph->chi2dof < 0)                            continue;
//-----------------------------------------------------------------------------
// account for potential offset
//-----------------------------------------------------------------------------  
        cfr->chi2dof = cfr_ph->chi2dof;
        for (int i=0; i<3; i++) {
          cfr->p[i] = cfr_ph->p[i];
          cfr->e[i] = cfr_ph->e[i];
        }
        cfr->p[1]  += ref_dt_offset;
      }
    }
  }
  return rc;
}
