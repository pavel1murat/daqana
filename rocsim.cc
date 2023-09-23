//

#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFolder.h"
#include "TFile.h"

//-----------------------------------------------------------------------------
class THit : public TObject {
public:
  long int fTdc;
  double   fTime;
};

//-----------------------------------------------------------------------------
class TRocChannel : public TObject {
public:
  THit       fHit[20];                  // in a given event
  int        fNHits;                    // in a given event
  int        fNReadoutHits;             // in a given event
  int        fNReadoutHitsTot;          // total number of hits readout in a channel
};


//-----------------------------------------------------------------------------
class TRocSim: public TNamed {
public:
  enum { kNChannels = 96 } ;

  struct ChannelHist_t {
    TH1F* fNHits;
    TH1F* fTime;
    TH1F* fDt;                          // fTime - fTime(ref_channel)
  };

  struct EventHist_t {
    TH1F* fNReadoutHitsTot;
    TH1F* fNHitsTotVsIch;
    TH1F* fNHitsTotVsAdc;
  };


  struct Hist_t {
    ChannelHist_t   fChannel[kNChannels] ;
    EventHist_t     fEvent;
  } fHist;
  
  int           fAdcIndex  [kNChannels];
  int           fInverseMap[kNChannels];
  TRocChannel   fChannel   [kNChannels];
  double        fChOffset  [kNChannels];

  int           fF0;                    // run number-dependent
  int           fHb;
  
  double        fEventWindow;
  double        fFreq0;                 // oscillator frequency
  double        fGenFreq;               // generator frequency
  double        fDeltaT;                // time distance between the two generated pulses
  double        fTimeBin;
  double        fFpgaOffset;

  int           fRefChannel;            // reference channel

  TRandom3      fRn3;

  int           fMaxNHitsPerEvent;      // ROC buffer size
  int           fNHitsEvent;

  TFolder*      fHistFolder;

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TRocSim(const char* Name, int RunNumber);
  ~TRocSim();

  int  InitRun(int RunNumber);
  
  int Run(int NEvents = 100);

  int BookHistograms();
  int FillHistograms();
  int SimulateEvent (int EventNumber);
  int SaveHistograms(const char* Filename);
  
};




//-----------------------------------------------------------------------------
TRocSim::TRocSim(const char* Name, int RunNumber) : TNamed(Name,Name) {
  int adc_index[96] = {
    91, 85, 79, 73, 67, 61, 55, 49,
    43, 37, 31, 25, 19, 13,  7,  1,
    90, 84, 78, 72, 66, 60, 54, 48,
    
    42, 36, 30, 24, 18, 12,  6,  0,
    93, 87, 81, 75, 69, 63, 57, 51,
    45, 39, 33, 27, 21, 15,  9,  3,
    
    44, 38, 32, 26, 20, 14,  8,  2, 
    92, 86, 80, 74, 68, 62, 56, 50,
    47, 41, 35, 29, 23, 17, 11,  5,
      
    95, 89, 83, 77, 71, 65, 59, 53,
    46, 40, 34, 28, 22, 16, 10,  4,
    94, 88, 82, 76, 70, 64, 58, 52
  };
//-----------------------------------------------------------------------------
// offsets, in ns, wrt the first readout channel in a given FPGA
// determined using run 105038
// need to cross check their stability
// an offsetfor a fiven 'i' is the offset of a given channel wrt the channel
// read out first in the corresponding FPGA (first 48 channels - CAL, channels 48-95 - DIGI)
//-----------------------------------------------------------------------------
  double ch_offset[96] = {
    -1.91142,  0.83799, -3.60547,  2.07482, -1.33494, -1.40625, -0.85195,  1.34135, -2.98527,  1.99266, 
     1.09860,  0.70164, -1.91974,  1.16734, -4.06390,  2.50224, -0.00802, -0.24147,  0.25977,  0.46001, 
    -2.18188,  0.52085,  0.03532, -0.49709, -2.23997, -0.29669, -2.40610,  2.24535, -0.49264, -0.33865, 
    -0.88099,  3.30617, -2.13120,  1.84080,  0.87094, -0.29687, -2.52983,  1.70986, -3.77760,  0.73012, 
     0.62384,  0.47639, -0.21883,  1.27542, -3.12222,  0.26359, -0.26839, -1.11142, -2.26955,  0.04862, 
     0.06483,  2.21114, -1.47437, -4.15237, -0.45083,  2.03541,  1.28405,  1.54520, -0.11049, -2.22120, 
    -3.48292,  1.72673, -0.28377,  2.06617, -0.09606, -4.40612, -2.03001,  2.13547, -0.53361,  2.34388, 
    -0.43629, -2.55260, -1.43804,  0.27971, -0.91122,  2.65964, -1.91070, -3.50933, -0.34008,  3.26338, 
     0.09934,  1.74737,  0.86814, -2.80922, -2.51025,  1.22813,  0.28790,  1.95993, -0.30327, -4.36486, 
    -2.34925,  0.00000, -0.34988,  1.59718,  0.00000, -2.29823
  };

  for (int i=0; i<kNChannels; i++) {
    fAdcIndex[i]              = adc_index[i];
    fInverseMap[adc_index[i]] = i;

    fChannel[i].fNHits           = 0;
    fChannel[i].fNReadoutHits    = 0;
    fChannel[i].fNReadoutHitsTot = 0;

    fChOffset[i] = ch_offset[i]*1.0e-9; // all in seconds
  }

  fFreq0            = 31.29e6;

  InitRun(RunNumber);

  if      (fF0 ==  60) fGenFreq = fFreq0/(pow(2,9)+1); // exact value of the frequency, sec
  else if (fF0 == 250) fGenFreq = fFreq0/(pow(2,7)+1);
  
  fDeltaT           = 1./fGenFreq;         // distance between the pulses
  fEventWindow      = 25.e-9*fHb;          // 25 microseconds in seconds
  fTimeBin          = (5./256)*1.e-9;      // TDC time bin, sec
  fMaxNHitsPerEvent = 255;
  fRefChannel       = 42;

  fHistFolder       = new TFolder("Hist","Hist");
}

//-----------------------------------------------------------------------------
TRocSim::~TRocSim() {
  delete fHistFolder;
}

//-----------------------------------------------------------------------------
int TRocSim::InitRun(int RunNumber) {

  if (RunNumber == 281) {
    fF0         = 60;
    fHb         = 2000;
    fFpgaOffset = 15030.e-9;
  }
  else if (RunNumber  == 105023) {
    fF0         = 250;
    fHb         = 1000;
    fFpgaOffset = 0.;   // can't tell - nothing seen on HV side
  }
  else if (RunNumber  == 105026) {
    fF0         = 60;
    fHb         = 2000;
    fFpgaOffset = 10400.e-9;
  }
  else if (RunNumber  == 105038) {
    fF0         = 60;
    fHb         = 1000;
    fFpgaOffset = 7965.e-9;
  }
  else if (RunNumber  == 105041) {
    fF0 = 60;
    fHb = 1400;
    fFpgaOffset = 1128.e-9;
  }
  else if (RunNumber  == 105042) {
    fF0 = 60;
    fHb = 1600;
    fFpgaOffset = 1128.e-9;
  }
  else if (RunNumber  == 105043) {
    fF0 = 60;
    fHb = 2400;
    fFpgaOffset = 1128.e-9;
  }
  else if (RunNumber  == 105044) {
    fF0 = 60;
    fHb = 2200;
    fFpgaOffset = 1128.e-9;
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::BookHistograms() {

  fHist.fEvent.fNHitsTotVsIch    = new TH1F(Form("nh_vs_ich"),"nhits vs ich",100,0,100);
  fHist.fEvent.fNHitsTotVsAdc    = new TH1F(Form("nh_vs_adc"),"nhits vs adc",100,0,100);
  fHist.fEvent.fNReadoutHitsTot  = new TH1F(Form("nrh_tot"  ),"nrhits tot"  ,500,0,500);

  fHistFolder->Add(fHist.fEvent.fNHitsTotVsIch);
  fHistFolder->Add(fHist.fEvent.fNHitsTotVsAdc);
  fHistFolder->Add(fHist.fEvent.fNReadoutHitsTot);
  
  for (int i=0; i<kNChannels; i++) {
    fHist.fChannel[i].fNHits = new TH1F(Form("nh_%02i"  ,i),Form("nhits ch # %02i",i),  20,0,20);
    fHist.fChannel[i].fTime  = new TH1F(Form("time_%02i",i),Form("time  ch # %02i",i),1000,0,fEventWindow);
    fHist.fChannel[i].fDt    = new TH1F(Form("dt_%02i"  ,i),Form("dt    ch # %02i",i),1000,-fEventWindow/2,fEventWindow/2);

    TFolder* f = fHistFolder->AddFolder(Form("ch_%02i",i),Form("ch_%02i_folder",i));
    f->Add(fHist.fChannel[i].fNHits);
    f->Add(fHist.fChannel[i].fTime);
    f->Add(fHist.fChannel[i].fDt);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::FillHistograms() {

  double tref = 1.e12;
  
  if (fChannel[fRefChannel].fNHits > 0) tref = fChannel[fRefChannel].fHit[0].fTime;
  
  for (int ich=0; ich<kNChannels; ich++) {
    TRocChannel* ch = &fChannel[ich];
    int nh = ch->fNReadoutHits;
    fHist.fChannel[ich].fNHits->Fill(nh);

    fHist.fEvent.fNHitsTotVsIch->Fill(ich,nh);
    int adc = fInverseMap[ich];
    fHist.fEvent.fNHitsTotVsAdc->Fill(adc,nh);

    for (int ih=0;ih<nh; ih++) {
      double time = ch->fHit[ih].fTime;
      fHist.fChannel[ich].fTime->Fill(time);

      if (ih == 0) {
        fHist.fChannel[ich].fDt->Fill(time-tref);
      }
    }
  }

  fHist.fEvent.fNReadoutHitsTot->Fill(fNHitsEvent);
  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::SimulateEvent(int EventNumber) {

  // simulate hits in all channels, start from assuming that 

//-----------------------------------------------------------------------------
// the time of the first DIGI#0 pulse
// if it is > fEventWindow, no pulses in this channel
//-----------------------------------------------------------------------------
  double t0 = fRn3.Rndm(EventNumber)*fDeltaT;   // 
//-----------------------------------------------------------------------------
// the time of the first DIGI#1 pulse
// fFpgaOffset < fDeltaT by definition
// if it is > fEventWindow:
//-----------------------------------------------------------------------------
  double t1 = t0+fFpgaOffset;
  while (t1-fDeltaT > 0) t1 = t1-fDeltaT;
  
  for (int i=0; i<kNChannels; i++) {
    TRocChannel* ch   = &fChannel[i];
    ch->fNHits        = 0;
    ch->fNReadoutHits = 0;
//-----------------------------------------------------------------------------
// now see how many pulses fit within the window
//-----------------------------------------------------------------------------
    double time = t0 + fChOffset[i];
    if (i >= 48) time  = t1;
    
    if (time < 0) {
      // no pulses in his channel for this event
      continue;
    }
//-----------------------------------------------------------------------------
// calculate the number of pulses in this event
//-----------------------------------------------------------------------------
    while (time < fEventWindow) {
      ch->fHit[ch->fNHits].fTime = time;
      ch->fNHits                += 1;
      time                      += fDeltaT;
    }
  }
//-----------------------------------------------------------------------------
// simulate the readout
// loop over channels in readout sequence
//-----------------------------------------------------------------------------
  fNHitsEvent =  0;
  for (int i=0; i<kNChannels; i++) {
    TRocChannel* ch = &fChannel[fAdcIndex[i]];
    int nh = fNHitsEvent + ch->fNHits;
    if (nh > fMaxNHitsPerEvent) {
      ch->fNReadoutHits = fMaxNHitsPerEvent-fNHitsEvent;
    }
    else {
//-----------------------------------------------------------------------------
// all hits in this channel are included into readout
//-----------------------------------------------------------------------------
      ch->fNReadoutHits     = ch->fNHits;
    }

    ch->fNReadoutHitsTot += ch->fNReadoutHits;
    fNHitsEvent          += ch->fNReadoutHits;
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::Run(int NEvents) {

  BookHistograms();
  
  for (int i=0; i<NEvents; i++) {
    SimulateEvent(i);
    FillHistograms();
  }

  return 0;
}


//-----------------------------------------------------------------------------
int TRocSim::SaveHistograms(const char* Filename) {
  TFile* f = new TFile(Filename,"recreate");

  // cd happens by default

  fHistFolder->Write();

  f->Close();
  
  return 0;
}

