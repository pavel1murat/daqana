//------------------------------------------------------------------------------
// example:
//
//  TRocSim * rs=new TRocSim("a",281);
//  rs->Run(1000);
//-----------------------------------------------------------------------------
#ifndef __daqana_ana_rocsim_hh__
#define __daqana_ana_rocsim_hh__

#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFolder.h"
#include "TFile.h"

//-----------------------------------------------------------------------------
class TRocSim: public TNamed {
public:
enum { kNChannels = 96 } ;

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
    int        fNReadoutHits_fpga1;             // in a given event
    int        fNReadoutHits_fpga2;             // in a given event
    int        fNReadoutHitsTot_fpga1;            
    int        fNReadoutHitsTot_fpga2;             
    int        fNReadoutHitsTot;          // total number of hits readout in a channel
  };

  struct ChannelHist_t {
    TH1F* fNHits;
    TH1F* fTime;
    TH1F* fDt;                          // fTime - fTime(ref_channel)
  };

  struct EventHist_t {
    TH1F* fNReadoutHitsTot;
    TH1F* fNHits_fpga1;
    TH1F* fNHits_fpga2;
    TH1F* fNHitsTot_fpga1;
    TH1F* fNHitsTot_fpga2;
    TH1F* fNReadoutHitsTot_fpga1;
    TH1F* fNReadoutHitsTot_fpga2;
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
  int           fNHitsEvent_fpga1;
  int           fNHitsEvent_fpga2;
  int           fEventNumber;
  TFolder*      fFolder;
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

#endif
