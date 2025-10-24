/////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////
#include "Stntuple/loop/TStnAna.hh"

#include "daqana/ana/rocsim.hh"
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
    -2.129703, 0.178027,-3.811677, 0.861855,-1.374878,-1.399761,-1.514854, 1.019253,-2.020230, 1.422935,
     1.143216, 0.171659,-2.541716, 1.047467,-4.131723, 1.406507,-0.261447,-0.075517,-0.551845, 0.212474,
    -2.377313, 0.391070, 0.057158,-0.705488,-2.450838,-0.219201,-2.419693, 1.790118,-1.148296,-0.839230,
    -1.183606, 2.331864,-2.183257, 2.073913, 0.113940,-0.698815,-2.883949, 1.185830,-3.837072, 0.980327,
     0.277410, 0.417103,-0.910458, 0.805904,-2.530011,-0.009944,-0.495195,-0.823185,-2.832128,-0.423499,
    -0.241145, 2.091442,-1.690023,-4.442818,-0.840601, 1.248325, 0.510533, 1.291149, 0.166627,-2.742950,
    -3.174086, 1.217493,-0.199138, 1.109670,-0.777815,-3.697247,-1.699748, 1.038269,-0.899651, 0.898945,
    -0.813779,-2.155251,-2.051586,-0.349616,-0.500565, 2.126341,-1.593504,-3.382061,-0.907907, 3.054471,
     0.270250, 1.368577, 0.024238,-2.378631,-3.113740, 1.135546,-0.086879, 1.400219,-0.417473,-4.387886,
    -2.371000, 0.0000  ,-0.676306, 0.663073, 0.000   ,-2.170244
   };

  for (int i=0; i<kNChannels; i++) {
    fAdcIndex[i]              = adc_index[i];
    fInverseMap[adc_index[i]] = i;

    fChannel[i].fNHits           = 0;
    fChannel[i].fNReadoutHits    = 0;
    fChannel[i].fNReadoutHitsTot = 0;

   fChOffset[i] = ch_offset[i]*1.0e-9; // all in seconds

  }
 
  printf("};");
  fFreq0            = 31.29e6;

  InitRun(RunNumber);

  if      (fF0 ==  60) fGenFreq = fFreq0/(pow(2,9)+1); // exact value of the frequency, sec
  else if (fF0 == 250) fGenFreq = fFreq0/(pow(2,7)+1);
  
  fDeltaT           = 1./fGenFreq;         // distance between the pulses
  fEventWindow      = 25.e-9*fHb;          // 25 microseconds in seconds
  fTimeBin          = (5./256)*1.e-9;      // TDC time bin, sec
  fMaxNHitsPerEvent = 255;
  fRefChannel       = 42;

  fFolder           = new TFolder("rocsim","rocsim");
  fHistFolder       = fFolder->AddFolder("Hist","Hist");
}

//-----------------------------------------------------------------------------
TRocSim::~TRocSim() {
  delete fFolder;
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
  else if (RunNumber == 1){
    fF0=60;
    fHb=2000;
    fFpgaOffset=0;
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::BookHistograms() {

  fHist.fEvent.fNHitsTotVsIch         = new TH1F(Form("nh_vs_ich"),"nhits vs ich",100,0,100);
  fHist.fEvent.fNHitsTotVsAdc         = new TH1F(Form("nh_vs_adc"),"nhits vs adc",100,0,100);
  fHist.fEvent.fNReadoutHitsTot       = new TH1F(Form("nrh_tot"  ),"nrhits tot"  ,500,0,500);
  fHist.fEvent.fNHits_fpga1           = new TH1F("nhits_fpga1","nhitsfpga1",  20,0,20);
  fHist.fEvent.fNHits_fpga2           = new TH1F("nhits_fpga2","nhitsfpga2",  20,0,20);
  fHist.fEvent.fNHitsTot_fpga1        = new TH1F("nhits_fpga1_tot","nhitsfpga1_tot",  500,0,500);
  fHist.fEvent.fNHitsTot_fpga2        = new TH1F("nhits_fpga2_tot","nhitsfpga2_tot",  500,0,500);
  fHist.fEvent.fNReadoutHitsTot_fpga1 = new TH1F("nhitsread_fpga1_tot","nhitsfpga1_tot",  500,0,500);
  fHist.fEvent.fNReadoutHitsTot_fpga2 = new TH1F("nhitsread_fpga2_tot","nhitsfpga2_tot",  500,0,500);

  fHistFolder->Add(fHist.fEvent.fNHitsTotVsIch);
  fHistFolder->Add(fHist.fEvent.fNHitsTotVsAdc);
  fHistFolder->Add(fHist.fEvent.fNReadoutHitsTot);
  fHistFolder->Add(fHist.fEvent.fNHits_fpga1);
  fHistFolder->Add(fHist.fEvent.fNHits_fpga2);
  fHistFolder->Add(fHist.fEvent.fNHitsTot_fpga1);
  fHistFolder->Add(fHist.fEvent.fNHitsTot_fpga2);
  fHistFolder->Add(fHist.fEvent.fNReadoutHitsTot_fpga1);
  fHistFolder->Add(fHist.fEvent.fNReadoutHitsTot_fpga2);

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
  int channels_1_fpga[48]={91,85,79,73,67,61,55,49,43,37,31,25,19,13,7,1,90,84,78,72,66,60,54,48,
			   42,36,30,24,18,12,6,0,93,87,81,75,69,63,57,51,45,39,33,27,21,15,9,3};
  int channels_2_fpga[48]={44, 38, 32, 26, 20, 14,8,2,92, 86, 80, 74, 68, 62, 56, 50,47, 41,35,29, 23, 17, 11, 5,
			   95, 89, 83, 77, 71, 65, 59, 53,46, 40, 34, 28, 22, 16, 10, 4,94, 88, 82, 76, 70, 64, 58, 52};
  double tref = 1.e12;

  if (fChannel[fRefChannel].fNHits > 0) tref = fChannel[fRefChannel].fHit[0].fTime;
  
  for (int ich=0; ich<kNChannels; ich++) {
    TRocChannel* ch = &fChannel[ich];

    int nh = ch->fNReadoutHits;
    fHist.fChannel[ich].fNHits->Fill(nh);

    for(int vv=0;vv<48;vv++){
     
      if(adc_index[ich]==channels_1_fpga[vv]){
	   
	fHist.fEvent.fNHits_fpga1->Fill(ch->fNReadoutHits_fpga1);
	break;
      }else if(adc_index[ich]==channels_2_fpga[vv]){
	fHist.fEvent.fNHits_fpga2->Fill(ch->fNReadoutHits_fpga2);
	break;
      }
    }
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
     
  if(fNHitsEvent >255){
    printf("EventNumber=%d NumerofHits=%d\n",fEventNumber,fNHitsEvent);
    
  }
   fHist.fEvent.fNReadoutHitsTot->Fill(fNHitsEvent);
   fHist.fEvent.fNReadoutHitsTot_fpga1->Fill(fNHitsEvent_fpga1);
   fHist.fEvent.fNReadoutHitsTot_fpga2->Fill(fNHitsEvent_fpga2);
 
  
  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::SimulateEvent(int EventNumber) {
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
//-----------------------------------------------------------------------------
// 'i' is the channel number in the readout order, the first channel to be 
// readout is channel adc_index[0] = 91
//-----------------------------------------------------------------------------
  for (int i=0; i<kNChannels; i++) {
    TRocChannel* ch   = &fChannel[adc_index[i]];
    ch->fNHits        = 0;
    ch->fNReadoutHits = 0;
//-----------------------------------------------------------------------------
// now see how many pulses fit within the window
//-----------------------------------------------------------------------------
    double time = t0 + fChOffset[adc_index[i]];
    if (i >= 48) time  = t1 + fChOffset[adc_index[i]];
    
    if (time < 0) {
      time +=fDeltaT;
    }
    if (time > fEventWindow){
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
  fNHitsEvent_fpga1=0;
  fNHitsEvent_fpga2=0;
  int channels_fpga1[48]={91,85,79,73,67,61,55,49,43,37,31,25,19,13,7,1,90,84,78,72,66,60,54,48,42,36,30,24,18,12,6,0,93,87,81,75,69,63,57,51,45,39,33,27,21,15,9,3};
  int channels_fpga2[48]={44, 38, 32, 26, 20, 14,8,2,92, 86, 80, 74, 68, 62, 56, 50,47, 41,35,29, 23, 17, 11, 5,95, 89, 83, 77, 71, 65, 59, 53,46, 40, 34, 28, 22, 16, 10, 4,94, 88, 82, 76, 70, 64, 58, 52};
  for (int i=0; i<kNChannels; i++) {
    TRocChannel* ch = &fChannel[fAdcIndex[i]];
    
    int nh = fNHitsEvent + ch->fNHits;
   
    //printf("nh=%d\n", nh);
    int nh_fpga1=fNHitsEvent_fpga1;
    int nh_fpga2=fNHitsEvent_fpga2;
    for(int vv=0;vv<48;vv++){
      if(adc_index[i]==channels_fpga1[vv]){
	nh_fpga1=nh_fpga1+ch->fNHits;
	break;
      }else if(adc_index[i]==channels_fpga2[vv]){
	nh_fpga2=nh_fpga2+ch->fNHits;
	break;
      }
    }
    if (nh > fMaxNHitsPerEvent) {//nh>255
      // printf("max=%d\n", fMaxNHitsPerEvent);
      ch->fNReadoutHits =fMaxNHitsPerEvent-fNHitsEvent;
      for(int vv=0;vv<48;vv++){
	if(adc_index[i]==channels_fpga1[vv]){
	  ch->fNReadoutHits_fpga1=fMaxNHitsPerEvent-fNHitsEvent;
	  break;
	}else if(adc_index[i]==channels_fpga2[vv]){
	  ch->fNReadoutHits_fpga2=fMaxNHitsPerEvent-fNHitsEvent;
	  break;
	}
      }
    }
    else {//nh in questo caso e' 232
//-----------------------------------------------------------------------------
// all hits in this channel are included into readout
//-----------------------------------------------------------------------------
      ch->fNReadoutHits     = ch->fNHits;
      for(int vv=0;vv<48;vv++){
	if(adc_index[i]==channels_fpga1[vv]){
	  ch->fNReadoutHits_fpga1=ch->fNHits;
	  break;
	}else if(adc_index[i]==channels_fpga2[vv]){
	  ch->fNReadoutHits_fpga2=ch->fNHits;
	  break;
	}
      }
    }
   
      ch->fNReadoutHitsTot += ch->fNReadoutHits;
      fNHitsEvent          += ch->fNReadoutHits;
    
    //printf("nhitsevent=%d\n",fNHitsEvent);
    // printf("chfNreadouthits=%d\n",ch->fNReadoutHitsTot);
     for(int vv=0;vv<48;vv++){
	if(adc_index[i]  ==channels_fpga1[vv]){
	  ch->fNReadoutHitsTot_fpga1+=ch->fNReadoutHits_fpga1;
	  fNHitsEvent_fpga1          += ch->fNReadoutHits_fpga1;
	  break;
	}else if(adc_index[i]==channels_fpga2[vv]){
	  ch->fNReadoutHitsTot_fpga2+=ch->fNReadoutHits_fpga2;
	  fNHitsEvent_fpga2          += ch->fNReadoutHits_fpga2;
	  break;
	}
	
     }
    
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::Run(int NEvents) {
  fEventNumber=0;
  BookHistograms();
  
  for (int i=0; i<NEvents; i++) {
    fEventNumber=i;
    SimulateEvent(i);
    FillHistograms();

  }

  return 0;
}

//-----------------------------------------------------------------------------
int TRocSim::SaveHistograms(const char* Filename) {

  TFile* f = new TFile(Filename,"recreate");
  TStnAna::SaveFolder(fFolder,f);
  f->Close();
  
  return 0;
}
