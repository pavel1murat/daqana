///////////////////////////////////////////////////////////////////////////////
// tracker DQM module
///////////////////////////////////////////////////////////////////////////////
#include "TRACE/tracemf.h"
#define  TRACE_NAME "TrackerDQM"

#include "TSystem.h"
#include "TROOT.h"

#include "daqana/mod/TrackerDQM_module.hh"
//-----------------------------------------------------------------------------
unsigned int reverseBits(unsigned int num) {
  unsigned int numOfBits = 10; // sizeof(num) * 8; // Number of bits in an unsigned int

  unsigned int reversedNum = 0;
  for (unsigned int i = 0; i < numOfBits; ++i) {
    if ((num & (1 << i)) != 0)
      reversedNum |= 1 << ((numOfBits - 1) - i);
  }

  return reversedNum;
}

//-----------------------------------------------------------------------------
unsigned int correctedTDC(unsigned int TDC) {
  uint32_t corrected_tdc = ((TDC & 0xFFFF00) + (0xFF  - (TDC & 0xFF)));
  return corrected_tdc;
}

//-----------------------------------------------------------------------------
// NWords : number of short (16-bit) words
//-----------------------------------------------------------------------------
void TrackerDQM::printFragment(const artdaq::Fragment* Frag, int NWords) {
//-----------------------------------------------------------------------------
// print fragments in HEX, for the tracker, the data has to be in 2-byte words
//-----------------------------------------------------------------------------
  ushort* buf = (ushort*) (Frag->dataBegin());

  int loc     = 0;
      
  for (int i=0; i<NWords; i++) {
    if (loc == 0) printf(" 0x%08x: ",i*2);

    ushort  word = buf[i];
    printf("0x%04x ",word);
    
    loc += 1;
    if (loc == 8) {
      printf("\n");
      loc = 0;
    }
  }
      
  if (loc != 0) printf("\n");
}

//-----------------------------------------------------------------------------
TrackerDQM::TrackerDQM(art::EDAnalyzer::Table<Config> const& conf) : 
  art::EDAnalyzer   (conf),

  _diagLevel       (conf().diagLevel       ()), 
  _minNBytes       (conf().minNBytes       ()), 
  _maxNBytes       (conf().maxNBytes       ()), 
  _dataHeaderOffset(conf().dataHeaderOffset()), 
  _activeLinks     (conf().activeLinks     ()),
  _refChCal        (conf().refChCal        ()),
  _refChHV         (conf().refChHV         ()),
  _dumpDTCRegisters(conf().dumpDTCRegisters()),
  _analyzeFragments(conf().analyzeFragments()),
  _maxFragmentSize (conf().maxFragmentSize ()),
  _pulserFrequency (conf().pulserFrequency ()),
  _nSamplesBL      (conf().nSamplesBL      ()),
  _minPulseHeight  (conf().minPulseHeight  ()),
  _port            (conf().port            ())
{
  // conf_             (conf()),
  // histType_         (conf().histType ()),
  // diagLevel_        (conf().diag     ()),
//-----------------------------------------------------------------------------
// for now, assume only one station, but implement data structures handling 
// full tracker
//-----------------------------------------------------------------------------
  _station = 0;
  _plane   = 0;

  double f0(31.29e6);                   // 31.29 MHz

  _timeWindow = conf().timeWindow()*25.;  // in ns

  if      (_pulserFrequency ==  60) _freq = f0/(pow(2,9)+1);     // ~ 60 kHz
  else if (_pulserFrequency == 250) _freq = f0/(pow(2,7)+1);     // ~250 kHz

  _dt   = 1/_freq*1.e9;               // in ns
//------------------------------------------------------------------------------
// default map, Richie says TS1 may have an old firmware with some bugs
//-----------------------------------------------------------------------------
  int adc_index_0[96] = {
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
// TS1: compared to _0, swap lanes 3 and 4, and in the new lane3 swap lines 1 and 3
//-----------------------------------------------------------------------------
  int adc_index_1[96] = {
    91, 85, 79, 73, 67, 61, 55, 49,
    43, 37, 31, 25, 19, 13,  7,  1,
    90, 84, 78, 72, 66, 60, 54, 48,
      
    42, 36, 30, 24, 18, 12,  6,  0,
    93, 87, 81, 75, 69, 63, 57, 51,
    45, 39, 33, 27, 21, 15,  9,  3,
      
    94, 88, 82, 76, 70, 64, 58, 52,
    46, 40, 34, 28, 22, 16, 10,  4,
    95, 89, 83, 77, 71, 65, 59, 53,
      
    44, 38, 32, 26, 20, 14,  8,  2, 
    92, 86, 80, 74, 68, 62, 56, 50,
    47, 41, 35, 29, 23, 17, 11,  5
  };


  for (int i=0; i<96; i++) {
    _adc_index_0[adc_index_0[i]] = i;
    _adc_index_1[adc_index_1[i]] = i;
  }

//-----------------------------------------------------------------------------
// initialize reference channels, at this point use channels 91 and 94 for all 
// ROC's (the readout order is defined in firmware and is the same for all channels)
//-----------------------------------------------------------------------------
  _nActiveLinks        = _activeLinks.size();

  for (int roc=0; roc<kMaxNLinks; roc++) {
    _referenceChannel[roc][0] = 91;
    _referenceChannel[roc][1] = 94;
    if (roc < _nActiveLinks) {
      _referenceChannel[roc][0] = _refChCal[roc];
      _referenceChannel[roc][1] = _refChHV [roc];
    }

    _event_data.rdata[_station][_plane][roc].ref_ch[0] = &_event_data.rdata[_station][_plane][roc].channel[_referenceChannel[roc][0]];
    _event_data.rdata[_station][_plane][roc].ref_ch[1] = &_event_data.rdata[_station][_plane][roc].channel[_referenceChannel[roc][1]];
  }

  _tdc_bin             = (5/256.*1e-3);       // TDC bin width (Richie), in us
  _tdc_bin_ns          = _tdc_bin*1e3;        // convert to ns

  _initialized         = 0;

  _plot_wf.link    = -1;
  _plot_wf.channel = -1;
  int nch_plot = _plotWaveforms.size();
  if (nch_plot == 2) {
    _plot_wf.link    = _plotWaveforms[0];
    _plot_wf.channel = _plotWaveforms[1];
  }

}

//-----------------------------------------------------------------------------
// I : channel number
//-----------------------------------------------------------------------------
void TrackerDQM::book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int I) {

  Hist->nhits   = Dir->make<TH1F>(Form("ch_%02i_nhits"  ,I),Form("run %06i: ch %02i nhits"  ,RunNumber,I), 20, -0.5, 19.5);
    
  Hist->time[0] = Dir->make<TH1F>(Form("ch_%02i_time0"  ,I),Form("run %06i: ch %02i time[0]",RunNumber,I),1000, 0., 100.);  // us
  Hist->time[1] = Dir->make<TH1F>(Form("ch_%02i_time1"  ,I),Form("run %06i: ch %02i time[0]",RunNumber,I),1000, 0., 100.);  // us

  Hist->t0  [0] = Dir->make<TH1F>(Form("ch_%02i_t0_0"   ,I),Form("run %06i: ch %02i t0[0]  ",RunNumber,I),1000,-20., 80.);  // ns
  Hist->t0  [1] = Dir->make<TH1F>(Form("ch_%02i_t0_1"   ,I),Form("run %06i: ch %02i t0[1]  ",RunNumber,I),1000,-20., 80.);  // ns

  Hist->t1  [0] = Dir->make<TH1F>(Form("ch_%02i_t1_0"   ,I),Form("run %06i: ch %02i t1[0]  ",RunNumber,I),1000,-20., 80.);  // ns
  Hist->t1  [1] = Dir->make<TH1F>(Form("ch_%02i_t1_1"   ,I),Form("run %06i: ch %02i t1[1]  ",RunNumber,I),1000,-20., 80.);  // ns
    
  Hist->tot [0] = Dir->make<TH1F>(Form("ch_%02i_tot0"   ,I),Form("run %06i: ch %02i tot[0]" ,RunNumber,I), 100, 0., 100.);
  Hist->tot [1] = Dir->make<TH1F>(Form("ch_%02i_tot1"   ,I),Form("run %06i: ch %02i tot[1]" ,RunNumber,I), 100, 0., 100.);
  Hist->pmp     = Dir->make<TH1F>(Form("ch_%02i_pmp"    ,I),Form("run %06i: ch %02i pmp"    ,RunNumber,I), 100, 0.,  10.);
  Hist->dt0     = Dir->make<TH1F>(Form("ch_%02i_dt0"    ,I),Form("run %06i: ch %02i T0(i+1)-T0(i)",RunNumber,I)      ,50000,  0.,50);
  Hist->dt1     = Dir->make<TH1F>(Form("ch_%02i_dt1"    ,I),Form("run %06i: ch %02i T1(i+1)-T1(i)",RunNumber,I)      ,50000,  0.,50);
  Hist->dt2     = Dir->make<TH1F>(Form("ch_%02i_dt2"    ,I),Form("run %06i: ch %02i T2(i+1)-T2(i)",RunNumber,I)      ,50000,  0.,50);
  Hist->dt0r    = Dir->make<TH1F>(Form("ch_%02i_dt0r_0" ,I),Form("run %06i: ch %02i T0(ich,0)-T0(ref,0)[0]",RunNumber,I),20000,-10.,10);
  Hist->dt1r    = Dir->make<TH1F>(Form("ch_%02i_dt1r_0" ,I),Form("run %06i: ch %02i T1(ich,0)-T1(ref,0)[0]",RunNumber,I),20000,-10.,10);
//-----------------------------------------------------------------------------
// waveform parameters
//-----------------------------------------------------------------------------
  Hist->fsample = Dir->make<TH1F>(Form("ch_%02i_fs"    ,I),Form("run %06i: ch %02i first sample"   ,RunNumber,I), 30,-0.5, 29.5);
  Hist->tmean   = Dir->make<TH1F>(Form("ch_%02i_tm"    ,I),Form("run %06i: ch %02i tmean"          ,RunNumber,I),1500,  0, 30  );
  Hist->bline   = Dir->make<TH1F>(Form("ch_%02i_bl"    ,I),Form("run %06i: ch %02i WF baseline"    ,RunNumber,I),200,0,300);
  Hist->pheight = Dir->make<TH1F>(Form("ch_%02i_ph"    ,I),Form("run %06i: ch %02i WF pulse height",RunNumber,I),500,0,500);
  Hist->q       = Dir->make<TH1F>(Form("ch_%02i_q"     ,I),Form("run %06i: ch %02i WF charge"      ,RunNumber,I),500,0,500);//
  Hist->qt      = Dir->make<TH1F>(Form("ch_%02i_qt"    ,I),Form("run %06i: ch %02i WF tail charge" ,RunNumber,I),500,0,500);
  Hist->qtq     = Dir->make<TH1F>(Form("ch_%02i_qtq"   ,I),Form("run %06i: ch %02i WF Qt/Q"        ,RunNumber,I),200,0,1);
//-----------------------------------------------------------------------------
// waveform histograms, assume number of samples < 30
//-----------------------------------------------------------------------------
  for (int j=0; j<kMaxNHitsPerChannel; j++) {
    Hist->wf[j] = Dir->make<TH1F>(Form("wf_ch_%02i_%i",I,j),Form("run %06i: ch [%02i][%i] waveform",RunNumber,I,j),30, 0.,30.);
    Hist->raw_wf[j] = Dir->make<TH1F>(Form("raw_wf_ch_%02i_%i",I,j),Form("run %06i: ch [%02i][%i] raw waveform",RunNumber,I,j),30, 0.,30.);
  }
}


//-----------------------------------------------------------------------------
void TrackerDQM::book_roc_histograms(art::TFileDirectory* Dir, int RunNumber, RocHist_t* Hist, int Link) {
  Hist->nbytes          = Dir->make<TH1F>("nbytes"  , Form("run %06i: n bytes"     ,RunNumber),10000,    0., 10000.);
  Hist->npackets        = Dir->make<TH1F>("npackets", Form("run %06i: n packets"   ,RunNumber), 1000,    0.,  1000.);
  Hist->nhits           = Dir->make<TH1F>("nhits"   , Form("run %06i: n hits"      ,RunNumber),  300,    0.,   300.);
  Hist->valid           = Dir->make<TH1F>("valid"   , Form("run %06i: valid"       ,RunNumber),    2,    0.,     2.);

  Hist->nh_vs_ch        = Dir->make<TH2F>("nh_vs_ch"  , Form("run %06i: nh vs ch"  ,RunNumber),  100,0.,100., 10,0,10);
  Hist->nh_vs_adc1      = Dir->make<TH2F>("nh_vs_adc1", Form("run %06i: nh vs adc" ,RunNumber),  100,0.,100., 10,0,10);

  Hist->dt0r_vs_ch      = Dir->make<TH2F>("dt0r_vs_ch_0", Form("run %06i: dt0r vs ch[0]",RunNumber),  100,0.,100.,2500,-25,25);
  Hist->dt1r_vs_ch      = Dir->make<TH2F>("dt1r_vs_ch_0", Form("run %06i: dt1r vs ch[0]",RunNumber),  100,0.,100.,2500,-25,25);

  Hist->dt0r01          = Dir->make<TH1F>("dt0r01", Form("run %06i: dt0r01"             ,RunNumber), 40000,-20000,20000);
  Hist->dt1r01          = Dir->make<TH1F>("dt1r01", Form("run %06i: dt1r01"             ,RunNumber), 40000,-20000,20000);

  Hist->nhits_vs_ich    = Dir->make<TH1F>("nh_vs_ich"  , Form("run %06i: nh vs ich"    ,RunNumber),  100, 0.,   100.);
  Hist->nhits_vs_adc[0] = Dir->make<TH1F>("nh_vs_adc_0", Form("run %06i: nh vs adc_0",RunNumber),  100, 0.,   100.);
  Hist->nhits_vs_adc[1] = Dir->make<TH1F>("nh_vs_adc_1", Form("run %06i: nh vs adc_1",RunNumber),  100, 0.,   100.);

  Hist->sum_error_vs_ch = Dir->make<TH1F>("sum_err_vs_ch", Form("run %06i: sum err vs ch"     ,RunNumber),  100,0,100);
  Hist->dt0rc_vs_ch[0]  = Dir->make<TH2F>("dt0rc_vs_ch_0", Form("run %06i: dt0rc vs ch[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  Hist->dt0rc_vs_ch[1]  = Dir->make<TH2F>("dt0rc_vs_ch_1", Form("run %06i: dt0rc vs ch[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  Hist->dt1rc_vs_ch[0]  = Dir->make<TH2F>("dt1rc_vs_ch_0", Form("run %06i: dt1rc vs ch[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  Hist->dt1rc_vs_ch[1]  = Dir->make<TH2F>("dt1rc_vs_ch_1", Form("run %06i: dt1rc vs ch[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  Hist->dt0rc_vs_adc[0] = Dir->make<TH2F>("dt0rc_vs_adc_0", Form("run %06i: dt0rc vs adc[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  Hist->dt0rc_vs_adc[1] = Dir->make<TH2F>("dt0rc_vs_adc_1", Form("run %06i: dt0rc vs adc[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  Hist->dt1rc_vs_adc[0] = Dir->make<TH2F>("dt1rc_vs_adc_0", Form("run %06i: dt1rc vs adc[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  Hist->dt1rc_vs_adc[1] = Dir->make<TH2F>("dt1rc_vs_adc_1", Form("run %06i: dt1rc vs adc[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  Hist->fs_vs_ich       = Dir->make<TProfile>("fs_vs_ich"  , Form("run %06i: fs vs ich"    ,RunNumber),  100, 0.,   100.,0,  30);
  Hist->bl_vs_ich       = Dir->make<TProfile>("bl_vs_ich"  , Form("run %06i: bl vs ich"    ,RunNumber),  100, 0.,   100.,0, 500);
  Hist->ph_vs_ich       = Dir->make<TProfile>("ph_vs_ich"  , Form("run %06i: ph vs ich"    ,RunNumber),  100, 0.,   100.,0, 500);
  Hist->q_vs_ich        = Dir->make<TProfile>("q_vs_ich"   , Form("run %06i: Q vs ich"     ,RunNumber),  100, 0.,   100.,0,1500);
  Hist->qt_vs_ich       = Dir->make<TProfile>("qt_vs_ich"  , Form("run %06i: Qt vs ich"    ,RunNumber),  100, 0.,   100.,0, 500);
  Hist->qtq_vs_ich      = Dir->make<TProfile>("qtq_vs_ich" , Form("run %06i: Qt/Q vs ich"  ,RunNumber),  100, 0.,   100.,0, 1);
  
  for (int i=0; i<kNChannels; i++) {
    art::TFileDirectory chan_dir = Dir->mkdir(Form("ch_%02i",i));
    book_channel_histograms(&chan_dir,RunNumber,&Hist->channel[i],i);
  }
}


//-----------------------------------------------------------------------------
void TrackerDQM::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
  Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"      ,RunNumber), 1000,   0.,   1000.);
  Hist->nbtot           = Dir->make<TH1F>("nbtot"      , Form("run %06i: nbytes total"     ,RunNumber), 1000,   0., 100000.);
  Hist->nfrag           = Dir->make<TH1F>("nfrag"      , Form("run %06i: n fragments"      ,RunNumber),  100,   0.,    100.);
  Hist->fsize           = Dir->make<TH1F>("fsize"      , Form("run %06i: fragment size"    ,RunNumber), 1000,   0., 100000.);
  Hist->error           = Dir->make<TH1F>("error"      , Form("run %06i: error code"       ,RunNumber),  512,   0.,    512.);
  Hist->error_rate      = Dir->make<TH1F>("erate"      , Form("run %06i: error rate"       ,RunNumber),  100,   0.,    100.);
  Hist->n_nb_errors     = Dir->make<TH1F>("n_nb_err"   , Form("run %06i: N(nb errors)"     ,RunNumber),  500,   0.,    500.);
  Hist->n_nwfs_errors   = Dir->make<TH1F>("n_nwfs_err" , Form("run %06i: N(nwfs errors)"   ,RunNumber),  500,   0.,    500.);
  Hist->n_linkid_errors = Dir->make<TH1F>("n_lnid_err" , Form("run %06i: N(link id errors)",RunNumber),  500,   0.,    500.);
  Hist->n_chid_errors   = Dir->make<TH1F>("n_chid_err" , Form("run %06i: N(chid errors)"   ,RunNumber),  500,   0.,    500.);
  Hist->n_nchh_errors   = Dir->make<TH1F>("n_nchh_err" , Form("run %06i: N(chh errors)"    ,RunNumber),  500,   0.,    500.);

  Hist->valid           = Dir->make<TH1F>("valid"      , Form("run %06i: valid code"       ,RunNumber),  100,   0.,    100.);
}

//-----------------------------------------------------------------------------
void TrackerDQM::book_histograms(int RunNumber) {
  art::ServiceHandle<art::TFileService> tfs;
    
  TLOG(TLVL_INFO) << "starting";

  art::TFileDirectory top_dir = tfs->mkdir("trk");

  book_event_histograms(&top_dir,RunNumber,&_hist.event);

  int station = 0;
  int plane   = 0;

  for (int i=0; i<_nActiveLinks; i++) {
    int link  = _activeLinks[i];
    art::TFileDirectory roc_dir = top_dir.mkdir(Form("roc_%i",link));
    book_roc_histograms(&roc_dir,RunNumber,&_hist.roc[station][plane][link],link);
  }
  
  TLOG(TLVL_INFO) << Form("pointer to the module: 0x%8p\n",(void*) this);
}

//-----------------------------------------------------------------------------
void TrackerDQM::beginJob() {
  TLOG(TLVL_INFO) << "starting";

  //  int           tmp_argc(2);
  int           tmp_argc(0);
  char**        tmp_argv(nullptr);

  // tmp_argv    = new char*[2];
  // tmp_argv[0] = new char[100];
  // tmp_argv[1] = new char[100];

  // strcpy(tmp_argv[0],"-b");
  // strcpy(tmp_argv[1],Form("--web=server:%d",_port));

  _app = new TApplication("TrackerDQM", &tmp_argc, tmp_argv);
  gROOT->SetWebDisplay(Form("server:%d",_port));

  // app->Run()
  // app_->Run(true);
  delete [] tmp_argv;

}

//-----------------------------------------------------------------------------
void TrackerDQM::endJob() {
  delete _canvas[0];
  delete _canvas[1];
  delete _canvas[2];
}

//-----------------------------------------------------------------------------
void TrackerDQM::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();

  if (_initialized != 0) return;
  _initialized = 1;

  _canvas[0] = new TCanvas("canvas_000");
  _canvas[0]->Divide(2,2);

  _canvas[1] = new TCanvas("canvas_001");
  _canvas[1]->Divide(2,2);

  _canvas[2] = new TCanvas("canvas_002");
  _canvas[2]->Divide(2,2);

  _browser   = new TBrowser();

  book_histograms(rn);

  _canvas[0]->cd(1);
  _hist.event.nbtot->Draw();
  _canvas[0]->cd(2);
  _hist.event.error_rate->Draw();
  _canvas[0]->cd(3);
  _hist.event.nfrag->Draw();
  _canvas[0]->cd(4);
  _hist.event.fsize->Draw();

  _canvas[1]->cd(1);
  _hist.event.error->Draw();
  _canvas[1]->cd(2);
  _hist.event.nhits->Draw();
  _canvas[1]->cd(3);
  _hist.event.n_chid_errors->Draw();
  _canvas[1]->cd(4);
  _hist.event.n_nchh_errors->Draw();
//-----------------------------------------------------------------------------
// up to four waveforms in a given channel
//-----------------------------------------------------------------------------
  
  _plot_wf.link    = 0;
  _plot_wf.channel = 4;

  if (_plot_wf.channel >= 0) {
    int station = 0;
    int plane   = 0;

    for (int i=0; i<4; i++) {
      _canvas[2]->cd(i+1);
      _hist.roc[station][plane][_plot_wf.link].channel[_plot_wf.channel].wf[i]->Draw();
    }
  }
}

//-----------------------------------------------------------------------------
// in the end of the run need to save histograms - will figure that out
//-----------------------------------------------------------------------------
void TrackerDQM::endRun(const art::Run& aRun) {
  //  int rn  = aRun.run();

  _initialized = 0;
}

// //-----------------------------------------------------------------------------
//   void TrackerDQM::fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data) {
//     Hist->nhits->Fill(Data->nhits);
//   }

//-----------------------------------------------------------------------------
void TrackerDQM::unpack_adc_waveform(mu2e::TrackerDataDecoder::TrackerDataPacket* Hit, float* Wf, WfParam_t* Wp) {

    int n_adc_packets = Hit->NumADCPackets;

    Wf[ 0] = reverseBits(Hit->ADC00);
    Wf[ 1] = reverseBits(Hit->ADC01A + (Hit->ADC01B << 6));
    Wf[ 2] = reverseBits(Hit->ADC02);

    for (int i=0; i<n_adc_packets; i++) {
      mu2e::TrackerDataDecoder::TrackerADCPacket* ahit = (mu2e::TrackerDataDecoder::TrackerADCPacket*) (((uint16_t*) Hit)+8+8*i);
      int loc = 12*i+2;

      Wf[loc+ 1] = reverseBits(ahit->ADC0);
      Wf[loc+ 2] = reverseBits(ahit->ADC1A + (ahit->ADC1B << 6));
      Wf[loc+ 3] = reverseBits(ahit->ADC2);
      Wf[loc+ 4] = reverseBits(ahit->ADC3);
      Wf[loc+ 5] = reverseBits(ahit->ADC4A + (ahit->ADC4B << 6));
      Wf[loc+ 6] = reverseBits(ahit->ADC5);
      Wf[loc+ 7] = reverseBits(ahit->ADC6);
      Wf[loc+ 8] = reverseBits(ahit->ADC7A + (ahit->ADC7B << 6));
      Wf[loc+ 9] = reverseBits(ahit->ADC8);
      Wf[loc+10] = reverseBits(ahit->ADC9);
      Wf[loc+11] = reverseBits(ahit->ADC10A + (ahit->ADC10B << 6));
      Wf[loc+12] = reverseBits(ahit->ADC11);
    }
//-----------------------------------------------------------------------------
// waveform processing
// 1. determine the baseline
//-----------------------------------------------------------------------------
    Wp->bl = 0;
    for (int i=0; i<_nSamplesBL; i++) {
      Wp->bl += Wf[i];
    }
    Wp->bl = Wp->bl/_nSamplesBL;
//-----------------------------------------------------------------------------
// 2. subtract the baseline and calculate the charge
//-----------------------------------------------------------------------------
    int nsamples = 15+12*(Hit->NumADCPackets-1);
    for (int i=0; i<nsamples; i++) {
      Wf[i] = Wf[i]-Wp->bl;
    }

    int tail  =  0;
    Wp->fs    = -1;
    Wp->q     =  0;
    Wp->qt    =  0;
    Wp->ph    = -1;
    Wp->q_x_i =  0;
    Wp->ns    =  0;

    for (int i=_nSamplesBL; i<nsamples; i++) {
      if (Wf[i] > _minPulseHeight) {
        if (tail == 0) {
                                        // first sample above the threshold
          if (Wp->fs < 0) Wp->fs = i;

          Wp->q += Wf[i];
                                        // as integrating the charge, also calculate the time
          Wp->q_x_i += Wf[i]*i;
          Wp->ns    += 1;

          if (Wf[i] > Wp->ph) {
            Wp->ph = Wf[i];
          }
        }
      }
      else if (Wf[i] < 0) {
        if (Wp->ph > 0) {
          tail  = 1;
        }
        if (tail == 1) Wp->qt -= Wf[i];
      }
    }
    Wp->tm = Wp->q_x_i/(Wp->q+1.e-12);
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    if (Wp->q < 100) {
      TLOG(TLVL_DEBUG+10) << "event=" << _event->run() << ":" << _event->subRun() << ":" << _event->event() 
                          << " Q=" << Wp->q;
    }
  }

//-----------------------------------------------------------------------------
  void TrackerDQM::fill_roc_histograms(RocHist_t* Hist, RocData_t* Rd) {
    
    Hist->nbytes->Fill  (Rd->nbytes);
    Hist->npackets->Fill(Rd->npackets);
    Hist->nhits->Fill   (Rd->nhits);
    Hist->valid->Fill   (Rd->valid);

    Hist->dt0r01->Fill  (Rd->dt0r01);
    Hist->dt1r01->Fill  (Rd->dt1r01);

//-----------------------------------------------------------------------------
// fill channel histograms
//-----------------------------------------------------------------------------
    for (int ich=0; ich<kNChannels; ich++) {
      ChannelData_t* chd = &Rd->channel[ich];
      ChannelHist_t* hch = &Hist->channel[ich];

      int fpga = _adc_index_1[ich] / 48;

      hch->nhits->Fill(chd->nhits);
      hch->dt0r->Fill(chd->dt0r);
      hch->dt1r->Fill(chd->dt1r);

      if (chd->error > 0) {
//-----------------------------------------------------------------------------
// channel error = 0 or 1(too many hits), count the number of "channel overflows"
//-----------------------------------------------------------------------------
        Hist->sum_error_vs_ch->Fill(ich);
      }
            
      Hist->dt0r_vs_ch->Fill(ich,chd->dt0r);
      Hist->dt1r_vs_ch->Fill(ich,chd->dt1r);

      Hist->dt0rc_vs_ch[fpga]->Fill(ich,chd->dt0r_c);
      Hist->dt1rc_vs_ch[fpga]->Fill(ich,chd->dt1r_c);
            
      int iadc = _adc_index_1[ich];
      Hist->dt0rc_vs_adc[fpga]->Fill(iadc,chd->dt0r_c);
      Hist->dt1rc_vs_adc[fpga]->Fill(iadc,chd->dt1r_c);
//-----------------------------------------------------------------------------
// there are hits in this channel
//-----------------------------------------------------------------------------
      for (int ih=0; ih<chd->nhits; ih++) {
        mu2e::TrackerDataDecoder::TrackerDataPacket* hit = chd->hit[ih];

        uint32_t corr_tdc0 = correctedTDC(hit->TDC0());
        uint32_t corr_tdc1 = correctedTDC(hit->TDC1());

        hch->time[0]->Fill(corr_tdc0*_tdc_bin); // in us
        hch->time[1]->Fill(corr_tdc1*_tdc_bin); // in us

        hch->t0  [0]->Fill(corr_tdc0*_tdc_bin_ns);
        hch->t0  [1]->Fill(corr_tdc1*_tdc_bin_ns);

        hch->t1  [0]->Fill(_timeWindow-corr_tdc0*_tdc_bin_ns);
        hch->t1  [1]->Fill(_timeWindow-corr_tdc1*_tdc_bin_ns);

        hch->tot [0]->Fill(hit->TOT0);
        hch->tot [1]->Fill(hit->TOT1);
        hch->pmp    ->Fill(hit->PMP);
//-----------------------------------------------------------------------------
// waveforms in a given channel
// one more sanity check: if hit reports too many samples, just skip it
//-----------------------------------------------------------------------------
        int nsamples = 15+12*(hit->NumADCPackets-1);
        if (nsamples > kMaxNSamples) {
          _event_data.error         |= kNWfsErrorBit;
          _event_data.n_nwfs_errors += 1;

          TLOG(TLVL_DEBUG+8) << "005 event " << _event->run() << ":" << _event->subRun() << ":" << _event->event()
                             << " too many samples:" << nsamples;
        }
        else {
          float     wform[kMaxNSamples];

          WfParam_t* wpar = &chd->wp[ih];
          unpack_adc_waveform(hit,wform,wpar);

          hch->wf[ih]->Reset();
          for (int is=0; is<nsamples; is++) {
            hch->wf[ih]->Fill(is,wform[is]);
            hch->raw_wf[ih]->Fill(is,wform[is]+wpar->bl);
          }
                                        // also set bin errors to zero
          int nb =  hch->wf[ih]->GetNbinsX();
          for (int ib=0; ib<nb; ib++) {
            hch->wf[ih]->SetBinError(ib+1,0);
            hch->wf[ih]->SetOption("HIST");
          }
//-----------------------------------------------------------------------------
// reconstructed waveform parameters
//-----------------------------------------------------------------------------
          hch->fsample->Fill(wpar->fs);
          hch->tmean->Fill(wpar->tm);
          hch->bline->Fill(wpar->bl);
          hch->pheight->Fill(wpar->ph);
          hch->q->Fill(wpar->q);
          hch->qt->Fill(wpar->qt);
          hch->qtq->Fill(wpar->qt/(wpar->q+1e-12));

          Hist->fs_vs_ich->Fill(ich,wpar->fs);
          Hist->bl_vs_ich->Fill(ich,wpar->bl);
          Hist->ph_vs_ich->Fill(ich,wpar->ph);
          Hist->q_vs_ich->Fill(ich,wpar->q);
          Hist->qt_vs_ich->Fill(ich,wpar->qt);
          Hist->qtq_vs_ich->Fill(ich,wpar->qt/(wpar->q+1e-12));
        }
      }
//-----------------------------------------------------------------------------
// time distance between the two sequential hits - need at least two
//-----------------------------------------------------------------------------
      for (int ih=1; ih<chd->nhits; ih++) {
        int corr_tdc0_ih  = (int) correctedTDC(chd->hit[ih  ]->TDC0());
        int corr_tdc1_ih  = (int) correctedTDC(chd->hit[ih  ]->TDC1());
        int corr_tdc0_ih1 = (int) correctedTDC(chd->hit[ih-1]->TDC0());
        int corr_tdc1_ih1 = (int) correctedTDC(chd->hit[ih-1]->TDC1());

        double dt0        = (corr_tdc0_ih-corr_tdc0_ih1)*_tdc_bin;
        double dt1        = (corr_tdc1_ih-corr_tdc1_ih1)*_tdc_bin;
        double dt2        = (dt0+dt1)/2;

        Hist->channel[ich].dt0->Fill(dt0);
        Hist->channel[ich].dt1->Fill(dt1);
        Hist->channel[ich].dt2->Fill(dt2);
      }
    }
//-----------------------------------------------------------------------------
// make sure the number of hits doesn't exceed the maximum - although it shouldn't 
//-----------------------------------------------------------------------------
    for (int ich=0; ich<kNChannels; ich++) {
      int ind_0 = _adc_index_0[ich];
      int ind_1 = _adc_index_1[ich];

      int nh    = Rd->channel[ich].nhits;
                                        // number of hits in a channel vs the channel number
      Hist->nh_vs_ch->Fill(ich,nh);
      Hist->nh_vs_adc1->Fill(ind_1,nh);

      for (int ihit=0; ihit<nh; ihit++) {
        Hist->nhits_vs_ich->Fill(ich);
        Hist->nhits_vs_adc[0]->Fill(ind_0);
        Hist->nhits_vs_adc[1]->Fill(ind_1);
      }
    }

  }

//-----------------------------------------------------------------------------
void TrackerDQM::fill_event_histograms(EventHist_t* Hist, EventData_t* Data) {

  Hist->nbtot->Fill(Data->nbtot);
  Hist->nhits->Fill(Data->nhtot);
  Hist->nfrag->Fill(Data->nfrag);

  Hist->n_nb_errors->Fill(Data->n_nb_errors);
  Hist->n_nwfs_errors->Fill(Data->n_nwfs_errors);
  Hist->n_linkid_errors->Fill(Data->n_linkid_errors);
  Hist->n_chid_errors->Fill(Data->n_chid_errors);
  Hist->n_nchh_errors->Fill(Data->n_nchh_errors);

  for (int i=0; i<kNErrorBits; i++) {
    if ((Data->error >> i) & 0x1) {
      Hist->error_rate->Fill(i);
    }
  }

  for (int i=0; i<Data->nfrag; i++) {
    int fsize = Data->fragments[i].nbytes;
    Hist->fsize->Fill(fsize);
  }
}

//-----------------------------------------------------------------------------
// fill_roc_histograms also fills the channel histograms
// if in error, only histogram the error code
//-----------------------------------------------------------------------------
int TrackerDQM::fill_histograms() {

  _hist.event.error->Fill(_event_data.error);
  _hist.event.valid->Fill(_event_data.valid);

  fill_event_histograms(&_hist.event,&_event_data);

  if ((_event_data.error & kNBytesErrorBit) != 0) return -1;
//-----------------------------------------------------------------------------
// event error = 0x0001 at this point means that either the fragment data cant be unpacked
// or there is nothing to unpack
//-----------------------------------------------------------------------------
  for (int ist=0; ist<1; ist++) {
    for (int plane=0; plane<2; plane++) {
      for (int ir=0; ir<_nActiveLinks; ir++) {
        int link = _activeLinks[plane][ir];
        fill_roc_histograms(&_hist.roc[_station][plane][link],&_event_data.rdata[_station][plane][link]);
      }
    }
  }

  return 0;
}
//-----------------------------------------------------------------------------
// a fragment may have multiple ROC blocks
//-----------------------------------------------------------------------------
  void TrackerDQM::analyze_fragment(const art::Event& Evt, const artdaq::Fragment* Fragment) {

    short* fdata = (short*) Fragment->dataBegin();

    _event_data.fragments.push_back(FragmentData_t());
    FragmentData_t* fdt = &_event_data.fragments.back();
//-----------------------------------------------------------------------------
// fragment size is specified in longs and includes service data, don't use
//-----------------------------------------------------------------------------
    fdt->nbytes  = fdata[0];
    if (fdata[0] > _maxFragmentSize) {
      _event_data.error = (_event_data.error | kNBytesErrorBit);

      TLOG(TLVL_DEBUG+8) << "001 event " << Evt.run() << ":" << Evt.subRun() << ":" << Evt.event()
                         << " error:" << _event_data.error << " fdt->nbytes:" << fdt->nbytes << " BAIL OUT";
      return;
    }
//-----------------------------------------------------------------------------
// start handling the ROC data
//-----------------------------------------------------------------------------
    short* first_address = fdata+_dataHeaderOffset; // offset is specified in 2-byte words
    short* last_address  = fdata+fdt->nbytes/2; // 

    while (first_address < last_address) {

      DtcDataHeaderPacket_t* dh = (DtcDataHeaderPacket_t*) first_address;
      int link      = dh->ROCID;
//-----------------------------------------------------------------------------
// check link number
//-----------------------------------------------------------------------------
      int found = 0;
      for (int i=0; i<_nActiveLinks; i++) {
        if (_activeLinks[i] == link) {
          found = 1;
          break;
        }
      }

      if (found == 0) {
//-----------------------------------------------------------------------------
// channel ID not found, assume a bit error, print a warning, and proceed
//-----------------------------------------------------------------------------
        _event_data.error           |= kLinkIDErrorBit;
        _event_data.n_linkid_errors += 1;

        TLOG(TLVL_DEBUG+8) << "002 event " << Evt.run() << ":" << Evt.subRun() << ":" << Evt.event()
                           << " error:" << _event_data.error << " link:" << link << " BAIL OUT";
        return;
      }

      RocData_t* rd = &_event_data.rdata[_station][_plane][link];
//-----------------------------------------------------------------------------
// for a given FPGA, a reference channel is the first channel in the readout order
//-----------------------------------------------------------------------------
      ChannelData_t* ref_ch[2];

      ref_ch[0]     = &rd->channel[_referenceChannel[link][0]];
      ref_ch[1]     = &rd->channel[_referenceChannel[link][1]];
        
      rd->nbytes    = dh->byteCount;
      rd->npackets  = dh->nPackets;

//-----------------------------------------------------------------------------
// for now, assume that all hits in the run have the same number of packets per hit
// take that from the first hit
//-----------------------------------------------------------------------------
      mu2e::TrackerDataDecoder::TrackerDataPacket* hit0 ;
      hit0     = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (fdata+_dataHeaderOffset+0x08);
      int n_adc_packets = hit0->NumADCPackets;

      rd->nhits     = dh->nPackets/(n_adc_packets+1);         //  printf("nhits : %3i\n",nhits);
      rd->valid     = dh->valid;
      rd->dt0r01    = -1.e12;
      rd->dt1r01    = -1.e12;
      
      _event_data.nhtot += rd->nhits;
      _event_data.valid += dh->valid*10;
//-----------------------------------------------------------------------------
// initialize channel data
//-----------------------------------------------------------------------------
      for (int i=0; i<kNChannels; i++) {
        rd->channel[i].nhits = 0;
        rd->channel[i].error = 0;
      }

      for (int ihit=0; ihit<rd->nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
        mu2e::TrackerDataDecoder::TrackerDataPacket* hit ;
        int offset = ihit*(8+8*n_adc_packets);
        hit     = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (fdata+offset+_dataHeaderOffset+0x08);
        int ich = hit->StrawIndex;
//------------------------------------------------------------------------------
// handle known channel ID errors:
// - channel with wrong cID adds 1 to teh event error code
// - overflown (N(hits) > N(max)) channel adds 100 to the event error code
// - event with negative error code never get here
//-----------------------------------------------------------------------------
        if (ich > 0x80) ich = ich-0x80;

        if (ich > 95) {
          _event_data.error         |= kChIDErrorBit;
          _event_data.n_chid_errors += 1;

          TLOG(TLVL_DEBUG+8) << "003 event " << Evt.run() << ":" << Evt.subRun() << ":" << Evt.event()
                             << " error:" << _event_data.error << " link:" << link << " ihit:" << ihit 
                             << " ich:0x" << std::hex << hit->StrawIndex;
        }
        else {

          ChannelData_t* chd = &rd->channel[ich];

          int nh = chd->nhits;
          if (nh >= kMaxNHitsPerChannel) {
            if (chd->error == 0) {
              chd->error                 = 1;
              _event_data.error         |= kNChHitsErrorBit;
              _event_data.n_nchh_errors += 1;
            }

            TLOG(TLVL_DEBUG+8) << "004 event " << Evt.run() << ":" << Evt.subRun() << ":" << Evt.event()
                               << " error:" << _event_data.error << " link:" << link 
                               << " ich:0x" << std::hex << hit->StrawIndex 
                               << " N(hits) > " << kMaxNHitsPerChannel;
          }
          else {
            chd->hit[nh]   = hit;
            chd->nhits    += 1;
          }
        }
      }
//-----------------------------------------------------------------------------
// hits in all channels counted
// time difference between a channel and a reference channel
//-----------------------------------------------------------------------------
      for (int i=0; i<kNChannels; i++) {
        ChannelData_t* chd = &rd->channel[i];

        int nh   = chd->nhits;
        int fpga = _adc_index_1[i] / 48;

        ChannelData_t* rch = ref_ch[fpga];
//-----------------------------------------------------------------------------
// in most cases, the number of hits in the reference channel should be greater 
// than the number of channels in any other channel of a given FPGA
//-----------------------------------------------------------------------------
        int iref = _referenceChannel[link][fpga];
        int nhr = rd->channel[iref].nhits;
        if ((nhr > 0) and (nh > 0)) {
//-----------------------------------------------------------------------------
// at least one hit in both reference and test channels
//-----------------------------------------------------------------------------
          int t0r = correctedTDC(rch->hit[0]->TDC0());
          int t1r = correctedTDC(rch->hit[0]->TDC1());
          int t0  = correctedTDC(chd->hit[0]->TDC0());
          int t1  = correctedTDC(chd->hit[0]->TDC1());
          
          float dt_over_2(_dt/2);
          
          chd->dt0r   = (t0-t0r)*_tdc_bin_ns;        // convert to ns  

          chd->dt0r_c = chd->dt0r;
          if (chd->dt0r >  dt_over_2/2) chd->dt0r_c = chd->dt0r + _gen_offset[i] - _dt;
          if (chd->dt0r < -dt_over_2/2) chd->dt0r_c = chd->dt0r + _gen_offset[i];
          
          chd->dt1r   = (t1-t1r)*_tdc_bin_ns;        // convert to ns

          chd->dt1r_c = chd->dt1r;
          if (chd->dt1r >  dt_over_2/2) chd->dt1r_c = chd->dt1r + _gen_offset[i] - _dt;
          if (chd->dt1r < -dt_over_2/2) chd->dt1r_c = chd->dt1r + _gen_offset[i];
        }
      }
//-----------------------------------------------------------------------------
// time offset between the two pulsers for the same ROC
//-----------------------------------------------------------------------------
      if ((rd->ref_ch[0]->nhits > 0) and (rd->ref_ch[1]->nhits > 0)) {
        int t0r0   = correctedTDC(rd->ref_ch[0]->hit[0]->TDC0());
        int t1r0   = correctedTDC(rd->ref_ch[0]->hit[0]->TDC1());
        
        int t0r1   = correctedTDC(rd->ref_ch[1]->hit[0]->TDC0());
        int t1r1   = correctedTDC(rd->ref_ch[1]->hit[0]->TDC1());
        
        rd->dt0r01 = (t0r0-t0r1)*_tdc_bin_ns;        // convert to ns  
        rd->dt1r01 = (t1r0-t1r1)*_tdc_bin_ns;        // convert to ns  
      }
//-----------------------------------------------------------------------------
// address in 2-byte words (N(data packets)+data header packet)
//-----------------------------------------------------------------------------
      first_address += (dh->nPackets + 1)*8;
    }
  }

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TrackerDQM::analyze(const art::Event& event) {

  _event = &event;

  _event_data.nbtot = 0;
  _event_data.nhtot = 0;
  _event_data.nfrag = 0;
  _event_data.valid = 0;
  _event_data.error = 0;
//-----------------------------------------------------------------------------
// reset error counters
//-----------------------------------------------------------------------------
  _event_data.n_nb_errors   = 0;
  _event_data.n_nwfs_errors = 0;
  _event_data.n_chid_errors = 0;
  _event_data.n_nchh_errors = 0;

  _event_data.fragments.clear();
//-----------------------------------------------------------------------------
// first get all fragments, select Tracker ones
//-----------------------------------------------------------------------------
  auto handle = event.getValidHandle<std::vector<artdaq::Fragment> >("daq:TRK");
//-----------------------------------------------------------------------------
// calculate the fragment size manually - big thank you to designers (:
//----------------------------------------------------------------------------- 
  _ifrag = 0;

  for (const artdaq::Fragment& frag : *handle) {
    ushort* buf    = (ushort*) (frag.dataBegin());
    int     nbytes = buf[0];
    int     fsize  = frag.sizeBytes();

    if (nbytes < 2) {
      _event_data.error       |= kNBytesErrorBit;
      _event_data.n_nb_errors += 1;

      TLOG(TLVL_DEBUG+8) << Form("event %6i:%8i:%8i : ERROR:%i nbytes=%i",
                                 event.run(),event.subRun(),event.event(),
                                 _event_data.error,nbytes);
    }

    _event_data.nfrag += 1;
    _event_data.nbtot += nbytes;        // including artdaq part

    if (_diagLevel > 2) {
      TLOG(TLVL_DEBUG+10) << Form("---------- TRK fragment # %3i nbytes: %5i fsize: %5i error: %5i\n",
                                  _ifrag,nbytes,fsize,_event_data.error);
      printFragment(&frag,nbytes/2);
    }
    _plane = _ifrag; // a hack, FIXME
    if ((_event_data.error == 0) and _analyzeFragments) analyze_fragment(event,&frag);

    _ifrag++;
  }
//-----------------------------------------------------------------------------
// proxy for event histograms
//-----------------------------------------------------------------------------
  if (_diagLevel > 1) {
    if ((_event_data.nbtot >= _minNBytes) and (_event_data.nbtot <= _maxNBytes)) {
      TLOG(TLVL_DEBUG+10) << Form(" Run : %5i subrun: %5i event: %8i nfrag: %3i nbytes: %5i\n", 
                               event.run(),event.subRun(),event.event(), _event_data.nfrag, _event_data.nbtot);
    }
  }
//-----------------------------------------------------------------------------
// go into interactive mode, 
// fInteractiveMode = 0 : do not stop
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
// TModule::analyze(event);
//-----------------------------------------------------------------------------
// event data un(re)packed , fill histograms
//-----------------------------------------------------------------------------
  if (_analyzeFragments != 0) {
    fill_histograms();
  }
//-----------------------------------------------------------------------------
// DTC registers
//-----------------------------------------------------------------------------
  if (_dumpDTCRegisters) {
    auto h = event.getValidHandle<std::vector<artdaq::Fragment>>("daq:TRKDTC");

    for (const artdaq::Fragment& frag : *h) {
      int *buf  = (int*) (frag.dataBegin());
      int nreg  = buf[0];
      int fsize = frag.sizeBytes();
      printf("%s: -------- DTC registers dump n(reg)=%5i size: %5i\n",__func__,nreg,fsize);
      printFragment(&frag,2+4*nreg);
    }
  }

  gSystem->ProcessEvents();
//-----------------------------------------------------------------------------
// in the end, print a line per event with errors
// or print the same line for all events if in debug mode
//----------------------------------------------------------------------------- 
  if (_event_data.error != 0) {
    TLOG(TLVL_ERROR) << Form("event %6i:%8i:%8i : ERROR:%i",
                             event.run(),event.subRun(),event.event(),
                             _event_data.error);
  }
  else {
    TLOG(TLVL_DEBUG+10) << Form("event %6i:%8i:%8i : ERROR:%i",
                                event.run(),event.subRun(),event.event(),
                                _event_data.error);
  }
}


DEFINE_ART_MODULE(TrackerDQM)
