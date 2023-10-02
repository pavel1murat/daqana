
#include "daqana/mod/TrkFragmentAna_module.hh"

namespace mu2e {

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
// NWords : the number of short words
//-----------------------------------------------------------------------------
  void TrkFragmentAna::printFragment(const artdaq::Fragment* Frag, int NWords) {
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

// ======================================================================

  TrkFragmentAna::TrkFragmentAna(fhicl::ParameterSet const& PSet) : 
    THistModule      (PSet                  ,"TrkFragmentAna"),
    _diagLevel       (PSet.get<int>         ("diagLevel"       )), 
    _minNBytes       (PSet.get<int>         ("minNBytes"       )), 
    _maxNBytes       (PSet.get<int>         ("maxNBytes"       )), 
    _dataHeaderOffset(PSet.get<int>         ("dataHeaderOffset")), 
    _trkfCollTag     (PSet.get<std::string> ("trkfCollTag"     )),
    _dumpDTCRegisters(PSet.get<int>         ("dumpDTCRegisters")),
    _analyzeFragments(PSet.get<int>         ("analyzeFragments"))
  {
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

    _referenceChannel[0] = 91;
    _referenceChannel[1] = 94;

  }

//-----------------------------------------------------------------------------
void TrkFragmentAna::book_histograms(int RunNumber) {
  art::ServiceHandle<art::TFileService> tfs;

  art::TFileDirectory top_dir     = tfs->mkdir("trk");
  art::TFileDirectory frag_dir    = tfs->mkdir("trk/frag_0");

  _Hist.event.nhits           = top_dir.make<TH1F>("nhits"    , Form("run %06i: nhits total"  ,RunNumber), 1000, 0.,  1000.);
  _Hist.event.nhits_vs_ich    = top_dir.make<TH1F>("nh_vs_ich", Form("run %06i: nh vs ich"    ,RunNumber),  100, 0.,   100.);
  _Hist.event.nhits_vs_adc[0] = top_dir.make<TH1F>("nh_vs_adc_0", Form("run %06i: nh vs adc_0",RunNumber),  100, 0.,   100.);
  _Hist.event.nhits_vs_adc[1] = top_dir.make<TH1F>("nh_vs_adc_1", Form("run %06i: nh vs adc_1",RunNumber),  100, 0.,   100.);
  _Hist.event.nbtot           = top_dir.make<TH1F>("nbtot"    , Form("run %06i: nbytes total",RunNumber),10000, 0., 10000.);
  _Hist.event.nfrag           = top_dir.make<TH1F>("nfrag"    , Form("run %06i: n fragments" ,RunNumber),  100, 0.,   100.);
  _Hist.event.error           = top_dir.make<TH1F>("error"    , Form("run %06i: error code"  ,RunNumber),  100, 0.,   100.);

  _Hist.frag.nbytes           = frag_dir.make<TH1F>("nbytes"  , Form("run %06i: n bytes"     ,RunNumber),10000,    0., 10000.);
  _Hist.frag.dsize            = frag_dir.make<TH1F>("dsize"   , Form("run %06i: size()-nb"   ,RunNumber),  200, -100.,   100.);
  _Hist.frag.npackets         = frag_dir.make<TH1F>("npackets", Form("run %06i: n packets"   ,RunNumber), 1000,    0.,  1000.);
  _Hist.frag.nhits            = frag_dir.make<TH1F>("nhits"   , Form("run %06i: n hits"      ,RunNumber),  300,    0.,   300.);
  _Hist.frag.valid            = frag_dir.make<TH1F>("valid"   , Form("run %06i: valid"       ,RunNumber),    2,    0.,     2.);

  _Hist.frag.nh_vs_ch         = frag_dir.make<TH2F>("nh_vs_ch"  , Form("run %06i: nh vs ch"  ,RunNumber),  100,0.,100., 10,0,10);

  _Hist.frag.dt0r_vs_ch       = frag_dir.make<TH2F>("dt0r_vs_ch_0", Form("run %06i: dt0r vs ch[0]",RunNumber),  100,0.,100.,5000,-25000,25000);
  _Hist.frag.dt1r_vs_ch       = frag_dir.make<TH2F>("dt1r_vs_ch_0", Form("run %06i: dt1r vs ch[0]",RunNumber),  100,0.,100.,5000,-25000,25000);

  _Hist.frag.dt0r01           = frag_dir.make<TH1F>("dt0r01", Form("run %06i: dt0r01"             ,RunNumber), 40000,-20000,20000);
  _Hist.frag.dt1r01           = frag_dir.make<TH1F>("dt1r01", Form("run %06i: dt1r01"             ,RunNumber), 40000,-20000,20000);

  _Hist.frag.dt0rc_vs_ch[0]   = frag_dir.make<TH2F>("dt0rc_vs_ch_0", Form("run %06i: dt0rc vs ch[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  _Hist.frag.dt0rc_vs_ch[1]   = frag_dir.make<TH2F>("dt0rc_vs_ch_1", Form("run %06i: dt0rc vs ch[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  _Hist.frag.dt1rc_vs_ch[0]   = frag_dir.make<TH2F>("dt1rc_vs_ch_0", Form("run %06i: dt1rc vs ch[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  _Hist.frag.dt1rc_vs_ch[1]   = frag_dir.make<TH2F>("dt1rc_vs_ch_1", Form("run %06i: dt1rc vs ch[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  _Hist.frag.dt0rc_vs_adc[0]  = frag_dir.make<TH2F>("dt0rc_vs_adc_0", Form("run %06i: dt0rc vs adc[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  _Hist.frag.dt0rc_vs_adc[1]  = frag_dir.make<TH2F>("dt0rc_vs_adc_1", Form("run %06i: dt0rc vs adc[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  _Hist.frag.dt1rc_vs_adc[0]  = frag_dir.make<TH2F>("dt1rc_vs_adc_0", Form("run %06i: dt1rc vs adc[0], ns",RunNumber),  100,0.,100.,1000,-10,10);
  _Hist.frag.dt1rc_vs_adc[1]  = frag_dir.make<TH2F>("dt1rc_vs_adc_1", Form("run %06i: dt1rc vs adc[1], ns",RunNumber),  100,0.,100.,1000,-10,10);

  printf("[mu2e::TrkFragmentAna] pointer to the module: 0x%8p\n",(void*) this);

  for (int i=0; i<kNChannels; i++) {
    art::TFileDirectory chan_dir = tfs->mkdir(Form("trk/frag_0/ch_%02i",i));
    _Hist.channel[i].nhits   = chan_dir.make<TH1F>(Form("ch_%02i_nhits",i),Form("run %06i: ch %02i nhits"  ,RunNumber,i), 100, 0., 100.);
    _Hist.channel[i].time[0] = chan_dir.make<TH1F>(Form("ch_%02i_time0",i),Form("run %06i: ch %02i time[0]",RunNumber,i),1000, 0., 100.);
    _Hist.channel[i].time[1] = chan_dir.make<TH1F>(Form("ch_%02i_time1",i),Form("run %06i: ch %02i time[0]",RunNumber,i),1000, 0., 100.);
    _Hist.channel[i].tot [0] = chan_dir.make<TH1F>(Form("ch_%02i_tot0" ,i),Form("run %06i: ch %02i time[0]",RunNumber,i), 100, 0., 100.);
    _Hist.channel[i].tot [1] = chan_dir.make<TH1F>(Form("ch_%02i_tot1" ,i),Form("run %06i: ch %02i time[0]",RunNumber,i), 100, 0., 100.);
    _Hist.channel[i].pmp     = chan_dir.make<TH1F>(Form("ch_%02i_pmp"  ,i),Form("run %06i: ch %02i pmp"    ,RunNumber,i), 100, 0.,  10.);
    _Hist.channel[i].dt0     = chan_dir.make<TH1F>(Form("ch_%02i_dt0"  ,i),Form("run %06i: ch %02i T0(i+1)-T0(i)",RunNumber,i)      ,50000,  0.,50);
    _Hist.channel[i].dt1     = chan_dir.make<TH1F>(Form("ch_%02i_dt1"  ,i),Form("run %06i: ch %02i T1(i+1)-T1(i)",RunNumber,i)      ,50000,  0.,50);
    _Hist.channel[i].dt2     = chan_dir.make<TH1F>(Form("ch_%02i_dt2"  ,i),Form("run %06i: ch %02i T2(i+1)-T2(i)",RunNumber,i)      ,50000,  0.,50);
    _Hist.channel[i].dt0r    = chan_dir.make<TH1F>(Form("ch_%02i_dt0r_0" ,i),Form("run %06i: ch %02i T0(ich,0)-T0(ref,0)[0]",RunNumber,i),20000,-10.,10);
    _Hist.channel[i].dt1r    = chan_dir.make<TH1F>(Form("ch_%02i_dt1r_0" ,i),Form("run %06i: ch %02i T1(ich,0)-T1(ref,0)[0]",RunNumber,i),20000,-10.,10);

    for (int j=0; j<kMaxNHitsPerChannel; j++) {
      _Hist.channel[i].wf[j] = chan_dir.make<TH1F>(Form("h_wf_ch_%02i_%i",i,j),Form("run %06i: ch [%02i][%i] waveform",RunNumber,i,j),20, 0.,20.);
    }
  }

}

//-----------------------------------------------------------------------------
void TrkFragmentAna::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();

  if (rn <= 285) _dataHeaderOffset =  0;
  else           _dataHeaderOffset = 32;

  book_histograms(rn);
//-----------------------------------------------------------------------------
// init timing offsets for known runs
// 'offset' is defined from the analysis of the dt0r_vs_ch_0
//-----------------------------------------------------------------------------
  double f0(31.29e6);                   // 31.29 MHz
  double offset;

  if      (rn ==    281) { _freq  = f0/(pow(2,9)+1); offset = 15030; } // 60 kHz
  else if (rn == 105023) { _freq  = f0/(pow(2,7)+1); offset =     0; } // 250 kHz, undefined - no hits in lanes 2 and 3 (HV)
  else if (rn == 105026) { _freq  = f0/(pow(2,9)+1); offset = 10400; } 
  else if (rn == 105038) { _freq  = f0/(pow(2,9)+1); offset =  7964; } 
  else if (rn == 105041) { _freq  = f0/(pow(2,9)+1); offset =  1128; } //  60 kHz
  else if (rn == 105042) { _freq  = f0/(pow(2,9)+1); offset =  1128; } 
  else if (rn == 105043) { _freq  = f0/(pow(2,9)+1); offset =  1128; }
  else if (rn == 105044) { _freq  = f0/(pow(2,9)+1); offset =  1129; } // 60 kHz
  else if (rn == 105060) { _freq  = f0/(pow(2,9)+1); offset = 11130; } // 60 kHz,1000x25 usec
  else if (rn == 105066) { _freq  = f0/(pow(2,9)+1); offset = 10580; } // 60 kHz, 700x25 usec

                                        // in nanoseconds
  _dt   = 1/_freq*1.e9;
//-----------------------------------------------------------------------------
// add 20 ns to have the HV lanes distinct on the plot
//-----------------------------------------------------------------------------
  for (int i=0; i<kNChannels; i++) {
    if (_adc_index_0[i] < 48) _gen_offset[i] = 0;
    else                      _gen_offset[i] = offset;
  }

}

//--------------------------------------------------------------------------------
void TrkFragmentAna::beginJob() {
}

//-----------------------------------------------------------------------------
void TrkFragmentAna::endJob() {
  printf("[mu2e::TrkFragmentAna] pointer to the module: 0x%8p\n",(void*) this);
}

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TrkFragmentAna::analyze(const art::Event& event) {
  //art::EventNumber_t eventNumber = event.event();

  int    nbtot = 0;
  int    nfrag = 0;

  for (int i=0; i<100; i++) _nwf[i] = 0;

  auto handle = event.getValidHandle<std::vector<artdaq::Fragment> >(_trkfCollTag);
//-----------------------------------------------------------------------------
// calculate the fragment size manually - big thank you to designers (:
//----------------------------------------------------------------------------- 
  _error = 0;

  int ifrag = 0;
  for (/* auto */ const artdaq::Fragment& frag : *handle) {
    ushort* buf = (ushort*) (frag.dataBegin());
    int nbytes  = buf[0];
    int fsize   = frag.sizeBytes();

    nfrag      += 1;
    nbtot      += nbytes;

    if (_analyzeFragments) analyze_fragment(&frag,&_Hist.frag);

    if (_diagLevel > 2) {
      printf("%s: ---------- TRK fragment # %3i nbytes: %5i fsize: %5i\n",__func__,ifrag,nbytes,fsize);
      printFragment(&frag,nbytes/2);
    }
    ifrag++;
  }

  if ((_error > 0) and (_error < 10)) {
    printf(" Run : %5i subrun: %5i event: %8i error: %5i\n", 
           event.run(),event.subRun(),event.event(), _error);
  }
//-----------------------------------------------------------------------------
// proxy for event histograms
//-----------------------------------------------------------------------------
  _Hist.event.nbtot->Fill(nbtot);
  _Hist.event.nfrag->Fill(nfrag);
  _Hist.event.error->Fill(_error);

  if (_diagLevel > 1) {
    if ((nbtot >= _minNBytes) and (nbtot <= _maxNBytes)) {
      printf(" Run : %5i subrun: %5i event: %8i nfrag: %3i nbytes: %5i\n", 
	     event.run(),event.subRun(),event.event(), nfrag, nbtot);
    }
  }
//-----------------------------------------------------------------------------
// go into interactive mode, 
// fInteractiveMode = 0 : do not stop
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
  TModule::analyze(event);
//-----------------------------------------------------------------------------
// DTC registers
//-----------------------------------------------------------------------------
  if (_dumpDTCRegisters) {
    auto h = event.getValidHandle<std::vector<artdaq::Fragment>>("daq:TRKDTC");

    for (/* auto */ const artdaq::Fragment& frag : *h) {
      int *buf  = (int*) (frag.dataBegin());
      int nreg  = buf[0];
      int fsize = frag.sizeBytes();
      printf("%s: -------- DTC registers dump n(reg)=%5i size: %5i\n",__func__,nreg,fsize);
      printFragment(&frag,2+4*nreg);
    }
  }
  
}

//-----------------------------------------------------------------------------
// void TrkFragmentAna::analyze_tracker(const mu2e::TrackerFragment& Fragment) {
  void TrkFragmentAna::analyze_fragment(const artdaq::Fragment* Fragment, FragmentHist_t* Hist) {

    double tdc_bin(5/256.*1e-3); // TDC bin width (Richie), in us

    int size   = Fragment->size();

    short* data = (short*) Fragment->dataBegin();

    int nbytes = data[0];

    DtcDataHeaderPacket_t* dh = (DtcDataHeaderPacket_t*) (data + _dataHeaderOffset);

    int npackets = dh->nPackets;
    int nhits    = npackets/2;

    //  printf("nhits : %3i\n",nhits);

    Hist->nbytes->Fill(nbytes);
    Hist->dsize->Fill(size-nbytes);
    Hist->npackets->Fill(npackets);
    Hist->nhits->Fill(nhits);
    Hist->valid->Fill(dh->valid);

    _error  += dh->valid*10;

    for (int i=0; i<kNChannels; i++) {
      _data.channel[i].nhits = 0;
    }

    for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
      TrackerFragment::DataPacket* hit ;
      hit     = (TrackerFragment::DataPacket*) (data+ihit*0x10+_dataHeaderOffset+0x08);
      int ich = hit->StrawIndex;

      if (ich > 128) ich = ich-128;

      if (ich > 95) {
        printf ("ERROR in %s: ich = %i, BAIL OUT\n",__func__,ich);
        _error = 1;
        return;
      }

      int nh = _data.channel[ich].nhits;

      if (nh >= kMaxNHitsPerChannel) {
        printf ("ERROR in %s: ich = %i, N(hits) >= %i BAIL OUT\n",__func__,ich,kMaxNHitsPerChannel);
        _error = 2;
        return;
      }

      _data.channel[ich].hit[nh] = hit;
      _data.channel[ich].nhits  += 1;

      uint32_t corr_tdc0 = correctedTDC(hit->TDC0());
      uint32_t corr_tdc1 = correctedTDC(hit->TDC1());

      _Hist.channel[ich].time[0]->Fill(corr_tdc0*tdc_bin);
      _Hist.channel[ich].time[1]->Fill(corr_tdc1*tdc_bin);

      _Hist.channel[ich].tot [0]->Fill(hit->TOT0);
      _Hist.channel[ich].tot [1]->Fill(hit->TOT1);
      _Hist.channel[ich].pmp    ->Fill(hit->PMP);
//-----------------------------------------------------------------------------
// waveforms in a given channel
//-----------------------------------------------------------------------------
      uint16_t adc[15];

      adc[ 0] = reverseBits(hit->_ADC00);
      adc[ 1] = reverseBits(hit->_ADC01A + (hit->_ADC01B << 6));
      adc[ 2] = reverseBits(hit->_ADC02);

      TrackerFragment::ADCPacket* ahit = (TrackerFragment::ADCPacket*) ((uint16_t*)hit+6);

      adc[ 3] = reverseBits(ahit->_ADC0);
      adc[ 4] = reverseBits(ahit->_ADC1A + (ahit->_ADC1B << 6));
      adc[ 5] = reverseBits(ahit->_ADC2);
      adc[ 6] = reverseBits(ahit->_ADC3);
      adc[ 7] = reverseBits(ahit->_ADC4A + (ahit->_ADC4B << 6));
      adc[ 8] = reverseBits(ahit->_ADC5);
      adc[ 9] = reverseBits(ahit->_ADC6);
      adc[10] = reverseBits(ahit->_ADC7A + (ahit->_ADC7B << 6));
      adc[11] = reverseBits(ahit->_ADC5);
      adc[12] = reverseBits(ahit->_ADC6);
      adc[13] = reverseBits(ahit->_ADC10A + (ahit->_ADC10B << 6));
      adc[14] = reverseBits(ahit->_ADC11);

      _Hist.channel[ich].wf[nh]->Reset();
      for (int is=0; is<15; is++) {
        _Hist.channel[ich].wf[nh]->Fill(is,adc[is]);
      }
    }
//-----------------------------------------------------------------------------
// time difference between  two sequential hits
//-----------------------------------------------------------------------------
    for (int i=0; i<kNChannels; i++) {
      int nh = _data.channel[i].nhits;
      _Hist.channel[i].nhits->Fill(nh);
      for (int ih=1; ih<nh; ih++) {
        int corr_tdc0_ih  = (int) correctedTDC(_data.channel[i].hit[ih  ]->TDC0());
        int corr_tdc1_ih  = (int) correctedTDC(_data.channel[i].hit[ih  ]->TDC1());
        int corr_tdc0_ih1 = (int) correctedTDC(_data.channel[i].hit[ih-1]->TDC0());
        int corr_tdc1_ih1 = (int) correctedTDC(_data.channel[i].hit[ih-1]->TDC1());

        float dt0         = (corr_tdc0_ih-corr_tdc0_ih1)*tdc_bin;
        float dt1         = (corr_tdc1_ih-corr_tdc1_ih1)*tdc_bin;
        float dt2         = (dt0+dt1)/2;

        _Hist.channel[i].dt0->Fill(dt0);
        _Hist.channel[i].dt1->Fill(dt1);
        _Hist.channel[i].dt2->Fill(dt2);
      }
    }
//-----------------------------------------------------------------------------
// time difference between a channel and a reference channel
//-----------------------------------------------------------------------------
    ChData_t* ref_ch[2];

    ref_ch[0] = &_data.channel[_referenceChannel[0]];
    ref_ch[1] = &_data.channel[_referenceChannel[1]];
    
    for (int i=0; i<kNChannels; i++) {
      ChData_t* ch = &_data.channel[i];

      int nh = ch->nhits;
      _Hist.channel[i].nhits->Fill(nh);

      int fpga = _adc_index_1[i] / 48;

      ChData_t* rch = ref_ch[fpga];
//-----------------------------------------------------------------------------
// two times corresponding to two ends of the straw 
//-----------------------------------------------------------------------------
      int t0r = correctedTDC(rch->hit[0]->TDC0());
      int t1r = correctedTDC(rch->hit[0]->TDC1());
//-----------------------------------------------------------------------------
// in most cases, the number of hits in the reference channel should be greater 
// than the number of channels in any other channel of a given FPGA
//-----------------------------------------------------------------------------
      int nhr = _data.channel[_referenceChannel[fpga]].nhits;
      if (nhr == 0) {
        if (nh > 0) { 
          _Hist.channel[i].dt0r->Fill(-1.e6);
          _Hist.channel[i].dt1r->Fill(-1.e6);
        }
      }
      else {
        if (nh > 0) {
//-----------------------------------------------------------------------------
// at least one hit in both reference and test channels
//-----------------------------------------------------------------------------
          int t0 = correctedTDC(ch->hit[0]->TDC0());
          int t1 = correctedTDC(ch->hit[0]->TDC1());
          
          float dt_over_2(_dt/2);
          
          float dt0r   = (t0-t0r)*tdc_bin*1.e3;        // convert to ns  

          float dt0r_c = dt0r;
          if (dt0r >  dt_over_2/2) dt0r_c = dt0r + _gen_offset[i] - _dt;
          if (dt0r < -dt_over_2/2) dt0r_c = dt0r + _gen_offset[i];
          
          float dt1r   = (t1-t1r)*tdc_bin*1.e3;        // convert to ns

          float dt1r_c = dt1r;
          if (dt1r >  dt_over_2/2) dt1r_c = dt1r + _gen_offset[i] - _dt;
          if (dt1r < -dt_over_2/2) dt1r_c = dt1r + _gen_offset[i];
          
          _Hist.channel[i].dt0r->Fill(dt0r);
          _Hist.channel[i].dt1r->Fill(dt1r);
          
          _Hist.frag.dt0r_vs_ch->Fill(i,dt0r);
          _Hist.frag.dt1r_vs_ch->Fill(i,dt1r);
          
          _Hist.frag.dt0rc_vs_ch[fpga]->Fill(i,dt0r_c);
          _Hist.frag.dt1rc_vs_ch[fpga]->Fill(i,dt1r_c);

          int iadc = _adc_index_1[i];
          _Hist.frag.dt0rc_vs_adc[fpga]->Fill(iadc,dt0r_c);
          _Hist.frag.dt1rc_vs_adc[fpga]->Fill(iadc,dt1r_c);
        }
      }
    }
//-----------------------------------------------------------------------------
// time offset between the two pulsers
//-----------------------------------------------------------------------------
    if ((ref_ch[0]->nhits > 0) and (ref_ch[1]->nhits > 0)) {
      int t0r0 = correctedTDC(ref_ch[0]->hit[0]->TDC0());
      int t1r0 = correctedTDC(ref_ch[0]->hit[0]->TDC1());

      int t0r1 = correctedTDC(ref_ch[1]->hit[0]->TDC0());
      int t1r1 = correctedTDC(ref_ch[1]->hit[0]->TDC1());

      float dt0r01 = (t0r0-t0r1)*tdc_bin*1.e3;        // convert to ns  
      float dt1r01 = (t1r0-t1r1)*tdc_bin*1.e3;        // convert to ns  

      _Hist.frag.dt0r01->Fill(dt0r01);
      _Hist.frag.dt1r01->Fill(dt1r01);
    }
//-----------------------------------------------------------------------------
// this part needs to be changed 
//-----------------------------------------------------------------------------
    _Hist.event.nhits->Fill(nhits);

    for (int ich=0; ich<kNChannels; ich++) {
      int ind_0 = _adc_index_0[ich];
      int ind_1 = _adc_index_1[ich];
      int nh    = _data.channel[ich].nhits;
      for (int ihit=0; ihit<nh; ihit++) {
        _Hist.event.nhits_vs_ich->Fill(ich);
        _Hist.event.nhits_vs_adc[0]->Fill(ind_0);
        _Hist.event.nhits_vs_adc[1]->Fill(ind_1);
      }
    }
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TrkFragmentAna)
