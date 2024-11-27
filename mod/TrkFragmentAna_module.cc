///////////////////////////////////////////////////////////////////////////////
// event error codes: 
// 0: OK
// 1:
// 2:
// 3: fragment too long
// 4: wrong active link ID
// 5: wrong N waveform packets
//
// use of debug bits :
// bit 00: print all errors
// bit 02: events with given error code
// bit 03: 
//         0x01: events with the total number of errors > _minNErrors
//         0x10: total event dump
// bit 04: print problematic hits
///////////////////////////////////////////////////////////////////////////////
#include "TRACE/tracemf.h"
#define TRACE_NAME "TrkFragmentAna"

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
  int TrkFragmentAna::unpack_adc_waveform(TrackerDataDecoder::TrackerDataPacket* Hit, float* Wf, WfParam_t* Wp) {

    Wf[ 0] = reverseBits(Hit->ADC00);
    Wf[ 1] = reverseBits(Hit->ADC01A + (Hit->ADC01B << 6));
    Wf[ 2] = reverseBits(Hit->ADC02);

    for (int i=0; i<_nADCPackets; i++) {
      TrackerDataDecoder::TrackerADCPacket* ahit = (TrackerDataDecoder::TrackerADCPacket*) (((uint16_t*) Hit)+8+8*i);
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

    int   tail  = 0;
    Wp->fs = -1;
    Wp->q  = 0;
    Wp->qt = 0;
    Wp->ph = -1;
    for (int i=_nSamplesBL; i<nsamples; i++) {
      if (Wf[i] > _minPulseHeight) {
        if (tail == 0) {
                                        // first sample above the threshold
          if (Wp->fs < 0) Wp->fs = i;

          Wp->q += Wf[i];
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
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    if (Wp->q < 100) {
      TLOG(TLVL_DEBUG+1) << "event=" << _edata._event->run() << ":" << _edata._event->subRun() << ":" << _edata._event->event() 
                         << " Q=" << Wp->q;
    }
    return 0;
  }
  
//-----------------------------------------------------------------------------
unsigned int correctedTDC(unsigned int TDC) {
  uint32_t corrected_tdc = ((TDC & 0xFFFF00) + (0xFF  - (TDC & 0xFF)));
  return corrected_tdc;
}

// ======================================================================

  TrkFragmentAna::TrkFragmentAna(fhicl::ParameterSet const& PSet) : 
    THistModule         (PSet                      ,"TrkFragmentAna"     ) ,

    _diagLevel          (PSet.get<int>             ("diagLevel"          )), 
    _minNBytes          (PSet.get<int>             ("minNBytes"          )), 
    _maxNBytes          (PSet.get<int>             ("maxNBytes"          )), 
    _dataHeaderOffset   (PSet.get<int>             ("dataHeaderOffset"   )),
    _activeLinks_0      (PSet.get<std::vector<int>>("activeLinks_0"      )),
    _activeLinks_1      (PSet.get<std::vector<int>>("activeLinks_1"      )),
    _refChCal           (PSet.get<std::vector<int>>("refChCal"           )),
    _refChHV            (PSet.get<std::vector<int>>("refChHV"            )),
    _trkfCollTag        (PSet.get<std::string>     ("trkfCollTag"        )),
    _dumpDTCRegisters   (PSet.get<int>             ("dumpDTCRegisters"   )),
    _analyzeFragments   (PSet.get<int>             ("analyzeFragments"   )),
    _maxFragmentSize    (PSet.get<int>             ("maxFragmentSize"    )),
    _pulserFrequency    (PSet.get<int>             ("pulserFrequency"    )),
    _nADCPackets        (PSet.get<int>             ("nADCPackets"        )),
    _nSamplesBL         (PSet.get<int>             ("nSamplesBL"         )),
    _minPulseHeight     (PSet.get<float>           ("minPulseHeight"     )),
    _minNErrors         (PSet.get<int>             ("minNErrors"         )),
    _errorCode          (PSet.get<int>             ("errorCode"          )),
    _validateAdcPatterns(PSet.get<int>             ("validateAdcPatterns"))
  {
    _activeLinks[0] = &_activeLinks_0;
    _activeLinks[1] = &_activeLinks_1;

    double f0(31.29e6);                   // 31.29 MHz

    _timeWindow = PSet.get<int>("timeWindow")*25.;  // in ns

    if      (_pulserFrequency ==  60) _freq = f0/(pow(2,9)+1);     // ~ 60 kHz
    else if (_pulserFrequency == 250) _freq = f0/(pow(2,7)+1);     // ~250 kHz

    _dt   = 1/_freq*1.e9;               // in ns
//------------------------------------------------------------------------------
// default map, Richie says TS1 may have an old firmware with some bugs
//-----------------------------------------------------------------------------
    int adc_index_0[96] = {
      91, 85, 79, 73, 67, 61, 55, 49,          // lane 0
      43, 37, 31, 25, 19, 13,  7,  1,
      90, 84, 78, 72, 66, 60, 54, 48,
      
      42, 36, 30, 24, 18, 12,  6,  0,          // lane 1
      93, 87, 81, 75, 69, 63, 57, 51,
      45, 39, 33, 27, 21, 15,  9,  3,
      
      44, 38, 32, 26, 20, 14,  8,  2,          // lane 2
      92, 86, 80, 74, 68, 62, 56, 50,
      47, 41, 35, 29, 23, 17, 11,  5,
      
      95, 89, 83, 77, 71, 65, 59, 53,          // lane 3
      46, 40, 34, 28, 22, 16, 10,  4,
      94, 88, 82, 76, 70, 64, 58, 52
    };
//-----------------------------------------------------------------------------
// TS1: compared to _0, swap lanes 3 and 4, and in the new lane3 swap lines 1 and 3
// the firmware was seemingly fixed since then
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
    _nActiveLinks[0]      = _activeLinks[0]->size();
    _nActiveLinks[1]      = _activeLinks[1]->size();

    for (int ilink=0; ilink<kMaxNLinks; ilink++) {
      _referenceChannel[ilink][0] = 91;
      _referenceChannel[ilink][1] = 94;
    }

    _nStations = 1; // for now
    for (int ist=0; ist<_nStations; ist++) {
      for (int idtc=0; idtc<2; idtc++) {
        for (int ilink=0; ilink<kMaxNLinks; ilink++) {
          RocData_t* rd = &_edata.station[ist].roc[idtc][ilink];
          rd->link      = ilink;

          if (ilink < _nActiveLinks[idtc]) {
            _referenceChannel[ilink][0] = _refChCal[ilink];
            _referenceChannel[ilink][1] = _refChHV [ilink];
          }

          rd->ref_ch[0] = &rd->channel[_referenceChannel[ilink][0]];
          rd->ref_ch[1] = &rd->channel[_referenceChannel[ilink][1]];
        }
      }
    }

    _tdc_bin             = (5/256.*1e-3);       // TDC bin width (Richie), in us
    _tdc_bin_ns          = _tdc_bin*1e3;        // convert to ns

    _initialized         = 0;
  }

//-----------------------------------------------------------------------------
// I : channel number
//-----------------------------------------------------------------------------
  void TrkFragmentAna::book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int Link, int I) {

    Hist->nhits   = Dir->make<TH1F>(Form("ch_%i_%02i_nhits" ,Link,I),Form("run %06i: link %i ch %02i nhits"  ,RunNumber,Link,I), 50, -0.5, 49.5);
    Hist->time[0] = Dir->make<TH1F>(Form("ch_%i_%02i_time0" ,Link,I),Form("run %06i: link %i ch %02i time[0]",RunNumber,Link,I),1000, 0., 100.);  // us
    Hist->time[1] = Dir->make<TH1F>(Form("ch_%i_%02i_time1" ,Link,I),Form("run %06i: link %i ch %02i time[0]",RunNumber,Link,I),1000, 0., 100.);  // us
    Hist->t0  [0] = Dir->make<TH1F>(Form("ch_%i_%02i_t0_0"  ,Link,I),Form("run %06i: link %i ch %02i t0[0]  ",RunNumber,Link,I),1000,-20., 80.);  // ns
    Hist->t0  [1] = Dir->make<TH1F>(Form("ch_%i_%02i_t0_1"  ,Link,I),Form("run %06i: link %i ch %02i t0[1]  ",RunNumber,Link,I),1000,-20., 80.);  // ns
    Hist->t1  [0] = Dir->make<TH1F>(Form("ch_%i_%02i_t1_0"  ,Link,I),Form("run %06i: link %i ch %02i t1[0]  ",RunNumber,Link,I),1000,-20., 80.);  // ns
    Hist->t1  [1] = Dir->make<TH1F>(Form("ch_%i_%02i_t1_1"  ,Link,I),Form("run %06i: link %i ch %02i t1[1]  ",RunNumber,Link,I),1000,-20., 80.);  // ns
    Hist->tot [0] = Dir->make<TH1F>(Form("ch_%i_%02i_tot0"  ,Link,I),Form("run %06i: link %i ch %02i tot[0]" ,RunNumber,Link,I), 100, 0., 100.);
    Hist->tot [1] = Dir->make<TH1F>(Form("ch_%i_%02i_tot1"  ,Link,I),Form("run %06i: link %i ch %02i tot[1]" ,RunNumber,Link,I), 100, 0., 100.);
    Hist->pmp     = Dir->make<TH1F>(Form("ch_%i_%02i_pmp"   ,Link,I),Form("run %06i: link %i ch %02i pmp"    ,RunNumber,Link,I), 100, 0.,  10.);
    Hist->dt01[0] = Dir->make<TH1F>(Form("ch_%i_%02i_dt01_0",Link,I),Form("run %06i: link %i ch %02i T0(i)-T1(i) [0],ns",RunNumber,Link,I) ,500, -25,25);
    Hist->dt01[1] = Dir->make<TH1F>(Form("ch_%i_%02i_dt01_1",Link,I),Form("run %06i: link %i ch %02i T0(i)-T1(i) [1],ns",RunNumber,Link,I) ,500, -2500,2500);
    Hist->dt0     = Dir->make<TH1F>(Form("ch_%i_%02i_dt0"   ,Link,I),Form("run %06i: link %i ch %02i T0(i+1)-T0(i)",RunNumber,Link,I)      ,50000,  0.,50);
    Hist->dt1     = Dir->make<TH1F>(Form("ch_%i_%02i_dt1"   ,Link,I),Form("run %06i: link %i ch %02i T1(i+1)-T1(i)",RunNumber,Link,I)      ,50000,  0.,50);
    Hist->dt2     = Dir->make<TH1F>(Form("ch_%i_%02i_dt2"   ,Link,I),Form("run %06i: link %i ch %02i T2(i+1)-T2(i)",RunNumber,Link,I)      ,50000,  0.,50);
    Hist->dt0r    = Dir->make<TH1F>(Form("ch_%i_%02i_dt0r_0",Link,I),Form("run %06i: link %i ch %02i T0(ich,0)-T0(ref,0)[0]",RunNumber,Link,I),20000,-10.,10);
    Hist->dt1r    = Dir->make<TH1F>(Form("ch_%i_%02i_dt1r_0",Link,I),Form("run %06i: link %i ch %02i T1(ich,0)-T1(ref,0)[0]",RunNumber,Link,I),20000,-10.,10);
//-----------------------------------------------------------------------------
// waveform parameters
//-----------------------------------------------------------------------------
    Hist->fsample = Dir->make<TH1F>(Form("ch_%i_%02i_fs"    ,Link,I),Form("run %06i: link %i ch %02i first sample"   ,RunNumber,Link,I), 30,-0.5, 29.5);
    Hist->bline   = Dir->make<TH1F>(Form("ch_%i_%02i_bl"    ,Link,I),Form("run %06i: link %i ch %02i WF baseline"    ,RunNumber,Link,I),250,0,500);
    Hist->pheight = Dir->make<TH1F>(Form("ch_%i_%02i_ph"    ,Link,I),Form("run %06i: link %i ch %02i WF pulse height",RunNumber,Link,I),500,0,500);
    Hist->q       = Dir->make<TH1F>(Form("ch_%i_%02i_q"     ,Link,I),Form("run %06i: link %i ch %02i WF charge"      ,RunNumber,Link,I),500,0,500);//
    Hist->qt      = Dir->make<TH1F>(Form("ch_%i_%02i_qt"    ,Link,I),Form("run %06i: link %i ch %02i WF tail charge" ,RunNumber,Link,I),500,0,500);
    Hist->qtq     = Dir->make<TH1F>(Form("ch_%i_%02i_qtq"   ,Link,I),Form("run %06i: link %i ch %02i WF Qt/Q"        ,RunNumber,Link,I),200,0,1);
//-----------------------------------------------------------------------------
// waveform histograms, assume number of samples < 30
//-----------------------------------------------------------------------------
    for (int j=0; j<kMaxNHWfPerChannel; j++) {
      Hist->raw_wf[j] = Dir->make<TH1F>(Form("raw_wf_ch_%i_%02i_%i",Link,I,j),Form("run %06i: link %i ch [%02i][%i] raw_waveform",RunNumber,Link,I,j),30, 0.,30.);
      Hist->wf    [j] = Dir->make<TH1F>(Form("wf_ch_%i_%02i_%i"    ,Link,I,j),Form("run %06i: link %i ch [%02i][%i] waveform"    ,RunNumber,Link,I,j),30, 0.,30.);
    }
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::book_dtc_histograms(art::TFileDirectory* Dir, int RunNumber, DtcHist_t* Hist, int IStation, int IDtc) {
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::book_roc_histograms(art::TFileDirectory* Dir, int RunNumber, RocHist_t* Hist, int Station, int Dtc, int Link) {
    Hist->nbytes          = Dir->make<TH1F>("nbytes"  ,        Form("run %06i: link_%02i:%i:%i n bytes"      ,RunNumber,Station,Dtc,Link),10000,    0., 10000.);
    Hist->npackets        = Dir->make<TH1F>("npackets",        Form("run %06i: link %02i:%i:%i n packets"    ,RunNumber,Station,Dtc,Link), 1000,    0.,  1000.);
    Hist->nhits           = Dir->make<TH1F>("nhits"   ,        Form("run %06i: link %02i:%i:%i n hits"       ,RunNumber,Station,Dtc,Link),  300,    0.,   300.);
    Hist->valid           = Dir->make<TH1F>("valid"   ,        Form("run %06i: link %02i:%i:%i valid"        ,RunNumber,Station,Dtc,Link),    2,    0.,     2.);
    Hist->error_code      = Dir->make<TH1F>("errcode" ,        Form("run %06i: link %02i:%i:%i errcode"      ,RunNumber,Station,Dtc,Link),  512,    0.,   512.);

    Hist->n_empty         = Dir->make<TH1F>("nempt"      , Form("run %06i: N(empty)"     ,RunNumber),  100, 0.,    100.);
    Hist->n_invalid_dr    = Dir->make<TH1F>("ninvr"      , Form("run %06i: N(invalid DR)",RunNumber),  100, 0.,    100.);
    Hist->n_corrupt       = Dir->make<TH1F>("ncorr"      , Form("run %06i: N(corrupt)"   ,RunNumber),  100, 0.,    100.);
    Hist->n_timeouts      = Dir->make<TH1F>("ntmo"       , Form("run %06i: N(timeouts)"  ,RunNumber),  100, 0.,    100.);
    Hist->n_overflows     = Dir->make<TH1F>("nover"      , Form("run %06i: N(overflows)" ,RunNumber),  100, 0.,    100.);

    Hist->nerr_tot        = Dir->make<TH1F>("nerr_tot"   , Form("run %06i: link %02i:%i:%i N errors"     ,RunNumber,Station,Dtc,Link), 1000, 0.,   2000.);
    Hist->nerr_vs_evt     = Dir->make<TH1F>("nerr_vs_evt", Form("run %06i: link %02i:%i:%i N err vs evt" ,RunNumber,Station,Dtc,Link), 1000, 0.,   5000000.);
    Hist->eflg_vs_evt     = Dir->make<TH1F>("eflg_vs_evt", Form("run %06i: link %02i:%i:%i eflag vs evt" ,RunNumber,Station,Dtc,Link), 1000, 0.,   5000000.);

    Hist->nh_vs_ch        = Dir->make<TH2F>("nh_vs_ch"  ,      Form("run %06i: link %02i:%i:%i nh vs ch"     ,RunNumber,Station,Dtc,Link),  100,0.,100., 50,0,50);
    Hist->nh_vs_adc1      = Dir->make<TH2F>("nh_vs_adc1",      Form("run %06i: link %02i:%i:%i nh vs adc"    ,RunNumber,Station,Dtc,Link),  100,0.,100., 50,0,50);

    Hist->dt0r_vs_ch      = Dir->make<TH2F>("dt0r_vs_ch_0",    Form("run %06i: link %02i:%i:%i dt0r vs ch[0]",RunNumber,Station,Dtc,Link),  100,0.,100.,2500,-25,25);
    Hist->dt1r_vs_ch      = Dir->make<TH2F>("dt1r_vs_ch_0",    Form("run %06i: link %02i:%i:%i dt1r vs ch[0]",RunNumber,Station,Dtc,Link),  100,0.,100.,2500,-25,25);

    Hist->dt0r01          = Dir->make<TH1F>("dt0r01",          Form("run %06i: link %02i:%i:%i dt0r01"       ,RunNumber,Station,Dtc,Link), 40000,-20000,20000);
    Hist->dt1r01          = Dir->make<TH1F>("dt1r01",          Form("run %06i: link %02i:%i:%i dt1r01"       ,RunNumber,Station,Dtc,Link), 40000,-20000,20000);

    Hist->nhits_vs_ich    = Dir->make<TH1F>("nh_vs_ich"  ,     Form("run %06i: link %02i:%i:%i nh vs ich"    ,RunNumber,Station,Dtc,Link),  100, 0.,   100.);
    Hist->nhits_vs_adc[0] = Dir->make<TH1F>("nh_vs_adc_0",     Form("run %06i: link %02i:%i:%i nh vs adc_0"  ,RunNumber,Station,Dtc,Link),  100, 0.,   100.);
    Hist->nhits_vs_adc[1] = Dir->make<TH1F>("nh_vs_adc_1",     Form("run %06i: link %02i:%i:%i nh vs adc_1"  ,RunNumber,Station,Dtc,Link),  100, 0.,   100.);

    Hist->dt0rc_vs_ch[0]  = Dir->make<TH2F>("dt0rc_vs_ch_0",   Form("run %06i: link %02i:%i:%i dt0rc vs ch[0], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);
    Hist->dt0rc_vs_ch[1]  = Dir->make<TH2F>("dt0rc_vs_ch_1",   Form("run %06i: link %02i:%i:%i dt0rc vs ch[1], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);

    Hist->dt1rc_vs_ch[0]  = Dir->make<TH2F>("dt1rc_vs_ch_0",   Form("run %06i: link %02i:%i:%i dt1rc vs ch[0], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);
    Hist->dt1rc_vs_ch[1]  = Dir->make<TH2F>("dt1rc_vs_ch_1",   Form("run %06i: link %02i:%i:%i dt1rc vs ch[1], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);

    Hist->dt0rc_vs_adc[0] = Dir->make<TH2F>("dt0rc_vs_adc_0",  Form("run %06i: link %02i:%i:%i dt0rc vs adc[0], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);
    Hist->dt0rc_vs_adc[1] = Dir->make<TH2F>("dt0rc_vs_adc_1",  Form("run %06i: link %02i:%i:%i dt0rc vs adc[1], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);

    Hist->dt1rc_vs_adc[0] = Dir->make<TH2F>("dt1rc_vs_adc_0",  Form("run %06i: link %02i:%i:%i dt1rc vs adc[0], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);
    Hist->dt1rc_vs_adc[1] = Dir->make<TH2F>("dt1rc_vs_adc_1",  Form("run %06i: link %02i:%i:%i dt1rc vs adc[1], ns",RunNumber,Station,Dtc,Link),  100,0.,100.,1000,-10,10);

    Hist->fs_vs_ich       = Dir->make<TProfile>("fs_vs_ich"  , Form("run %06i: link %02i:%i:%i fs vs ich"    ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0,  30);
    Hist->bl_vs_ich       = Dir->make<TProfile>("bl_vs_ich"  , Form("run %06i: link %02i:%i:%i bl vs ich"    ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0, 500);
    Hist->ph_vs_ich       = Dir->make<TProfile>("ph_vs_ich"  , Form("run %06i: link %02i:%i:%i ph vs ich"    ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0, 500);
    Hist->q_vs_ich        = Dir->make<TProfile>("q_vs_ich"   , Form("run %06i: link %02i:%i:%i Q vs ich"     ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0,1500);
    Hist->qt_vs_ich       = Dir->make<TProfile>("qt_vs_ich"  , Form("run %06i: link %02i:%i:%i Qt vs ich"    ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0, 500);
    Hist->qtq_vs_ich      = Dir->make<TProfile>("qtq_vs_ich" , Form("run %06i: link %02i:%i:%i Qt/Q vs ich"  ,RunNumber,Station,Dtc,Link),  100, 0.,   100.,0, 1);

    for (int i=0; i<kNChannels; i++) {
      art::TFileDirectory chan_dir = Dir->mkdir(Form("ch_%02i",i));
      book_channel_histograms(&chan_dir,RunNumber,&Hist->channel[i],Link,i);
    }
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
   
    Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"  ,RunNumber), 1000, 0.,   5000.);
    Hist->nbtot           = Dir->make<TH1F>("nbtot"      , Form("run %06i: nbytes total" ,RunNumber), 1000, 0., 100000.);
    Hist->nfrag           = Dir->make<TH1F>("nfrag"      , Form("run %06i: n fragments"  ,RunNumber),  100, 0.,    100.);
    Hist->fsize           = Dir->make<TH1F>("fsize"      , Form("run %06i: fragment size",RunNumber), 1000, 0., 100000.);
    Hist->valid           = Dir->make<TH1F>("valid"      , Form("run %06i: valid code"   ,RunNumber),  100, 0.,    100.);

    Hist->n_empty         = Dir->make<TH1F>("nempt"      , Form("run %06i: N(empty)"     ,RunNumber),  100, 0.,    100.);
    Hist->n_invalid_dr    = Dir->make<TH1F>("ninvr"      , Form("run %06i: N(invalid dr)",RunNumber),  100, 0.,    100.);
    Hist->n_corrupt       = Dir->make<TH1F>("ncorr"      , Form("run %06i: N(corrupt)"   ,RunNumber),  100, 0.,    100.);
    Hist->n_timeouts      = Dir->make<TH1F>("ntmo"       , Form("run %06i: N(timeouts)"  ,RunNumber),  100, 0.,    100.);
    Hist->n_overflows     = Dir->make<TH1F>("nover"      , Form("run %06i: N(overflows)" ,RunNumber),  100, 0.,    100.);

    Hist->error_code      = Dir->make<TH1F>("error_code" , Form("run %06i: error code"   ,RunNumber),  512, 0.,    512.);
    Hist->nerr_tot        = Dir->make<TH1F>("nerr_tot"   , Form("run %06i: N errors"     ,RunNumber), 1000, 0.,   2000.);
    Hist->nerr_vs_evt     = Dir->make<TH1F>("nerr_vs_evt", Form("run %06i: N err vs evt" ,RunNumber), 1000, 0.,   5000000.);
    Hist->eflg_vs_evt     = Dir->make<TH1F>("eflg_vs_evt", Form("run %06i: eflag vs evt" ,RunNumber), 1000, 0.,   5000000.);
  }

//-----------------------------------------------------------------------------
// for now - make the interface work with one station only
//-----------------------------------------------------------------------------
  void TrkFragmentAna::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
    
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;		// all events
    book_event_histset[ 1] = 0;	        // events with the error code = 0
    char folder_name[100];
    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory top_dir = tfs->mkdir(folder_name);

        _hist.event[i] = new EventHist_t;
        book_event_histograms(&top_dir,RunNumber,_hist.event[i]);
      }
    }
    
    int book_station_histset[kNStationHistSets];
    for (int i=0; i<kNStationHistSets; i++) book_station_histset[i] = 0;

    book_station_histset[ 0] = 1;		// all events
    book_station_histset[ 1] = 0;	        // events with the error code = 0

    for (int ist=0; ist<kNStationHistSets; ist++) {
      if (book_station_histset[ist] != 0) {
        sprintf(folder_name,"stn_%02i",ist);
        art::TFileDirectory stn_dir = tfs->mkdir(folder_name);
        _hist.station[ist] = new StationHist_t;

        for (int idtc=0; idtc<2; idtc++) {
          art::TFileDirectory dtc_dir = stn_dir.mkdir(Form("dtc%i",idtc));
          book_dtc_histograms(&dtc_dir,RunNumber,&_hist.station[ist]->dtc[idtc],ist,idtc);
          for (int i=0; i<_nActiveLinks[idtc]; i++) {
            int link  = (_activeLinks[idtc])->at(i);        // this assumes one station
            art::TFileDirectory roc_dir = dtc_dir.mkdir(Form("roc%i",link));
            book_roc_histograms(&roc_dir,RunNumber,&_hist.station[ist]->dtc[idtc].roc[link],ist,idtc,link);
          }
        }
      }
    }
        
    printf("[mu2e::TrkFragmentAna] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
// for run<=285, wasn't saving the event header and subheader, 64 bytes in total
// or 32 2-byte words
// book histograms only once
//-----------------------------------------------------------------------------
void TrkFragmentAna::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();

  if (_initialized != 0) return;
  _initialized = 1;

  // if (rn <= 285) _dataHeaderOffset =  0;
  // else           _dataHeaderOffset = 32;
//-----------------------------------------------------------------------------
// init timing offsets for known runs
// 'offset' is defined from the analysis of the dt0r_vs_ch_0
//-----------------------------------------------------------------------------
  double offset(0);

  if      (rn ==    281) offset = 15030; 
  else if (rn == 105023) offset =     0;
  else if (rn == 105026) offset = 10400;
  else if (rn == 105038) offset =  7964; 
  else if (rn == 105041) offset =  1128;
  else if (rn == 105042) offset =  1128;
  else if (rn == 105043) offset =  1128;
  else if (rn == 105044) offset =  1129;
  else if (rn == 105059) offset =     0; 
  else if (rn == 105060) offset = 11130; 
  else if (rn == 105066) offset = 10580;
  // 105067-1105070: four runs with the same settings, 2 ROCs, no data file for 105067 
  else if ((rn >= 105067) and (rn <= 105070)) offset = 0;
//-----------------------------------------------------------------------------
// add 20 ns to have the HV lanes distinct on the plot
//-----------------------------------------------------------------------------
  for (int i=0; i<kNChannels; i++) {
    if (_adc_index_0[i] < 48) _gen_offset[i] = 0;
    else                      _gen_offset[i] = offset;
  }
//-----------------------------------------------------------------------------
// as a last step, book histograms - need to know the number of active links
//-----------------------------------------------------------------------------
  book_histograms(rn);
}

//--------------------------------------------------------------------------------
  void TrkFragmentAna::beginJob() {
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::endJob() {
    printf("[mu2e::TrkFragmentAna] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data) {
    Hist->nhits->Fill(Data->nhits());
  }

//-----------------------------------------------------------------------------
  void TrkFragmentAna::fill_dtc_histograms(DtcHist_t* Hist, StationData_t* Sd, int IDtc) {
  }
    
//-----------------------------------------------------------------------------
  void TrkFragmentAna::fill_roc_histograms(RocHist_t* Hist, RocData_t* Rd) {
    
    Hist->nbytes->Fill      (Rd->nbytes);
    Hist->npackets->Fill    (Rd->npackets);
    Hist->nhits->Fill       (Rd->nhits);
    Hist->valid->Fill       (Rd->valid);
    Hist->error_code->Fill  (Rd->error_code);

    Hist->n_empty->Fill     (Rd->n_empty);
    Hist->n_invalid_dr->Fill(Rd->n_invalid_dr);
    Hist->n_corrupt->Fill   (Rd->n_corrupt);
    Hist->n_timeouts->Fill  (Rd->n_timeouts);
    Hist->n_overflows->Fill (Rd->n_overflows);

    Hist->nerr_tot->Fill    (Rd->nerr_tot);
    Hist->nerr_vs_evt->Fill (_edata._event->event(),Rd->nerr_tot);
    int eflg = (Rd->nerr_tot > 0);
    Hist->eflg_vs_evt->Fill (_edata._event->event(),eflg        );

    if (Rd->nhits > 0) {
      Hist->dt0r01->Fill  (Rd->dt0r01);
      Hist->dt1r01->Fill  (Rd->dt1r01);
    }
//-----------------------------------------------------------------------------
// fill channel histograms
//-----------------------------------------------------------------------------
    for (int ich=0; ich<kNChannels; ich++) {
      ChannelHist_t* hch = &Hist->channel[ich];
      ChannelData_t* chd = &Rd->channel[ich];

      int nh = chd->nhits();
      hch->nhits->Fill(nh);
      if (nh == 0)                                          continue;
      
      int fpga  = _adc_index_1[ich] / 48;
      hch->dt0r->Fill(chd->dt0r);
      hch->dt1r->Fill(chd->dt1r);
            
      Hist->dt0r_vs_ch->Fill(ich,chd->dt0r);
      Hist->dt1r_vs_ch->Fill(ich,chd->dt1r);

      Hist->dt0rc_vs_ch[fpga]->Fill(ich,chd->dt0r_c);
      Hist->dt1rc_vs_ch[fpga]->Fill(ich,chd->dt1r_c);
            
      int iadc = _adc_index_1[ich];
      Hist->dt0rc_vs_adc[fpga]->Fill(iadc,chd->dt0r_c);
      Hist->dt1rc_vs_adc[fpga]->Fill(iadc,chd->dt1r_c);
//-----------------------------------------------------------------------------
// there are hits in this channel, store not more than kMaxNHWfPerChannel,
// but want to handle all hits correctly
//-----------------------------------------------------------------------------
      for (int ih=0; ih<nh; ih++) {
        TrackerDataDecoder::TrackerDataPacket* hit = chd->hit[ih];

        uint32_t corr_tdc0 = correctedTDC(hit->TDC0());
        uint32_t corr_tdc1 = correctedTDC(hit->TDC1());

        float t0_us = corr_tdc0*_tdc_bin;
        float t1_us = corr_tdc1*_tdc_bin;
        
        float t0_ns = corr_tdc0*_tdc_bin_ns;
        float t1_ns = corr_tdc1*_tdc_bin_ns;
        
        hch->time[0]->Fill(t0_us);       // in us
        hch->time[1]->Fill(t1_us);       // in us

        hch->dt01[0]->Fill(t0_ns-t1_ns); // in ns
        hch->dt01[1]->Fill(t0_ns-t1_ns); // in ns

        hch->t0  [0]->Fill(t0_ns);
        hch->t0  [1]->Fill(t1_ns);

        hch->t1  [0]->Fill(_timeWindow-t0_ns);
        hch->t1  [1]->Fill(_timeWindow-t1_ns);

        hch->tot [0]->Fill(hit->TOT0);
        hch->tot [1]->Fill(hit->TOT1);
        hch->pmp    ->Fill(hit->PMP);
//-----------------------------------------------------------------------------
// assume that nsamples is known and try to process the waveform        
//-----------------------------------------------------------------------------
        int   nsamples = 15+12*(_nADCPackets-1);
        float wform[50];
          
        WfParam_t wpar;
        unpack_adc_waveform(hit,wform,&wpar);
          
        chd->wp.push_back(wpar);
//-----------------------------------------------------------------------------
// fill waveform histograms only for the firstkMaxNHWfPerChannel  channels
//-----------------------------------------------------------------------------
        if (ih < kMaxNHWfPerChannel) {
          hch->raw_wf[ih]->Reset();
          hch->wf    [ih]->Reset();
          for (int is=0; is<nsamples; is++) {
            hch->raw_wf[ih]->Fill(is,wform[is]+wpar.bl);
            hch->wf    [ih]->Fill(is,wform[is]);
          }
            // also set bin errors to zero
          int nb =  hch->wf[ih]->GetNbinsX();
          for (int ib=0; ib<nb; ib++) {
            hch->raw_wf[ih]->SetBinError(ib+1,0);
            hch->raw_wf[ih]->SetOption("HIST");
            hch->wf[ih]->SetBinError(ib+1,0);
            hch->wf[ih]->SetOption("HIST");
          }
        }
        else {
          int nmax = kMaxNHWfPerChannel;
          if (DebugBit(0) != 0) {
            TLOG(TLVL_WARNING) << "DTC:"     << Rd->dtc_id  << " link:" << Rd->link
                               << " ih="     << ih 
                               << " ch:0x"   << std::hex << ich 
                               << " hits > " << std::dec << nmax << " do not make extra waveform histograms"
                               << std::endl;
          }
        }
//-----------------------------------------------------------------------------
// reconstructed waveform parameters
//-----------------------------------------------------------------------------
        hch->fsample->Fill(wpar.fs);
        hch->bline->Fill(wpar.bl);
        hch->pheight->Fill(wpar.ph);
        hch->q->Fill(wpar.q);
        hch->qt->Fill(wpar.qt);
        hch->qtq->Fill(wpar.qt/(wpar.q+1e-12));
        
        Hist->fs_vs_ich->Fill(ich,wpar.fs);
        Hist->bl_vs_ich->Fill(ich,wpar.bl);
        Hist->ph_vs_ich->Fill(ich,wpar.ph);
        Hist->q_vs_ich->Fill(ich,wpar.q);
        Hist->qt_vs_ich->Fill(ich,wpar.qt);
        Hist->qtq_vs_ich->Fill(ich,wpar.qt/(wpar.q+1e-12));
      }
//-----------------------------------------------------------------------------
// time distance between the two sequential hits - need at least two
//-----------------------------------------------------------------------------
      for (int ih=1; ih<nh; ih++) {
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
// 
//-----------------------------------------------------------------------------
    for (int ich=0; ich<kNChannels; ich++) {
      int ind_0 = _adc_index_0[ich];
      int ind_1 = _adc_index_1[ich];

      int nh    = Rd->channel[ich].nhits();
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
  void TrkFragmentAna::fill_event_histograms(EventHist_t* Hist, EventData_t* Ed) {

    Hist->error_code->Fill(Ed->error_code);
    Hist->nerr_tot->Fill(Ed->nerr_tot);
    Hist->valid->Fill(Ed->valid);

    Hist->nbtot->Fill(Ed->nbtot);
    Hist->nhits->Fill(Ed->nhtot);
    Hist->nfrag->Fill(Ed->nfrag);

    Hist->n_empty->Fill     (Ed->n_empty);
    Hist->n_invalid_dr->Fill(Ed->n_invalid_dr);
    Hist->n_corrupt->Fill   (Ed->n_corrupt);
    Hist->n_timeouts->Fill  (Ed->n_timeouts);
    Hist->n_overflows->Fill (Ed->n_overflows);

    for (int i=0; i<Ed->nfrag; i++) {
      int fsize = Ed->fragments[i].nbytes;
      Hist->fsize->Fill(fsize);
    }

    Hist->nerr_vs_evt->Fill  (Ed->_event->event(),Ed->nerr_tot);
    int eflg = (Ed->nerr_tot > 0);
    Hist->eflg_vs_evt->Fill  (Ed->_event->event(),eflg        );
  }

//-----------------------------------------------------------------------------
// to handle more than one station, this function will need to evolve
//-----------------------------------------------------------------------------
  int TrkFragmentAna::fill_station_histograms(StationHist_t* Hist, EventData_t* Data) {
    
    for (int idtc=0; idtc<2; idtc++) {
      for (int ir=0; ir<_nActiveLinks[idtc]; ir++) {
        int link = _activeLinks[idtc]->at(ir);
        fill_roc_histograms(&Hist->dtc[idtc].roc[link],&Data->station[0].roc[idtc][link]);
      }
    }
    
    return 0;
  }
    
//-----------------------------------------------------------------------------
// fill_roc_histograms also fills the channel histograms
// if in error, only histogram the error code
//-----------------------------------------------------------------------------
  int TrkFragmentAna::fill_histograms() {

//-----------------------------------------------------------------------------
// later, all error handling will move to validate_roc_data()
//-----------------------------------------------------------------------------
    fill_event_histograms(_hist.event[0],&_edata);
    // if (_edata.error_code == 0) fill_event_histograms(_hist.event[1],&_edata);


    fill_station_histograms(_hist.station[0],&_edata);
    // if (_edata.error_code == 0) fill_station_histograms(_hist.station[1],&_edata);
    
    return 0;
  }

  
//-----------------------------------------------------------------------------
// a fragment may have multiple ROC blocks
//-----------------------------------------------------------------------------
  int TrkFragmentAna::validate_roc_data(RocDataHeaderPacket_t* Rdh) {
    int rc(0);

    return rc;
  }

//-----------------------------------------------------------------------------
// a fragment may have multiple ROC blocks
//-----------------------------------------------------------------------------
  void TrkFragmentAna::analyze_fragment(const art::Event& Evt, const artdaq::Fragment* Fragment) {

    short* fdata = (short*) Fragment->dataBegin();

    _edata.fragments.push_back(FragmentData_t());
//-----------------------------------------------------------------------------
// pointer to the last one
//-----------------------------------------------------------------------------
    FragmentData_t* fdt = &_edata.fragments.back();
//-----------------------------------------------------------------------------
// fragment size is specified in longs and includes service data, don't use
//-----------------------------------------------------------------------------
    fdt->nbytes  = fdata[0];
    if (fdata[0] > _maxFragmentSize) {
      _edata.error_code |= kNBytesErrorBit;
      _edata.nerr_tot   += 1;
      if (DebugBit(0) != 0) { 
        TLOG(TLVL_ERROR) << Form("event %6i:%8i:%8i : ERROR_CODE:0x%04x in %s: fdt->nbytes= %i, BAIL OUT\n",
                                 Evt.run(),Evt.subRun(),Evt.event(),_edata.error_code,__func__,fdt->nbytes);
      }
    }
//-----------------------------------------------------------------------------
// somewhere here handle the DTC ID
//-----------------------------------------------------------------------------
    _station          = 0;
    StationData_t* sd = &_edata.station[_station];
//-----------------------------------------------------------------------------
// start handling the ROC data
// need mapping of the DTC ID to the plane number
// for now, just label the DTCs by their PCIE address
//-----------------------------------------------------------------------------
    SubEventHeader_t* sh = (SubEventHeader_t*) fdata;
    int dtc_id           = sh->dtcID;
    
    short* first_address = fdata+sizeof(SubEventHeader_t)/2; // _dataHeaderOffset; // offset is specified in 2-byte words
    short* last_address  = fdata+fdt->nbytes/2;     // 

    while (first_address < last_address) {
//-----------------------------------------------------------------------------
// next ROC
//-----------------------------------------------------------------------------
      RocDataHeaderPacket_t* dh = (RocDataHeaderPacket_t*) first_address;
      int link      = dh->linkID;
      validate_roc_data(dh);
//-----------------------------------------------------------------------------
// check link number
//-----------------------------------------------------------------------------
      int found = 0;
      for (int i=0; i<_nActiveLinks[_idtc]; i++) {
        if (_activeLinks[_idtc]->at(i) == link) {
          found = 1;
          break;
        }
      }

      RocData_t* rd = &_edata.station[_station].roc[_idtc][link];

      if (found == 0) {
//-----------------------------------------------------------------------------
// in some cases, this one could be just a warning
//-----------------------------------------------------------------------------
        _edata.error_code |= kLinkIDErrorBit;
        _edata.nerr_tot   += 1;
        rd->nerr_tot      += 1;
        if (DebugBit(0) != 0) {
          TLOG(TLVL_ERROR) << Form("event %6i:%8i:%8i : ERROR_CODE:0x%04x in %s: link=%i not defined as active, BAIL OUT\n",
                                   Evt.run(),Evt.subRun(),Evt.event(),kLinkIDErrorBit,__func__,link);
        }
      }
//-----------------------------------------------------------------------------
// check ROC status
//-----------------------------------------------------------------------------
      if (dh->empty()     ) rd->n_empty      += 1;
      if (dh->invalid_dr()) rd->n_invalid_dr += 1;
      if (dh->corrupt()   ) rd->n_corrupt    += 1;
      if (dh->timeout()   ) rd->n_timeouts   += 1;
      if (dh->overflow()  ) rd->n_overflows  += 1;

      rd->error_code = dh->error_code();
//-----------------------------------------------------------------------------
// for a given FPGA, a reference channel is the first channel in the readout order
//-----------------------------------------------------------------------------
      ChannelData_t* ref_ch[2];

      ref_ch[0]     = &rd->channel[_referenceChannel[link][0]];
      ref_ch[1]     = &rd->channel[_referenceChannel[link][1]];
        
      rd->nbytes    = dh->byteCount;
      rd->npackets  = dh->packetCount;
      rd->dtc_id    = dtc_id;
//-----------------------------------------------------------------------------
// for now, assume that all hits in the run have the same number of packets per hit
// take that from the first hit
//-----------------------------------------------------------------------------
      rd->nhits     = dh->packetCount/(_nADCPackets+1);         //  printf("nhits : %3i\n",nhits);
      rd->valid     = dh->valid;
      rd->dt0r01    = -1.e12;
      rd->dt1r01    = -1.e12;
      
      _edata.nhtot += rd->nhits;
      _edata.valid += dh->valid*10;
      
      for (int ihit=0; ihit<rd->nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
        TrackerDataDecoder::TrackerDataPacket* hit ;
        int offset = ihit*(8+8*_nADCPackets);
        hit     = (TrackerDataDecoder::TrackerDataPacket*) (first_address+0x08+offset);

        if (DebugBit(5) != 0) {
          printf("offset : 0x%04x\n",offset);
          print_hit(hit);
        }
        
        if (hit->ErrorFlags != 0) {  // 4 bits
          if (DebugBit(0) == 1) {
          TLOG(TLVL_ERROR) << "DTC:" << rd->dtc_id  << " link:" << rd->link
                           << " ih=" << ihit
                           << " first_address:0x" << std::hex << first_address
                           << " offset:0x" << std::hex << offset
                           << " HIT ERROR FLAG:0x" << std::hex << hit->ErrorFlags
                           << std::endl;
          }
          _edata.error_code |= kHitErrorBit;
          _edata.nerr_tot   += 1;

          rd->error_code    |= kHitErrorBit;
          rd->nerr_tot      += 1;
        }
        
        if (hit->NumADCPackets != _nADCPackets) {
          // this is an error which doesn't allow to proceed looping over the hits
          if (DebugBit(0) != 0) {
            TLOG(TLVL_ERROR) << "DTC:" << rd->dtc_id  << " link:" << rd->link
                             << " ih=" << ihit
                             << " first_address:0x" << std::hex << first_address
                             << " offset:0x" << std::hex << offset
                             << " WRONG NUMBER OF ADC PACKETS: " << std::dec << hit->NumADCPackets
                             << " BAIL OUT" << std::endl;
          }
          
          if (DebugBit(4) != 0) print_hit(hit);

          _edata.error_code |= kNWfsErrorBit;
          _edata.nerr_tot   += 1;
          rd->error_code    |= kNWfsErrorBit;
          rd->nerr_tot      += 1;
          break;
        }

        int ich = hit->StrawIndex;

        if (ich > 128) ich = ich-128;

        if (ich > 95) {
//-----------------------------------------------------------------------------
// non existing channel ID : flag an error, don't save the hit, but continue
//-----------------------------------------------------------------------------
          _edata.error_code |= kChIDErrorBit;
          _edata.nerr_tot   += 1;
          rd->error_code    |= kChIDErrorBit;
          rd->nerr_tot      += 1;

          if (DebugBit(0) != 0) {
            TLOG(TLVL_ERROR) << Form("event %i:%i:%i : ERROR_CODE:0x%04x in %s: link = %i offset: 0x%04x hit->StrawIndex = 0x%04x\n",
                                     Evt.run(),Evt.subRun(),Evt.event(),kChIDErrorBit,__func__,link,offset,hit->StrawIndex);
          }
        }
        else {
          ChannelData_t* chd = &rd->channel[ich];
          chd->hit.push_back(hit);
        }
//-----------------------------------------------------------------------------
// hit pattern words should repeat
// for the 'checkerboard', pattern #4, it is as follows
//-----------------------------------------------------------------------------
        if (_validateAdcPatterns >= 0) {
          uint16_t pattern[4] = {0x56aa, 0x2aa5, 0xa955, 0x155a};
          uint16_t* w = (uint16_t*) hit;
          int nw = 2+_nADCPackets*8;
                                        // finding first word - they can move
          int offset = -1;
          int loc    = 6;
          for (int i=0; i<4; i++) {
            if (w[loc] == pattern[i]) {
              offset = i;
              break;
            }
          }

          if (offset == -1) {
            TLOG(TLVL_ERROR) << Form("event %i:%i:%i", Evt.run(),Evt.subRun(),Evt.event())
                             << Form(", : WRONG ADC PATTERN 0x%04x in %s: ",kAdcPatternErrorBit,__func__)
                             << Form(" link = %i offset: 0x%04x hit->StrawIndex = 0x%04x\n",link,offset,hit->StrawIndex);
            
            _edata.error_code |= kAdcPatternErrorBit;
            _edata.nerr_tot   += 1;
            rd->error_code    |= kAdcPatternErrorBit;
            rd->nerr_tot      += 1;
          }
          else {
//-----------------------------------------------------------------------------
// starting pattern found, validating - the first word (#0) already checked,
// thus starting from 1 
//-----------------------------------------------------------------------------
            for (int i=1; i<nw; i++) {
              int ipat = (offset+i) % 4;
              if (w[loc+i] == pattern[ipat]) continue;
//-----------------------------------------------------------------------------
// in trouble
//-----------------------------------------------------------------------------
            TLOG(TLVL_ERROR) << Form("event %i:%i:%i", Evt.run(),Evt.subRun(),Evt.event())
                             << Form(", : WRONG ADC PATTERN 0x%04x in %s: ",kAdcPatternErrorBit,__func__)
                             << Form(" link = %i offset: 0x%04x hit->StrawIndex = 0x%04x\n",link,offset,hit->StrawIndex);
            
            _edata.error_code |= kAdcPatternErrorBit;
            _edata.nerr_tot   += 1;
            rd->error_code    |= kAdcPatternErrorBit;
            rd->nerr_tot      += 1;
            }
          }
        }
      }
//-----------------------------------------------------------------------------
// hits in all channels counted, modulo those with corrupted channel IDs
// time difference between this channel and a reference channel
// determined using the first hit only
//-----------------------------------------------------------------------------
      for (int i=0; i<kNChannels; i++) {
        ChannelData_t* chd = &rd->channel[i];

        int nh   = chd->nhits();
        int fpga = _adc_index_1[i] / 48;

        ChannelData_t* rch = ref_ch[fpga];
//-----------------------------------------------------------------------------
// in most cases, the number of hits in the reference channel should be greater 
// than the number of channels in any other channel of a given FPGA
//-----------------------------------------------------------------------------
        int iref = _referenceChannel[link][fpga];
        int nhr  = rd->channel[iref].nhits();
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
      if ((rd->ref_ch[0]->nhits() > 0) and (rd->ref_ch[1]->nhits() > 0)) {
        int t0r0   = correctedTDC(rd->ref_ch[0]->hit[0]->TDC0());
        int t1r0   = correctedTDC(rd->ref_ch[0]->hit[0]->TDC1());
        
        int t0r1   = correctedTDC(rd->ref_ch[1]->hit[0]->TDC0());
        int t1r1   = correctedTDC(rd->ref_ch[1]->hit[0]->TDC1());
        
        rd->dt0r01 = (t0r0-t0r1)*_tdc_bin_ns;        // convert to ns  
        rd->dt1r01 = (t1r0-t1r1)*_tdc_bin_ns;        // convert to ns  
      }
//-----------------------------------------------------------------------------
// update station counters
//-----------------------------------------------------------------------------
      sd->n_empty      += rd->n_empty; 
      sd->n_invalid_dr += rd->n_invalid_dr; 
      sd->n_corrupt    += rd->n_corrupt; 
      sd->n_timeouts   += rd->n_timeouts; 
      sd->n_overflows  += rd->n_overflows; 
//-----------------------------------------------------------------------------
// address in 2-byte words (N(data packets)+data header packet)
//-----------------------------------------------------------------------------
      first_address    += (dh->packetCount + 1)*8;
    }
  }

//-----------------------------------------------------------------------------
int TrkFragmentAna::init_event(const art::Event& AnEvent) {
  _edata._event = &AnEvent;

  _edata.nbtot      = 0;
  _edata.nhtot      = 0;
  _edata.nfrag      = 0;
  _edata.valid      = 0;
//-----------------------------------------------------------------------------
// reset error counters
//-----------------------------------------------------------------------------
  _edata.error_code    = 0;
  _edata.nerr_tot      = 0;

  _edata.n_nb_errors   = 0;
  _edata.n_nwfs_errors = 0;
  _edata.n_chid_errors = 0;
  _edata.n_nchh_errors = 0;

  _edata.n_empty       = 0;
  _edata.n_invalid_dr  = 0;
  _edata.n_corrupt     = 0;
  _edata.n_timeouts    = 0;
  _edata.n_overflows   = 0;

  _edata.fragments.clear();
//-----------------------------------------------------------------------------
// clear all counters
//-----------------------------------------------------------------------------
  for (int ist=0; ist<kMaxStations; ist++) {
    StationData_t* sd = &_edata.station[ist];
    
    sd->n_empty       = 0;
    sd->n_invalid_dr  = 0;
    sd->n_corrupt     = 0;
    sd->n_timeouts    = 0;
    sd->n_overflows   = 0;
    
    for (int dtc=0; dtc<2; dtc++) {
      for (int lnk=0; lnk<6; lnk++) {
        RocData_t* rd    = &sd->roc[dtc][lnk];
        rd->n_empty      = 0;
        rd->n_invalid_dr = 0;
        rd->n_corrupt    = 0;
        rd->n_timeouts   = 0;
        rd->n_overflows  = 0;
        rd->nerr_tot     = 0;
        rd->error_code   = 0;
        for (int ich=0; ich<kNChannels; ich++) {
          ChannelData_t* chd = &rd->channel[ich];
          chd->hit.clear();
          chd->wp.clear();
        }
      }
    }
  }
  return 0;
}


//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TrkFragmentAna::analyze(const art::Event& AnEvent) {

  init_event(AnEvent);

  auto handle = AnEvent.getValidHandle<std::vector<artdaq::Fragment> >(_trkfCollTag);
//-----------------------------------------------------------------------------
// proxy for event histograms
//-----------------------------------------------------------------------------
  if (_diagLevel > 0) {
    printf(" Event : %06i:%08i%08i\n", AnEvent.run(),AnEvent.subRun(),AnEvent.event());
  }
  
  _idtc = 0;
  for (const artdaq::Fragment& frag : *handle) {
//-----------------------------------------------------------------------------
// different fragments correspond to different DTCs, somewhere there should be the DTC ID
//-----------------------------------------------------------------------------
    ushort* buf = (ushort*) (frag.dataBegin());
    int nbytes  = buf[0];
    int fsize   = frag.sizeBytes();

    if (nbytes < 2) {
      _edata.error_code  |= kNBytesErrorBit;
      _edata.nerr_tot    += 1;
      _edata.n_nb_errors += 1;
      
      TLOG(TLVL_DEBUG) << Form("event %6i:%8i:%8i : ERROR_CODE:0x%04x nbytes=%i EMPTY_FRAGMENT",
                               AnEvent.run(),AnEvent.subRun(),AnEvent.event(),
                               kNBytesErrorBit,nbytes);
    }

    _edata.nfrag += 1;
    _edata.nbtot += nbytes;        // including the artdaq part

    if (_diagLevel > 2) {
      printf("%s\n",Form("---------- fragment # %3i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x\n",
                         _idtc,nbytes,fsize,_edata.error_code));
      print_fragment(&frag,nbytes/2);
    }

    //    if ((_edata.error == 0) and _analyzeFragments) analyze_fragment(event,&frag);
    if (_analyzeFragments) analyze_fragment(AnEvent,&frag);
    _idtc++;
  }
//-----------------------------------------------------------------------------
// update per-event counters
//-----------------------------------------------------------------------------
  for (int ist=0; ist<kMaxStations; ist++) {
    StationData_t* sd    = &_edata.station[ist];
    
    _edata.n_empty      += sd->n_empty; 
    _edata.n_invalid_dr += sd->n_invalid_dr; 
    _edata.n_corrupt    += sd->n_corrupt; 
    _edata.n_timeouts   += sd->n_timeouts; 
    _edata.n_overflows  += sd->n_overflows; 
  }
//-----------------------------------------------------------------------------
// proxy for event histograms
//-----------------------------------------------------------------------------
  if (_diagLevel > 1) {
    if ((_edata.nbtot >= _minNBytes) and (_edata.nbtot <= _maxNBytes)) {
      TLOG(TLVL_DEBUG) << Form(" Run : %5i subrun: %5i event: %8i nfrag: %3i nbytes: %5i\n", 
                               AnEvent.run(),AnEvent.subRun(),AnEvent.event(),
                               _edata.nfrag, _edata.nbtot);
    }
  }
//-----------------------------------------------------------------------------
// print debug information
//-----------------------------------------------------------------------------
  debug(AnEvent);
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
    auto h = AnEvent.getValidHandle<std::vector<artdaq::Fragment>>("daq:TRKDTC");

    for (const artdaq::Fragment& frag : *h) {
      int *buf  = (int*) (frag.dataBegin());
      int nreg  = buf[0];
      int fsize = frag.sizeBytes();
      printf("%s: -------- DTC registers dump n(reg)=%5i size: %5i\n",__func__,nreg,fsize);
      print_fragment(&frag,2+4*nreg);
    }
  }
//-----------------------------------------------------------------------------
// finally, if requested, go into interactive mode, 
// fInteractiveMode = 0 : do not stop (default)
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
  TModule::analyze(AnEvent);
}

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TrkFragmentAna::print_message(const char* Message) {
  printf("TrkFragmentAna: event %6i:%8i%8i %s\n",
         _edata._event->run(),
         _edata._event->subRun(),
         _edata._event->event(),
         Message);
}

//-----------------------------------------------------------------------------
// NWords : the number of short words
//-----------------------------------------------------------------------------
  void TrkFragmentAna::print_fragment(const artdaq::Fragment* Frag, int NWords) {
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
// only need to print a couple of packets
//-----------------------------------------------------------------------------
void TrkFragmentAna::print_hit(const TrackerDataDecoder::TrackerDataPacket* Hit) {
  uint16_t* dat = (uint16_t*) Hit;

  int loc(0);
  for (int i=0; i<16; i++) {
    printf(" 0x%04x", dat[i]);
    int res = (loc+1) % 8;
    if (res == 0) {
      printf("\n");
      loc = 0;
    }
    else {
      loc = loc+1;
    }
  }
  if (loc != 0) printf("\n");
    
}
  
void TrkFragmentAna::debug(const art::Event& AnEvent) {
  
  auto handle = AnEvent.getValidHandle<std::vector<artdaq::Fragment> >(_trkfCollTag);

  int ifrag = 0;
  for (const artdaq::Fragment& frag : *handle) {
    ushort* buf = (ushort*) (frag.dataBegin());
    int nbytes  = buf[0];
    int fsize   = frag.sizeBytes();
    
    if (DebugBit(0) == 1) {
      print_message(Form("bit:000: fragment # %3i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i",
                         ifrag,nbytes,fsize,_edata.error_code,_edata.nerr_tot));
      print_fragment(&frag,nbytes/2);
    }

    if ((DebugBit(3) & 0x1) and (_edata.nerr_tot > _minNErrors)) {
      print_message(Form("bit:003: fragment # %3i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i\n",
                         ifrag,nbytes,fsize,_edata.error_code,_edata.nerr_tot));

      if (DebugBit(3) & 0x2) print_fragment(&frag,nbytes/2);
    }
    ifrag++;
  }

  if ((DebugBit(2) == 1) and (_edata.error_code == _errorCode)) {
    print_message(Form("bit:002: ERROR_CODE: 0x%04x NERR_TOT: %5i",
                       _edata.error_code,_edata.nerr_tot));
  }

}
  
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TrkFragmentAna)
