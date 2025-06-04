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
#define TRACE_NAME "StationAna"

#include "daqana/mod/StationAna_module.hh"
#include "daqana/mod/TrkPanelMap_t.hh"

namespace {
  const TrkPanelMap_t* _panel_map[200][6];   // indexing - offline: [plane][panel]
};

namespace mu2e {
  
// ======================================================================

  StationAna::StationAna(fhicl::ParameterSet const& PSet) : 
    THistModule            (PSet,PSet.get<fhicl::ParameterSet>("THistModule"),"StationAna") ,
    _debugMode             (PSet.get<int>             ("debugMode"             )),
    _shCollTag             (PSet.get<art::InputTag>   ("shCollTag"             )),
    _maxDt                 (PSet.get<float>           ("maxDt"                 )),
    _minEDep               (PSet.get<float>           ("minEDep"               ))
  {
    _initialized         = 0;
  }

//-----------------------------------------------------------------------------
// I : channel number
//-----------------------------------------------------------------------------
  void StationAna::book_straw_histograms(art::TFileDirectory* Dir, int RunNumber, StrawHist_t* Hist, int Mnid, int I) {

    // int ipanel       = 6*Dtc+Link;
    // const char* name = _panelName[ipanel].data();
    
    Hist->nsht   = Dir->make<TH1F>(Form("ch_%02i_nhits",I),Form("run %06i: MN%3i ch %02i nhits"  ,RunNumber,Mnid,I), 300,  -0.5, 299.5);
    Hist->tcal    = Dir->make<TH1F>(Form("ch_%02i_tcal" ,I),Form("run %06i: MN%3i ch %02i edep "  ,RunNumber,Mnid,I),1000,   0 , 100000);
    Hist->dtch    = Dir->make<TH1F>(Form("ch_%02i_dtch" ,I),Form("run %06i: MN%3i ch %02i dtCH "  ,RunNumber,Mnid,I), 400, -100, 100);
    Hist->edep    = Dir->make<TH1F>(Form("ch_%02i_edep" ,I),Form("run %06i: MN%3i ch %02i edep "  ,RunNumber,Mnid,I), 100, -0.002, 0.008);

    Hist->nshg    = Dir->make<TH1F>(Form("ch_%02i_nhitsg",I),Form("run %06i: MN%3i ch %02i nhitsG"  ,RunNumber,Mnid,I),  300,  -0.5, 299.5);
    Hist->tcalg    = Dir->make<TH1F>(Form("ch_%02i_tcalg" ,I),Form("run %06i: MN%3i ch %02i edepG "  ,RunNumber,Mnid,I),1000,   0 , 100000);
    Hist->dtchg    = Dir->make<TH1F>(Form("ch_%02i_dtchg" ,I),Form("run %06i: MN%3i ch %02i dtCHG "  ,RunNumber,Mnid,I), 400, -100, 100);
    Hist->edepg    = Dir->make<TH1F>(Form("ch_%02i_edepg" ,I),Form("run %06i: MN%3i ch %02i edepG "  ,RunNumber,Mnid,I), 100, -0.002, 0.008);
  }

//-----------------------------------------------------------------------------
  void StationAna::book_panel_histograms(art::TFileDirectory* Dir, int RunNumber, PanelHist_t* Hist, int Mnid) {
    
    Hist->nsht    = Dir->make<TH1F>(Form("nsht" ),Form("run %06i: MN%3i nshT"  ,RunNumber,Mnid), 300,  -0.5, 299.5);
    Hist->tcal    = Dir->make<TH1F>(Form("tcal" ),Form("run %06i: MN%3i tcal"  ,RunNumber,Mnid), 1000, 0, 100000);
    Hist->dtch    = Dir->make<TH1F>(Form("dtch" ),Form("run %06i: MN%3i dtch"  ,RunNumber,Mnid), 400, -100, 100);
    Hist->edep    = Dir->make<TH1F>(Form("edep" ),Form("run %06i: MN%3i edep"  ,RunNumber,Mnid), 100, -0.002, 0.008);
    Hist->occup   = Dir->make<TH1F>(Form("occup"),Form("run %06i: MN%3i occup" ,RunNumber,Mnid), 100,  0,   100);

    Hist->nshg    = Dir->make<TH1F>(Form("nshg"  ),Form("run %06i: MN%3i nshG"  ,RunNumber,Mnid), 100,  -0.5,  99.5);
    Hist->tcalg   = Dir->make<TH1F>(Form("tcalg" ),Form("run %06i: MN%3i tcalG"  ,RunNumber,Mnid), 1000, 0, 100000);
    Hist->dtchg   = Dir->make<TH1F>(Form("dtchg" ),Form("run %06i: MN%3i dtchG"  ,RunNumber,Mnid), 400, -100, 100);
    Hist->edepg   = Dir->make<TH1F>(Form("edepg" ),Form("run %06i: MN%3i edepG"  ,RunNumber,Mnid), 100, -0.002, 0.008);
    Hist->occupg  = Dir->make<TH1F>(Form("occupg"),Form("run %06i: MN%3i occupG" ,RunNumber,Mnid), 100,  0,   100);


    for (int i=0; i<kNStraws; i++) {
      art::TFileDirectory straw_dir = Dir->mkdir(Form("str_%02i",i));
      book_straw_histograms(&straw_dir,RunNumber,&Hist->straw[i],Mnid,i);
    }
  }

//-----------------------------------------------------------------------------
  void StationAna::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
    Hist->evt     = Dir->make<TH1F>("evt" ,Form("run %06i: nevents",RunNumber), 1000,  0.,    1e7);
    Hist->nsht    = Dir->make<TH1F>("nsht",Form("run %06i: nshT"   ,RunNumber),  300, -0.5, 299.5);
    Hist->nshg    = Dir->make<TH1F>("nshg",Form("run %06i: nshG"   ,RunNumber),  100, -0.5,  99.5);
    Hist->nshdt   = Dir->make<TH1F>("nshdt",Form("run %06i: nshDT" ,RunNumber),  100, -0.5,  99.5);
  }

//-----------------------------------------------------------------------------
// for now - make the interface work with one station only
//-----------------------------------------------------------------------------
  void StationAna::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
    
    TH1::AddDirectory(kFALSE);
    
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
    
    int book_panel_histset[kNPanelHistSets];
    for (int iset=0; iset<kNPanelHistSets; iset++) book_panel_histset[iset] = 0;

    book_panel_histset[ 0] = 1;		// all events
    
    for (int iset=0; iset<kNPanelHistSets; ++iset) {
      if (book_panel_histset[iset] != 0) {
        sprintf(folder_name,"pnlset_%02i",iset);
        art::TFileDirectory set_dir = tfs->mkdir(folder_name);

        for (int ip=0; ip<12; ++ip) {
          //_hist.panel_set[iset].panel[ip] = new PanelHist_t;
          int mnid = _edata.panel[ip].mnid;
          art::TFileDirectory panel_dir = set_dir.mkdir(Form("MN%i",mnid));
          book_panel_histograms(&panel_dir,RunNumber,&_hist.panel_set[iset].panel[ip],mnid);
        }
      }
    }
        
    printf("[mu2e::StationAna] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
// for run<=285, wasn't saving the event header and subheader, 64 bytes in total
// or 32 2-byte words
// book histograms only once
//-----------------------------------------------------------------------------
  void StationAna::beginRun(const art::Run& ArtRun) {
    int rn  = ArtRun.run();
    
    if (_initialized != 0) return;
    _initialized = 1;
//-----------------------------------------------------------------------------
// art may call beginRun on each input file
// the following are the steps which need to be taken just once
// initialize the panel map
//-----------------------------------------------------------------------------
    for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
      int plane = tpm->plane;
      int panel = tpm->panel;
      _panel_map[plane][panel] = tpm;

      int idtc   = tpm->dtc % 2;
      int ilink  = tpm->link;
      int ipanel = 6*idtc+ilink;
      
      _edata.panel[ipanel].mnid = tpm->mnid;
    }
//-----------------------------------------------------------------------------
// as a last step, book histograms - need to know the number of active links
//-----------------------------------------------------------------------------
    book_histograms(rn);
  }
  
//--------------------------------------------------------------------------------
  void StationAna::beginJob() {
  }

//-----------------------------------------------------------------------------
// make <nhits>_vs_evt for each link
//----------------------------------------------------------------------------- 
  void StationAna::endJob() {
    printf("[mu2e::StationAna] pointer to the module: 0x%8p\n",(void*) this);

  }

//-----------------------------------------------------------------------------
  void StationAna::fill_straw_histograms(StrawHist_t* Hist, StrawData_t* Data) {

    int nh = Data->list_of_hits.size();
    Hist->nsht->Fill(nh);
    for (int ih=0; ih<nh; ++ih) {
      const StrawHit* sh = Data->list_of_hits.at(ih);
      Hist->tcal->Fill(sh->time(mu2e::StrawEnd::cal));
      Hist->dtch->Fill(sh->dt());       // cal-hv
      Hist->edep->Fill(sh->energyDep());
    }

    int nhg = Data->list_of_good_hits.size();
    Hist->nshg->Fill(nh);
    for (int ih=0; ih<nhg; ++ih) {
      const StrawHit* sh = Data->list_of_good_hits.at(ih);
      Hist->tcalg->Fill(sh->time(mu2e::StrawEnd::cal));
      Hist->dtchg->Fill(sh->dt());       // cal-hv
      Hist->edepg->Fill(sh->energyDep());
    }
  }



//-----------------------------------------------------------------------------
// to handle more than one station, this function will need to evolve
//-----------------------------------------------------------------------------
  void StationAna::fill_panel_histograms(PanelHist_t* Hist, PanelData_t* Data) {

    Hist->nsht->Fill(Data->nsht);
    Hist->nshg->Fill(Data->nshg);

    for (int is=0; is<96; ++is) {
      StrawData_t* sd = &Data->straw_data[is];
      int nh = sd->list_of_hits.size();
      for (int ih=0; ih<nh; ++ih) {
        const StrawHit* sh = sd->list_of_hits.at(ih);
        Hist->tcal->Fill(sh->time(mu2e::StrawEnd::cal));
        Hist->dtch->Fill(sh->dt());       // cal-hv
        Hist->edep->Fill(sh->energyDep());
        Hist->occup->Fill(sh->strawId().straw());
      }

      int nhg = sd->list_of_good_hits.size();
      for (int ih=0; ih<nhg; ++ih) {
        const StrawHit* sh = sd->list_of_good_hits.at(ih);
        Hist->tcalg->Fill(sh->time(mu2e::StrawEnd::cal));
        Hist->dtchg->Fill(sh->dt());       // cal-hv
        Hist->edepg->Fill(sh->energyDep());
        Hist->occupg->Fill(sh->strawId().straw());
      }
      
      fill_straw_histograms(&Hist->straw[is],sd);
    }
  }
  
//-----------------------------------------------------------------------------
  void StationAna::fill_event_histograms(EventHist_t* Hist, EventData_t* Data) {
    Hist->evt->Fill(Data->evt_number);
    Hist->nsht->Fill(Data->nsht);
    Hist->nshg->Fill(Data->nshg);
    Hist->nshdt->Fill(Data->nshdt);
  }

//-----------------------------------------------------------------------------
// fill_roc_histograms also fills the channel histograms
// if in error, only histogram the error code
//-----------------------------------------------------------------------------
  int StationAna::fill_histograms() {

//-----------------------------------------------------------------------------
// later, all error handling will move to analyze_roc_data()
//-----------------------------------------------------------------------------
    fill_event_histograms(_hist.event[0]    ,&_edata);

    for (int ip=0; ip<12; ++ip) {
      fill_panel_histograms(&_hist.panel_set[0].panel[ip],&_edata.panel[ip]);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  void mu2e::StationAna::print_(const std::string& Message, const std::source_location& location) {
    if (_edata.event) {
      std::cout << std::format(" event:{}:{}:{}",
                               _edata.event->run(),_edata.event->subRun(),_edata.event->event());
    }
    std::cout << " " << location.file_name() << ":" << location.line()
      //            << location.function_name()
              << ": " << Message << std::endl;
  }

 //-----------------------------------------------------------------------------
int mu2e::StationAna::getData(const art::Event& ArtEvent) {
  int rc(0);

//-----------------------------------------------------------------------------
// tracker
//-----------------------------------------------------------------------------
//  art::Handle<mu2e::StrawDigiCollection>            sdch;
  // art::Handle<mu2e::StrawDigiADCWaveformCollection> sdawfch;
  art::Handle<mu2e::StrawHitCollection>             shch;

  // _nstrawdigis   = 0;
  // _nstrawhits    = 0;
  // _ntimeclusters = 0;
  // _ncalodigis  = 0;
  // _ncrvdigis   = 0;
  // _nstmdigis   = 0;

  // _sdc         = nullptr;
  // _sdawfc      = nullptr;
  _shc         = nullptr;
  
  bool ok = ArtEvent.getByLabel(_shCollTag,shch);
  if (ok) { 
    _shc          = shch.product();
    _n_straw_hits = _shc->size();
  }
  else {
    print_(std::format("ERROR: StrawHitCollection:{:s} is not available. Bail out",
                       _shCollTag.encode().data()));
    return -1;
  }

  // ok = ArtEvent.getByLabel(_sdCollTag,sdch);
  // if (ok) { 
  //   _sdc         = sdch.product();
  //   _nstrawdigis = _sdc->size();
  // }
  // else {
  //   print_(std::format("ERROR: StrawDigiCollection:{:s} is not available. Bail out\n",
  //                      _sdCollTag.encode().data()));
  //   return -1;
  // }

  // ok =  ArtEvent.getByLabel(_sdCollTag,sdawfch);
  // if (ok) { 
  //   _sdawfc = sdawfch.product();
  // }
  // else {
  //   print_(std::format("WARNING: StrawDigiADCWaveformCollection:{:s} is not available. Bail out\n",
  //                      _sdCollTag.encode().data()));
  //   return -1;
  // }

  // art::Handle<mu2e::TimeClusterCollection>             tcch;
  // ok =  ArtEvent.getByLabel(_tcCollTag,tcch);
  // if (ok) { 
  //   _tcc           = tcch.product();
  //   _ntimeclusters = _tcc->size();
  // }
  // else {
  //   print_(std::format("WARNING: TimeClusterCollection:{:s} is not available. Bail out\n",
  //                      _tcCollTag.encode().data()));
  //   return -1;
  // }

  return rc;
}

  //-----------------------------------------------------------------------------
int StationAna::init_event(const art::Event& ArtEvent) {
  _edata.event      = &ArtEvent;
  _edata.evt_number = ArtEvent.event();
  _edata.run_number = ArtEvent.run();
  _edata.srn_number = ArtEvent.subRun();

  _edata.nsht       = 0;
  _edata.nshg       = 0;
  _edata.nshdt      = 0;

  for (int ip=0; ip<12; ip++) {
    PanelData_t* pd = &_edata.panel[ip];
    pd->nsht      = 0;
    pd->nshg      = 0;
    for (int is=0; is<96; ++is) {
      StrawData_t* sd = &pd->straw_data[is];
      sd->list_of_hits.clear();
      sd->list_of_good_hits.clear();
    }
  }

  getData(ArtEvent);

  return 0;
}

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void StationAna::analyze(const art::Event& ArtEvent) {

  init_event(ArtEvent);

  if (_debugMode > 0) printf(" Event : %06i:%08i%08i\n", ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event());


  for (int i=0; i<_n_straw_hits; i++) {
    const mu2e::StrawHit* sh = &_shc->at(i);
    int pln = sh->strawId().plane();
    int pnl = sh->strawId().panel();
    int is  = sh->strawId().straw();
    
    const TrkPanelMap_t* tpm = _panel_map[pln][pnl];

    int pcie_addr = tpm->dtc % 2;  // this is a convention
    int ipanel    = pcie_addr*6+tpm->link;
    
    _edata.nsht               += 1;
    _edata.panel[ipanel].nsht += 1;
    std::vector<const mu2e::StrawHit*>* x = &_edata.panel[ipanel].straw_data[is].list_of_hits;
    x->push_back(sh);

    // good ("track") hits
    if (fabs(sh->dt()) < _maxDt) {
      _edata.nshdt += 1;
      if (sh->energyDep() > _minEDep) {
        _edata.nshg               += 1;
        _edata.panel[ipanel].nshg += 1;
        _edata.panel[ipanel].straw_data[is].list_of_good_hits.push_back(sh);
      }
    }
  }
  
//-----------------------------------------------------------------------------
// print debug information
//-----------------------------------------------------------------------------
  debug(ArtEvent);
  
//-----------------------------------------------------------------------------
// event data un(re)packed , fill histograms
//-----------------------------------------------------------------------------
  fill_histograms();
//-----------------------------------------------------------------------------
// finally, if requested, go into interactive mode, 
// fInteractiveMode = 0 : do not stop (default)
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
  TModule::analyze(ArtEvent);
}

// //--------------------------------------------------------------------------------
// // assume that we only have tracker fragment(s)
// //-----------------------------------------------------------------------------
// void StationAna::print_message(const char* Message) {
//   printf("StationAna: event %6i:%8i%8i %s",
//          _edata._event->run(),
//          _edata._event->subRun(),
//          _edata._event->event(),
//          Message);
// }


//-----------------------------------------------------------------------------
void StationAna::debug(const art::Event& AnEvent) {
  
  // auto handle = AnEvent.getValidHandle<std::vector<artdaq::Fragment> >(_trkfCollTag);

  // int ifrag = 0;
  // for (const artdaq::Fragment& frag : *handle) {
  //   ushort* buf          = (ushort*) (frag.dataBegin());
  //   int fsize            = frag.sizeBytes();
  //   SubEventHeader_t* sh = (SubEventHeader_t*) buf;
  //   int nbytes           = buf[0];
  //   int dtc_index        = dtcIndex(sh->dtcID);
    
  //   if (DebugBit(0) & 0x1) {
  //     print_message(Form("bit:000: fragment# %2i DTC_ID:%02i dtc_index:%i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i\n",
  //                        ifrag,sh->dtcID,dtc_index,nbytes,fsize,_edata.error_code,_edata.nerr_tot));
  //     if (DebugBit(0) & 0x2) print_fragment(&frag,nbytes/2);
  //   }

  //   if ((DebugBit(3) & 0x1) and (_edata.nerr_tot > _minNErrors)) {
  //     print_message(Form("bit:003: fragment# %2i DTC_ID:%02i dtc_index:%i nnbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i\n",
  //                        ifrag,sh->dtcID,dtc_index,nbytes,fsize,_edata.error_code,_edata.nerr_tot));

  //     if (DebugBit(3) & 0x2) print_fragment(&frag,nbytes/2);
  //   }

  //   ifrag++;
  // }

  // if ((DebugBit(2) == 1) and (_edata.error_code == _errorCode)) {
  //   print_message(Form("bit:002: ERROR_CODE: 0x%04x NERR_TOT: %5i\n",
  //                      _edata.error_code,_edata.nerr_tot));
  // }

}


  
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::StationAna)
