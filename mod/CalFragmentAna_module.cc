
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
// bit 00: print all events
// bit 02: events with given error code
// bit 03: 
//         0x01: events with the total number of errors > _minNErrors
//         0x10: total event dump
// bit 04: print problematic hits
///////////////////////////////////////////////////////////////////////////////
#include "TRACE/tracemf.h"
#define TRACE_NAME "CalFragmentAna"

#include "daqana/mod/CalFragmentAna_module.hh"

namespace mu2e {

  // ======================================================================

  // Constructor implementation
  CalFragmentAna::CalFragmentAna(fhicl::ParameterSet const& PSet) : 
    THistModule            (PSet,PSet.get<fhicl::ParameterSet>("THistModule","{}"),"CalFragmentAna") ,
    _diagLevel             (PSet.get<int>             ("diagLevel"             )), 
    _minNBytes             (PSet.get<int>             ("minNBytes"             )), 
    _maxNBytes             (PSet.get<int>             ("maxNBytes"             )), 
    _dataHeaderOffset      (PSet.get<int>             ("dataHeaderOffset"      )),
    _activeLinks_0         (PSet.get<std::vector<int>>("activeLinks_0"         )),
    _activeLinks_1         (PSet.get<std::vector<int>>("activeLinks_1"         )),
    _calCollTag            (PSet.get<std::string>     ("calCollTag"           )),
    _analyzeFragments      (PSet.get<int>             ("analyzeFragments"      )),
    // _maxFragmentSize       (PSet.get<int>             ("maxFragmentSize"       )),
    _minNErrors            (PSet.get<int>             ("minNErrors"            )),
    _errorCode             (PSet.get<int>             ("errorCode"             )),
    _validateAdcPatterns   (PSet.get<int>             ("validateAdcPatterns"   )),
    _fillHistograms        (PSet.get<int>             ("fillHistograms"        )),
    _rocDataFormat         (PSet.get<int>             ("rocDataFormat"         )),
    _initialized           (0)
  {
    _timeWindow = PSet.get<int>("timeWindow")*25.;  // in ns

    //-----------------------------------------------------------------------------
    // initialize reference channels, at this point use channels 91 and 94 for all 
    // ROC's (the readout order is defined in firmware and is the same for all channels)
    //-----------------------------------------------------------------------------
    for (int ist=0; ist<1; ist++) {
      for (int ilink=0; ilink<kMaxNLinks; ilink++) {
	RocData_t* rd = &_edata.caldtc[ist].roc[ilink];
	rd->link      = ilink;
      }
    }

    _initialized = 0;
  }

  //-----------------------------------------------------------------------------
  void CalFragmentAna::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
    Hist->nbtot           = Dir->make<TH1F>("nbtot", Form("run %06i: nbytes total", RunNumber), 1000, 0., 100000.);
    Hist->nhits           = Dir->make<TH1F>("nhits", Form("run %06i: nhits total", RunNumber), 1000, 0., 5000.);
    Hist->nfrag           = Dir->make<TH1F>("nfrag", Form("run %06i: n fragments", RunNumber), 100, 0., 100.);
    Hist->fsize           = Dir->make<TH1F>("fsize", Form("run %06i: fragment size", RunNumber), 1000, 0., 100000.);

    Hist->n_empty         = Dir->make<TH1F>("n_empty", Form("run %06i: N(empty)", RunNumber), 100, 0., 100.);
    Hist->n_invalid_dr    = Dir->make<TH1F>("n_invalid_dr", Form("run %06i: N(invalid dr)", RunNumber), 100, 0., 100.);
    Hist->n_corrupt       = Dir->make<TH1F>("n_corrupt", Form("run %06i: N(corrupt)", RunNumber), 100, 0., 100.);
    Hist->n_timeouts      = Dir->make<TH1F>("n_timeouts", Form("run %06i: N(timeouts)", RunNumber), 100, 0., 100.);
    Hist->n_overflows     = Dir->make<TH1F>("n_overflows", Form("run %06i: N(overflows)", RunNumber), 100, 0., 100.);

    Hist->error_code      = Dir->make<TH1F>("error_code", Form("run %06i: error code", RunNumber), 512, 0., 512.);
    Hist->nerr_tot        = Dir->make<TH1F>("nerr_tot", Form("run %06i: N errors", RunNumber), 1000, 0., 2000.);
    Hist->valid           = Dir->make<TH1F>("valid", Form("run %06i: valid code", RunNumber), 100, 0., 100.);
    Hist->eflg_vs_evt     = Dir->make<TH1F>("eflg_vs_evt", Form("run %06i: eflag vs evt", RunNumber), 1000, 0., 5000000.);


    Hist->nerr_vs_evt = Dir->make<TH1F>("nerr_vs_evt", Form("run %06i: N err vs evt", RunNumber), 1000, 0., 5000000.);

    // Initialize fragment_size_per_link histograms
    Hist->fragment_size_per_link.resize(kMaxNLinks);
    for (int link = 0; link < kMaxNLinks; ++link) {
      std::string hist_name = Form("fragment_size_link_%d", link);
      std::string hist_title = Form("Run %06i: Fragment Size for Link %d;Size;Counts", RunNumber, link);
      Hist->fragment_size_per_link[link] = Dir->make<TH1F>(hist_name.c_str(), hist_title.c_str(), 1000, 0., 100000.);
    }

    Hist->hits_per_channel = Dir->make<TH1F>(
					     "hits_per_channel",
					     Form("run %06i: Hits per Channel;Channel;Counts", RunNumber),
					     20, 0, 20 // 20 bins for channels 0 to 19
					     );
    for (int i = 0; i < 20; ++i) {
      Hist->hits_per_channel->GetXaxis()->SetBinLabel(i + 1, Form("Ch %d", i));
    }
 
    Hist->type_dist = Dir->make<TH1F>("type_dist", Form("run %06i: Type distribution;Type;Counts", RunNumber), 4, 0, 4);
    Hist->type_dist->GetXaxis()->SetBinLabel(1, "CAL"); // Label for bin 1
    Hist->type_dist->GetXaxis()->SetBinLabel(2, "CAP"); // Label for bin 2
    Hist->type_dist->GetXaxis()->SetBinLabel(3, "TRA"); // Label for bin 3
    Hist->type_dist->GetXaxis()->SetBinLabel(4, "LAS"); // Label for bin 4

    
    Hist->id_dist = Dir->make<TH1F>("id_dist", Form("run %06i: Board ID distribution;Board ID;Counts", RunNumber), 1000, 0, 1000);

    Hist->channel_dist_per_board.reserve(256);
    for (int board = 0; board < 256; ++board) {
      std::string hist_name = Form("channel_dist_board_%03d", board);
      std::string hist_title = Form("Run %06i: Channel Distribution for Board %03d;Channel;Counts", RunNumber, board);
      Hist->channel_dist_per_board.emplace_back(Dir->make<TH1F>(hist_name.c_str(), hist_title.c_str(), 20, 0, 20));
    }


    Hist->adc_values = Dir->make<TH1F>(
				       "adc_values",
				       Form("Run %06i: ADC Values;ADC;Counts", RunNumber),
				       10000, 0., 10000.     );

    Hist->hitn_adc_value = Dir->make<TH1F>(
					   "hitn_adc_value",
					   Form("Run %06i: Hit #n ADC Value from ROC #n;ADC Value;Counts", RunNumber),
					   20, 0., 1000.   );
    
  }



  //-----------------------------------------------------------------------------
  void CalFragmentAna::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
  
    TH1::AddDirectory(kFALSE);

    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;		// all events
    // book_event_histset[ 1] = 1;	        // events with the error code = 0

    char folder_name[100];
    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
	sprintf(folder_name,"evt_%i",i);
	art::TFileDirectory top_dir = tfs->mkdir(folder_name);

	_hist.event[i] = new EventHist_t;
	book_event_histograms(&top_dir,RunNumber,_hist.event[i]);
      }
    }
    printf("[mu2e::CalFragmentAna] pointer to the module: 0x%8p\n",(void*) this);
  }

  //-----------------------------------------------------------------------------
  void CalFragmentAna::beginRun(const art::Run& aRun) {
    int rn  = aRun.run();

    if (_initialized != 0) return;
    _initialized = 1;
    //-----------------------------------------------------------------------------
    // as a last step, book histograms - need to know the number of active links
    //-----------------------------------------------------------------------------
    if (_fillHistograms > 0) {
      book_histograms(rn);
    }
  }

  //--------------------------------------------------------------------------------
  void CalFragmentAna::beginJob() {
  }

  //-----------------------------------------------------------------------------
  void CalFragmentAna::endJob() {
    for (auto& histPair : _hitHistograms) {
      histPair.second->Write(); // Save each histogram to the output file
    }
    printf("[mu2e::CalFragmentAna] All histograms saved. Module pointer: 0x%8p\n", (void*)this);
  }


  //-----------------------------------------------------------------------------
  void CalFragmentAna::fill_event_histograms(EventHist_t* Hist, EventData_t* Ed) {

    // Fill general histograms
    Hist->error_code->Fill(Ed->error_code);
    Hist->nerr_tot->Fill(Ed->nerr_tot);
    Hist->valid->Fill(Ed->valid);

    Hist->nbtot->Fill(Ed->nbtot);
    Hist->nhits->Fill(Ed->nhtot);
    Hist->nfrag->Fill(Ed->nfrag);

    Hist->n_empty->Fill(Ed->n_empty);
    Hist->n_invalid_dr->Fill(Ed->n_invalid_dr);
    Hist->n_corrupt->Fill(Ed->n_corrupt);
    Hist->n_timeouts->Fill(Ed->n_timeouts);
    Hist->n_overflows->Fill(Ed->n_overflows);

    // Fill fragment size histograms
    for (int i = 0; i < Ed->nfrag; i++) {
      int fsize = Ed->fragments[i].nbytes;
      Hist->fsize->Fill(fsize);
    }

    // Fill error histograms vs events
    Hist->nerr_vs_evt->Fill(Ed->_event->event(), Ed->nerr_tot);
    int eflg = (Ed->nerr_tot > 0);
    Hist->eflg_vs_evt->Fill(Ed->_event->event(), eflg);

    // Loop over DTCs and ROCs
    for (int idtc = 0; idtc < kMaxDtcs; idtc++) {
      CalDtcData_t* dtcData = &Ed->caldtc[idtc];

	

      for (int link = 0; link < kMaxNLinks; link++) {
	RocData_t* rd = &dtcData->roc[link];

	if (link == 1) {
	  if (rd->hit_channel.size()) { // Ensure Hit #6 exists
	    uint16_t hitn_adc = rd->hit_channel[1];

	    // Fill the histogram for Hit #n ADC value
	    Hist->hitn_adc_value->Fill(hitn_adc);
	  }
	}

	    

	// Update fragment size per link
	if (link < static_cast<int>(Hist->fragment_size_per_link.size())) {
	  Hist->fragment_size_per_link[link]->Fill(rd->nbytes);
	} else {
	  printf("Warning: Link index %d out of range.\n", link);
	}

	// Fill ADC values histogram for this ROC
	for (const auto& adc : rd->hit_channel) {
	  Hist->adc_values->Fill(adc);
	}

	// For each hit in this ROC
	for (int ihit = 0; ihit < rd->nhits; ihit++) {
	  int chan = rd->hit_channel[ihit];

	  // Fill hits per channel histogram
	  Hist->hits_per_channel->Fill(chan);

	  // Fill channel distribution per board
	  int board_id = rd->id;
	  if (board_id >= 0 && board_id < static_cast<int>(Hist->channel_dist_per_board.size())) {
	    Hist->channel_dist_per_board[board_id]->Fill(chan);
	  } else {
	    printf("Warning: Board ID %d out of range (0-255). Channel %d not filled.\n", board_id, chan);
	  }
	}
      }
    }
  }



  //-----------------------------------------------------------------------------
  int CalFragmentAna::fill_histograms() {
    // Fill histograms for the first event hist set (all events)
    if (_fillHistograms > 0) {
      fill_event_histograms(_hist.event[0], &_edata);
    }
  
    return 0;
  }

  //-----------------------------------------------------------------------------
  int CalFragmentAna::init_event(const art::Event& AnEvent) {
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
    for (int ist=0; ist<kMaxDtcs; ist++) {
      CalDtcData_t* sd = &_edata.caldtc[ist];
    
      sd->n_empty       = 0;
      sd->n_invalid_dr  = 0;
      sd->n_corrupt     = 0;
      sd->n_timeouts    = 0;
      sd->n_overflows   = 0;
    
      for (int lnk=0; lnk<6; lnk++) {
	RocData_t* rd    = &sd->roc[lnk];
	rd->n_empty      = 0;
	rd->n_invalid_dr = 0;
	rd->n_corrupt    = 0;
	rd->n_timeouts   = 0;
	rd->n_overflows  = 0;
	rd->nerr_tot     = 0;
	rd->error_code   = 0;
	// for (int ich=0; ich<kNChannels; ich++) {
	// ChannelData_t* chd = &rd->channel[ich];
	// chd->hit.clear();
	// chd->wp.clear();
	// }
      }
    }
    return 0;
  }



  bool isHitEnd(uint16_t word1, uint16_t word2) {
    bool firstWordPattern = (word1 & 0xF000) == 0xF000;
    
    bool secondWordPattern = (word2 & 0x00FF) == 0x00FF;
    
    return firstWordPattern && secondWordPattern;
  }

  uint16_t adc_sample(uint16_t* data, int index) {
    int first_bit = 12 * index;
    int loc = first_bit / 16;

    // Combine bits from current and next word if necessary
    uint32_t combined = (static_cast<uint32_t>(data[loc]) | 
			 (static_cast<uint32_t>(data[loc + 1]) << 16));
    uint16_t result = (combined >> (first_bit % 16)) & 0xFFF; // Extract the 12-bit ADC value

    return result;
  }


  
  //-----------------------------------------------------------------------------
  void CalFragmentAna::analyze_dtc_fragment(const art::Event& Evt, const artdaq::Fragment* Fragment) {

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
    // if (fdata[0] > _maxFragmentSize) {
    //   _edata.error_code |= kNBytesErrorBit;
    //   _edata.nerr_tot   += 1;
    //   if (DebugBit(0) != 0) { 
    //     TLOG(TLVL_ERROR) << Form("event %6i:%8i:%8i : ERROR_CODE:0x%04x in %s: fdt->nbytes= %i, BAIL OUT\n",
    //                              Evt.run(),Evt.subRun(),Evt.event(),_edata.error_code,__func__,fdt->nbytes);
    //   }
    // }
    //-----------------------------------------------------------------------------
    // somewhere here handle the DTC ID
    //-----------------------------------------------------------------------------
    _node          = 0;
    //-----------------------------------------------------------------------------
    // start handling the ROC data
    // need mapping of the DTC ID to the plane number
    // for now, just label the DTCs by their PCIE address
    // sh->dtcID uniquely identifies the DTC, so it could be any number
    //-----------------------------------------------------------------------------
    SubEventHeader_t* sh = (SubEventHeader_t*) fdata;
    int dtc_index        = 0; // dtcIndex(sh->dtcID);
  
    short* first_address = fdata + sizeof(SubEventHeader_t)/2; // _dataHeaderOffset; // offset is specified in 2-byte words
    short* last_address  = fdata + fdt->nbytes/2;     // 

    while (first_address < last_address) {
      //-----------------------------------------------------------------------------
      // next ROC
      //-----------------------------------------------------------------------------
      RocDataHeaderPacket_t* dh = (RocDataHeaderPacket_t*) first_address;
      int link      = dh->linkID;
      RocData_t* rd = &_edata.caldtc[_node].roc[link];
      rd->dtc_id    = sh->dtcID;
      //-----------------------------------------------------------------------------
      // check link number
      //-----------------------------------------------------------------------------
      int found = 0;
      for (int i=0; i<_nActiveLinks[dtc_index]; i++) {
	if (_activeLinks[dtc_index]->at(i) == link) {
	  found = 1;
	  break;
	}
      }
    
      if (found == 0) {
	//-----------------------------------------------------------------------------
	// in some cases, this one could be just a warning
	//-----------------------------------------------------------------------------
	//        _edata.error_code |= kLinkIDErrorBit;
	_edata.nerr_tot   += 1;
	rd->nerr_tot      += 1;
	// if (DebugBit(0) != 0) {
	//   TLOG(TLVL_ERROR) << Form("event %6i:%8i:%8i : ERROR_CODE:0x%04x in %s: link=%i not defined as active, BAIL OUT\n",
	//                            Evt.run(),Evt.subRun(),Evt.event(),kLinkIDErrorBit,__func__,link);
	// }
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

      analyze_roc_data(dh, rd);
      //-----------------------------------------------------------------------------
      // update station error counters
      //-----------------------------------------------------------------------------
      // sd->n_empty      += rd->n_empty; 
      // sd->n_invalid_dr += rd->n_invalid_dr; 
      // sd->n_corrupt    += rd->n_corrupt; 
      // sd->n_timeouts   += rd->n_timeouts; 
      // sd->n_overflows  += rd->n_overflows; 
      //-----------------------------------------------------------------------------
      // address in 2-byte words (N(data packets)+data header packet)
      //-----------------------------------------------------------------------------
      first_address    += (dh->packetCount + 1)*8;
    }
  }


  //--------------------------------------------------------------------------------
  void CalFragmentAna::analyze(const art::Event& AnEvent) {

    init_event(AnEvent);

    auto handle = AnEvent.getValidHandle<std::vector<artdaq::Fragment> >(_calCollTag);
    //-----------------------------------------------------------------------------
    // proxy for event histograms
    //-----------------------------------------------------------------------------
    if (_diagLevel > 0) {
      printf(" Event : %06i:%08i%08i\n", AnEvent.run(),AnEvent.subRun(),AnEvent.event());
    }
  
    int ifrag = 0;
    for (const artdaq::Fragment& frag : *handle) {
      //-----------------------------------------------------------------------------
      // different fragments correspond to different DTCs, somewhere there should be the DTC ID
      //-----------------------------------------------------------------------------
      ushort* buf = (ushort*) (frag.dataBegin());
      int nbytes  = buf[0];
      int fsize   = frag.sizeBytes();

      if (nbytes < 2) {
	// _edata.error_code  |= kNBytesErrorBit;
	_edata.nerr_tot    += 1;
	_edata.n_nb_errors += 1;
      
	// TLOG(TLVL_DEBUG) << Form("event %i:%i:%i : ERROR_CODE:0x%04x nbytes=%i EMPTY_FRAGMENT",
	//                          AnEvent.run(),AnEvent.subRun(),AnEvent.event(),
	//                          kNBytesErrorBit,nbytes);
      }

      _edata.nfrag += 1;
      _edata.nbtot += nbytes;        // including the artdaq part

      if (_diagLevel > 2) {
	printf("%s\n",Form("---------- fragment # %3i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x\n",
			   ifrag,nbytes,fsize,_edata.error_code));
	print_fragment(&frag,nbytes/2);
      }

      //    if ((_edata.error == 0) and _analyzeFragments) analyze_fragment(event,&frag);
      if (_analyzeFragments) analyze_dtc_fragment(AnEvent,&frag);
      ifrag++;
    }
    //-----------------------------------------------------------------------------
    // proxy for event histograms
    //-----------------------------------------------------------------------------
    if (_diagLevel > 1) {
      if ((_edata.nbtot >= _minNBytes) and (_edata.nbtot <= _maxNBytes)) {
	TLOG(TLVL_DEBUG) << Form("AAAAAAAAAAAAAAAAAAAAAAAAAAAA Run : %5i subrun: %5i event: %8i nfrag: %3i nbytes: %5i\n", 
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
    // finally, if requested, go into interactive mode, 
    // fInteractiveMode = 0 : do not stop (default)
    // fInteractiveMode = 1 : stop after each event (event display mode)
    // fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
    //-----------------------------------------------------------------------------
    TModule::analyze(AnEvent);
  }


  //--------------------------------------------------------------------------------
  void CalFragmentAna::print_message(const char* Message) {
    printf("OOOOOOOOOOOOOOOOOOOOO CalFragmentAna: event %6i:%8i%8i %s\n",
	   _edata._event->run(),
	   _edata._event->subRun(),
	   _edata._event->event(),
	   Message);
  }

  //-----------------------------------------------------------------------------
  void CalFragmentAna::print_fragment(const artdaq::Fragment* Frag, int NWords) {
    //-----------------------------------------------------------------------------
    // print fragments in HEX, for the tracker, the data has to be in 2-byte words
    //-----------------------------------------------------------------------------
    ushort* buf = (ushort*) (Frag->dataBegin());

    int loc     = 0;
    
    for (int i=0; i<NWords; i++) {
      if (loc == 0) printf("0x%08x: ",i*2);

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
  void CalFragmentAna::debug(const art::Event& AnEvent) {
  
    auto handle = AnEvent.getValidHandle<std::vector<artdaq::Fragment> >(_calCollTag);

    int ifrag = 0;
    for (const artdaq::Fragment& frag : *handle) {
      ushort* buf          = (ushort*) (frag.dataBegin());
      int fsize            = frag.sizeBytes();
      //    SubEventHeader_t* sh = (SubEventHeader_t*) buf;
      int nbytes           = buf[0];
      int dtc_index        = 0; // dtcIndex(sh->dtcID);
    
      if (DebugBit(0) == 1) {
	print_message(Form("bit:000: fragment # %3i dtc_index:%i nbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i",
			   ifrag,dtc_index,nbytes,fsize,_edata.error_code,_edata.nerr_tot));
	print_fragment(&frag,nbytes/2);
      }

      if ((DebugBit(3) & 0x1) and (_edata.nerr_tot > _minNErrors)) {
	print_message(Form("bit:003: fragment # %3i dtc_index:%i nnbytes: %5i fsize: %5i ERROR_CODE: 0x%04x NERR_TOT: %5i",
			   ifrag,dtc_index,nbytes,fsize,_edata.error_code,_edata.nerr_tot));

	if (DebugBit(3) & 0x2) print_fragment(&frag,nbytes/2);
      }
      ifrag++;
    }

    if ((DebugBit(2) == 1) and (_edata.error_code == _errorCode)) {
      print_message(Form("bit:002: ERROR_CODE: 0x%04x NERR_TOT: %5i",
			 _edata.error_code,_edata.nerr_tot));
    }

  }

  void CalFragmentAna::analyze_roc_data(RocDataHeaderPacket_t* Dh, RocData_t* Rd) {
    // Initialize ROC Data
    Rd->nbytes = Dh->byteCount;
    Rd->npackets = Dh->packetCount;
    Rd->nhits = 0;
    Rd->valid = Dh->valid;
    Rd->dtc_id = Dh->dtcID;
    Rd->link = Dh->linkID;

    // Update ROC status counters
    if (Dh->empty()) Rd->n_empty += 1;
    if (Dh->invalid_dr()) Rd->n_invalid_dr += 1;
    if (Dh->corrupt()) Rd->n_corrupt += 1;
    if (Dh->timeout()) Rd->n_timeouts += 1;
    if (Dh->overflow()) Rd->n_overflows += 1;

    // Store error code
    Rd->error_code = Dh->error_code();

    int hitIndex = 0; // To keep track of the number of hits

    // Parse header words
    uint16_t* headerPtr = reinterpret_cast<uint16_t*>(Dh);
    uint16_t w0 = headerPtr[0];
    uint16_t w1 = headerPtr[1];
    uint16_t w2 = headerPtr[2];

    // Extract type, id, and channel from the header
    unsigned type = (w0 >> 13) & 0x7;
    unsigned id = (w0 >> 5) & 0xFF;
    unsigned channel = w0 & 0x1F;

    // Store type and id in ROC data
    Rd->type = type;
    Rd->id = id;

    // Fill histograms for type and id
    _hist.event[0]->type_dist->Fill(type);
    _hist.event[0]->id_dist->Fill(id);

    // Combine w1 and w2 to form a 32-bit reserved field
    uint32_t empty32 = (static_cast<uint32_t>(w1) << 16) | w2;
    printf("\n--- PARSED HEADER BITS ---\n");
    printf(" Detector Type    (3 bits)  = 0x%x  (%u decimal)\n", type, type);
    printf(" Board ID         (8 bits)  = 0x%x  (%u decimal)\n", id, id);
    printf(" Channel Number   (5 bits)  = 0x%x  (%u decimal)\n", channel, channel);
    printf(" Reserved         (32 bits) = 0x%08x (%u decimal)\n", empty32, empty32);

    // Pointer to the start of data after the header
    uint16_t* dataPtr = reinterpret_cast<uint16_t*>(Dh + 1);
    int remainingWords = Rd->npackets * 8;

    printf("\n[CalFragmentAna::analyze_roc_data] Parsing hits dynamically...\n");

    bool hasHits = false; // Flag to indicate if any hits are found

    while (remainingWords > 0) {
      uint16_t word1 = dataPtr[0];

      // Detect single-word termination pattern
      if ((word1 & 0x0FFF) == 0x0FFF) {
	printf("Single-word termination pattern detected: 0x%04x\n", word1);
	dataPtr++;
	remainingWords--;
	continue;
      }

      // Detect two-word termination pattern
      if (remainingWords >= 2 && isHitEnd(dataPtr[0], dataPtr[1])) {
	printf("Hit termination pattern detected at: 0x%04x 0x%04x\n", dataPtr[0], dataPtr[1]);
	dataPtr += 2;
	remainingWords -= 2;
	continue;
      }

      // Check if remaining words could be a footer (last 10-11 words with no hits)
      if (remainingWords <= 11) {
	printf("[CalFragmentAna::analyze_roc_data] Footer detected. Parsing footer...\n");

	// Extract footer information
	Rd->footerData.clear();
	for (int i = 0; i < remainingWords; ++i) {
	  Rd->footerData.push_back(dataPtr[i]);
	}

	// Debug footer content
	printf("Footer Words:\n");
	for (size_t i = 0; i < Rd->footerData.size(); ++i) {
	  printf(" 0x%04x", Rd->footerData[i]);
	  if ((i + 1) % 8 == 0 || i == Rd->footerData.size() - 1) {
	    printf("\n");
	  }
	}

	// Break out of the loop since the footer is the last part
	break;
      }

      // Parse hits as usual
      std::vector<uint16_t> currentHit;

      while (remainingWords > 0) {
	// Check for two-word termination pattern within hit
	if (remainingWords >= 2 && isHitEnd(dataPtr[0], dataPtr[1])) {
	  break;
	}

	// Check for single-word termination pattern within hit
	if ((dataPtr[0] & 0x0FFF) == 0x0FFF) {
	  printf("Single-word termination pattern within hit detected: 0x%04x\n", dataPtr[0]);
	  dataPtr++;
	  remainingWords--;
	  break;
	}

	// Add word to current hit
	currentHit.push_back(*dataPtr);
	if (_diagLevel > 2) { // Conditional logging based on diagnostic level
	  printf("Added Word to Hit: 0x%04x\n", *dataPtr);
	}
	dataPtr++;
	remainingWords--;
      }

      if (!currentHit.empty()) {
	hasHits = true;
	Rd->nhits += 1;
	hitIndex++;
	printf("[CalFragmentAna::analyze_roc_data] Hit #%d Length: %lu Words:\n", Rd->nhits, currentHit.size());

	for (size_t i = 0; i < currentHit.size(); ++i) {
	  printf(" 0x%04x", currentHit[i]);
	  if ((i + 1) % 8 == 0 || i == currentHit.size() - 1) {
	    printf("\n");
	  }
	}

	// Create a unique histogram name using hitIndex and ROC ID
	std::string histName = Form("Hit_%d_Board_%d_adc_values", hitIndex, Rd->id);

	// Check if the histogram already exists
	if (_hitHistograms.find(histName) == _hitHistograms.end()) {
	  _hitHistograms[histName] = new TH1F(
					      histName.c_str(),
					      Form("ADC Values for Hit #%d ROC #%d;ADC Value;Count", hitIndex, Rd->id),
					      100, 0., 100. // ADC values range
					      );
	}

	// Fill the histogram with ADC values from the current hit
	for (size_t i = 0; i < currentHit.size(); ++i) {
	  uint16_t adcValue = adc_sample(currentHit.data(), i);
	  // printf(" 0x%04x", *((uint16_t*)currentHit.data()+i)& 0x0FFF);
	  _hitHistograms[histName]->SetBinContent(i+1,adcValue);
	  _hitHistograms[histName]->SetBinError(i+1,0);
	  
	  
	}

	// Correctly associate the channel number from the header
	if (channel < 32) { // Ensure channel number is within expected range (5 bits)
	  Rd->hit_channel.push_back(channel); // Use channel from header
	  _hist.event[0]->hits_per_channel->Fill(channel);

	  int boardId = Rd->id;
	  if (boardId >= 0 && boardId < 256) {
	    _hist.event[0]->channel_dist_per_board[boardId]->Fill(channel);
	  } else {
	    printf("Warning: Board ID %d out of range (0-255). Channel %d not filled.\n", boardId, channel);
	  }
	} else {
	  printf("Warning: Channel number %d out of expected range (0-31).\n", channel);
	}
      }
    }

    if (!hasHits) {
      printf("[CalFragmentAna::analyze_roc_data] No hits detected for Board ID: 0x%x\n", Rd->id);
    }
  }



} // end namespace mu2e


DEFINE_ART_MODULE(mu2e::CalFragmentAna)
