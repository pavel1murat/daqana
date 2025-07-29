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
// -*- buffer-read-only:t -*- 
// ======================================================================
//
// ======================================================================
// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

#include "fhiclcpp/ParameterSet.h"

// #include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
// #include "artdaq-core/Data/Fragment.hh"

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>

#include "TRACE/tracemf.h"
#define TRACE_NAME "FileCatalog"

// ======================================================================
namespace mu2e {

  class FileCatalog : public THistModule {

  public:

    struct Range_t {
      int srn;
      int ev1;
    };
                                        // pointer to the raw event data
    struct EventData_t {
      const art::Event*     event;
      std::vector<Range_t>  range;
    } _edata;

//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    int              _debugMode;
    int              _debugBit[100];
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    int                             _initialized;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    explicit FileCatalog(fhicl::ParameterSet const& pset);
    // explicit FileCatalog(const art::EDAnalyzer::Table<Config>& config);
    virtual ~FileCatalog() {}
    
    virtual void analyze (const art::Event& ArtEvent) override;
    virtual void beginRun(const art::Run&   ArtRun  ) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    void         debug                  (const art::Event& event);
  
    void         print_(const std::string& Message,
                        const std::source_location& location = std::source_location::current());
  };


// ======================================================================

  FileCatalog::FileCatalog(fhicl::ParameterSet const& PSet) : 
    THistModule            (PSet,PSet.get<fhicl::ParameterSet>("THistModule"),"FileCatalog") ,
    _debugMode             (PSet.get<int>                     ("debugMode"  )              ) {
    _initialized         = 0;
  }

//-----------------------------------------------------------------------------
// for run<=285, wasn't saving the event header and subheader, 64 bytes in total
// or 32 2-byte words
// book histograms only once
//-----------------------------------------------------------------------------
  void FileCatalog::beginRun(const art::Run& ArtRun) {
    // int rn  = ArtRun.run();
    
    if (_initialized != 0) return;
    _initialized = 1;
  }
  
//--------------------------------------------------------------------------------
  void FileCatalog::beginJob() {
  }

//-----------------------------------------------------------------------------
// make <nhits>_vs_evt for each link
//----------------------------------------------------------------------------- 
  void FileCatalog::endJob() {
    printf("[mu2e::FileCatalog] pointer to the module: 0x%8p\n",(void*) this);

  }

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void FileCatalog::analyze(const art::Event& ArtEvent) {

  if (_debugMode > 0) printf(" Event : %06i:%08i:%08i\n", ArtEvent.run(),ArtEvent.subRun(),ArtEvent.event());

}

//-----------------------------------------------------------------------------
void FileCatalog::debug(const art::Event& AnEvent) {
  
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

DEFINE_ART_MODULE(mu2e::FileCatalog)
