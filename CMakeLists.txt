cet_build_plugin(CaloDigiAnalyzerExample art::module
    REG_SOURCE CaloDigiAnalyzerExample_module.cc
    LIBRARIES REG
     artdaq_core_mu2e::artdaq-core-mu2e_Data
     artdaq_core_mu2e::artdaq-core-mu2e_Data_dict
     artdaq_core_mu2e::artdaq-core-mu2e_Overlays
     artdaq::DAQdata
     art_root_io::TFileService_service
     Offline::DAQ
     ROOT::Hist
     ROOT::Tree
)
