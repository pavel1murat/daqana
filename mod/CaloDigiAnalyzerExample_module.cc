#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Principal/Run.h"

#include "messagefacility/MessageLogger/MessageLogger.h"



#include "canvas/Utilities/Exception.h"

#include "canvas/Utilities/InputTag.h"



#include "Offline/RecoDataProducts/inc/CaloDigi.hh"



#include "art_root_io/TFileService.h"



#include <iostream>

#include <iomanip>

#include <fstream> 

#include <sstream>

#include <vector>

#include <algorithm>

#include <cmath>

#include <map>



#include "TH1.h"

#include "TH2.h"

#include "TFile.h"

#include "TTree.h"

#include "TCanvas.h"



namespace mu2e {

  class CaloDigiAnalyzerExample : public art::EDAnalyzer

  {

  public:

    struct Config {

      fhicl::Atom<std::string> caloDigiModuleLabel {fhicl::Name("caloDigiModuleLabel" ) , fhicl::Comment("caloDigiModuleLabel"), ""};

    };



    explicit CaloDigiAnalyzerExample(const art::EDAnalyzer::Table<Config>& config);

    void analyze(art::Event const& event) override;

    void endJob() override;



  private:

    std::string   caloDigiModuleLabel_;

	

    TH1F* h1_example;



  };

}  // namespace mu2e



mu2e::CaloDigiAnalyzerExample::CaloDigiAnalyzerExample(const art::EDAnalyzer::Table<Config>& config)

  : art::EDAnalyzer{config},

    caloDigiModuleLabel_(config().caloDigiModuleLabel())

{

		//Something that you want to happen at the very beginning (constructor)

		art::ServiceHandle<art::TFileService> tfs;

        h1_example = tfs->make<TH1F>("hist_name","hist_title",100,0,10);

}



void mu2e::CaloDigiAnalyzerExample::analyze(art::Event const& event){

  //This will be executed once for every event



  //Get the vector of calodigis of this event

  const auto& caloDigis = *event.getValidHandle(consumes<mu2e::CaloDigiCollection>(caloDigiModuleLabel_));



  //Loop over the calo digis of this event

  for(uint ihit = 0; ihit<caloDigis.size(); ihit++){



    // int SiPMID = caloDigis[ihit].SiPMID();

    int t0 = caloDigis[ihit].t0();

    //int peakpos = caloDigis[ihit].peakpos();

    std::vector<int> waveform = caloDigis[ihit].waveform();



	  /* analysis code */

    h1_example->Fill(t0);
  }
  



}



void mu2e::CaloDigiAnalyzerExample::endJob(){

  //Something that you want to happen at the very end

}



DEFINE_ART_MODULE(mu2e::CaloDigiAnalyzerExample)
