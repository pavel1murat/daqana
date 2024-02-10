//-----------------------------------------------------------------------------
// Stntuple/mod/scripts/display_vst.C
// display tracker VST data 
//-----------------------------------------------------------------------------
namespace {
  TObjArray*          list_of_wf[6];       // to be accessed interactively
  TCanvas*            c(nullptr);
}

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

mu2e::MuHitDisplay* m_disp (nullptr); 
//-----------------------------------------------------------------------------
// Mode = 0: begin job (run)
// Mode = 1: event
// Mode = 2: end job (run)
//-----------------------------------------------------------------------------
void display_waveforms(int Mode, TModule* Module) {
  printf("display_vst called: Mode = %i, Module = %p\n",Mode,Module);

  // TStnVisManager* vm = TStnVisManager::Instance();
  
  m_disp = (mu2e::MuHitDisplay*) Module;

  if (Mode == 0) {  
//-----------------------------------------------------------------------------
// book histograms , perform other initializations
//-----------------------------------------------------------------------------

    c = new TCanvas("c","c",1500,800);
    c->Divide(3,2);

    for (int i=0; i<6; i++) {
      list_of_wf[i] = new TObjArray(16);
      list_of_wf[i]->SetOwner(kTRUE);
    }
  }
  else if (Mode == 1) {
//-----------------------------------------------------------------------------
// event entry point: fill histograms 
//-----------------------------------------------------------------------------
    const mu2e::StrawHitCollection*               shc = m_disp->GetShColl();
    const mu2e::StrawDigiCollection*              sdc = m_disp->GetSdColl();
    const mu2e::StrawDigiADCWaveformCollection*   swc = m_disp->GetSwColl();

    for (int ip=0; ip<6; ip++) list_of_wf[ip]->Delete();

    if (sdc and swc) {

      printf(" Straw DIGI waveforms \n");
      printf("  sID    t0     t1      tot0   tot1    eDep  ns  base ");
      printf(" -------------------------- samples ----------------------------");
      printf("----------------------------------------------------------------------\n");

      int nw = swc->size();
      for (int i=0; i<nw; i++) {
	const mu2e::StrawDigi*            sd = &sdc->at(i);
	const mu2e::StrawDigiADCWaveform* sw = &swc->at(i);
	const mu2e::StrawHit*             sh = &shc->at(i);

	int nsamples      = sw->samples().size();
	mu2e::StrawId sid = sd->strawId();
	int panel         = sid.panel();
	int straw         = sid.straw();

	printf("%5i %7.1f %7.1f %6.1f %6.1f %8.5f",
	       sid.asUint16(),
	       sh->time(mu2e::StrawEnd::hv),sh->time(mu2e::StrawEnd::cal),
	       sh->TOT (mu2e::StrawEnd::hv),sh->TOT (mu2e::StrawEnd::cal),
	       sh->energyDep());

	TH1F* h = new TH1F(Form("h_panel_%i_ch_%02i",panel,straw),
			   Form("panel_%i_ch_%02i"  ,panel,straw),
			   30,0,30);
//-----------------------------------------------------------------------------
// determine the baseline for this waveform
//-----------------------------------------------------------------------------
	float sum=0;
	for (int is=0; is<5; is++) {
	  sum += sw->samples().at(is);
	}

	float baseline = sum/5.;

	printf("%3i %5.1f",nsamples,baseline);
//-----------------------------------------------------------------------------
// fill histogram and print the waveform 
//-----------------------------------------------------------------------------
	for (int is=0; is<nsamples; is++) {
	  printf("%5i",sw->samples().at(is));
	  float adc = sw->samples().at(is)-baseline;
	  h->SetBinContent(is+1,adc);
	}
	printf("\n");

	h->GetYaxis()->SetRangeUser(-350,650);
//-----------------------------------------------------------------------------
// pad numbering starts from 1
//-----------------------------------------------------------------------------
	c->cd(panel+1);

	int nw_panel =  list_of_wf[panel]->GetEntriesFast();

	h->SetLineColor(nw_panel+1);

	if (nw_panel == 0) h->Draw("");
	else               h->Draw("sames");

	list_of_wf[panel]->Add(h);

	gPad->Modified();
	gPad->Update();
      }
      printf("n(waveforms):  %i\n",nw);
    }
  }
  else if (Mode == 2) {
//-----------------------------------------------------------------------------
// end run
//-----------------------------------------------------------------------------
  }
}
