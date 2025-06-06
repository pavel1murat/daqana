//-----------------------------------------------------------------------------
// c->Print("trk_fragment_ana_run_002386.figure_0001.pdf")
// c_nh_per_panel_dtc_1->Print("trk_fragment_ana_run_002386.figure_0102.pdf"
// c_nh_vs_ch_dtc_1->Print("trk_fragment_ana_run_002386.figure_0103.pdf")
// c_nh_vs_adc_0_dtc_1->Print("trk_fragment_ana_run_002386.figure_0104.pdf")
//-----------------------------------------------------------------------------


int MaxNHits   (-1);
int MaxChannel (-1);
int MaxEvent   (-1);
int MaxY       (-1);

namespace {
  //  const char* HistDir   = "/exp/mu2e/data/projects/tracker/vst/hist" ;
  const char* HistDir   = "/data/tracker/vst/hist" ;
  int         _Figure   = 0;
  TCanvas*    _Canvas(nullptr);
};

TFile* open_file(const char* Fn, int RunNumber) {
  TFile* f(nullptr);
  
  //  if (Fn == nullptr) f = TFile::Open(Form("%s/trk_fragment_ana.%06i.hist",HistDir,RunNumber));
  if (Fn == nullptr) f = TFile::Open(Form("%s/hst.mu2e.vst00s000r000n000.trk_fragment_ana.%06i_000001.root",HistDir,RunNumber));
  else               f = TFile::Open(Fn);

  return f;
}

//-----------------------------------------------------------------------------
int plot_totals(int RunNumber, const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure),Form("c_%05i",_Figure),1700,1000);
  c->Divide(3,2);

  c->cd(1);
  gPad->SetLogy(1);
  TH1* h1 = (TH1*) f->Get("//TrkFragmentAna/evt_0/nhits");
  h1->GetYaxis()->SetRangeUser(0.1,2*h1->GetEntries());
  h1->Draw();
  
  c->cd(2);
  gPad->SetLogy(1);
  TH1* h2 = (TH1*) f->Get("//TrkFragmentAna/evt_0/nbtot");
  h2->GetYaxis()->SetRangeUser(0.1,2*h2->GetEntries());
  h2->Draw();
  
  c->cd(3);
  gPad->SetLogy(1);
  TH1* h3 = (TH1*) f->Get("//TrkFragmentAna/evt_0/nfrag");
  h3->GetYaxis()->SetRangeUser(0.1,2*h3->GetEntries());
  h3->GetXaxis()->SetRangeUser(0,19.9);
  h3->Draw();
  
  c->cd(4);
  gPad->SetLogy(1);
  TH1* h4 = (TH1*) f->Get("//TrkFragmentAna/evt_0/fsize");
  h4->GetYaxis()->SetRangeUser(0.1,2*h4->GetEntries());
  h4->Draw();
  
  c->cd(5);
  gPad->SetLogy(1);
  TH1* h5 = (TH1*) f->Get("//TrkFragmentAna/evt_0/error_code");
  h5->GetYaxis()->SetRangeUser(0.1,2*h5->GetEntries());
  h5->Draw();
  
  c->cd(6);
  gPad->SetLogy(1);
  TH1* h6 = (TH1*) f->Get("//TrkFragmentAna/evt_0/nerr_tot");
  h6->GetXaxis()->SetRangeUser(0.,50);
  h6->GetYaxis()->SetRangeUser(0.1,2*h6->GetEntries());
  h6->Draw();

  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
int plot_nhits_per_panel(int RunNumber, int Station, int Dtc, const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/nhits",Station,Dtc,i));
    h->Draw();
  }

  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
int plot_nh_vs_ch_per_panel(int RunNumber, int Station, int Dtc, int MaxChannel = -1, int MaxNHits = -1,
                              const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(0);

    TH2* h = (TH2*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/nh_vs_ch",Station,Dtc,i));

    if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->SetFillColor(kBlue+2);
    h->Draw("box");
  }
  
  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
// these are 1D histograms
//-----------------------------------------------------------------------------
int plot_nh_vs_ich_per_panel(int RunNumber, int Station, int Dtc, int MaxChannel = -1, int MaxNHits = -1,
                             const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(0);

    TH1* h = (TH2*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/nh_vs_ich",Station,Dtc,i));

    if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->Draw("");
  }
  
  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
// these are 1D histograms
//-----------------------------------------------------------------------------
int plot_ph_per_panel(int RunNumber, int Station, int Dtc, int MaxChannel = -1, int MaxNHits = -1,
                             const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH1* h = (TH2*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/ph",Station,Dtc,i));

    // if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    // if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->Draw("");
  }
  
  _Canvas = c;
  return 0;
}

int plot_nh100_per_panel(int RunNumber, int Station, int Dtc, int MaxChannel = -1, int MaxNHits = -1,
                         const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH1* h = (TH2*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/nh100",Station,Dtc,i));
    h->GetXaxis()->SetRangeUser(0,18.8);

    // if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    // if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->Draw("");
  }
  
  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
// use default readout map
//-----------------------------------------------------------------------------
int plot_nh_vs_adc_per_panel(int RunNumber, int Station, int Dtc,
                             int MaxChannel=-1, int MaxNHits = -1, const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+Dtc*100),Form("c_%05i",_Figure+Dtc*100),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(0);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/nh_vs_adc0",Station,Dtc,i));

    if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->SetFillColor(kBlue+2);
    h->Draw("box");
  }
  
  _Canvas = c;
  return 0;
}


//-----------------------------------------------------------------------------
// use default readout map
//-----------------------------------------------------------------------------
int plot_eflg_vs_evt_panels(int RunNumber, int Station, int Dtc, int MaxEvent,
                     const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+100*Dtc),Form("c_%05i",_Figure+100*Dtc),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(0);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/eflg_vs_evt",Station,Dtc,i));
                                        // dont' need the error bars
    int nbx = h->GetNbinsX();
    for (int ib=0; ib<nbx; ib++) {
      h->SetBinError(ib+1,0);
    }
    
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);

    if (MaxY     > 0) h->GetYaxis()->SetRangeUser(0,MaxY-0.0001);
    if (MaxEvent > 0) h->GetXaxis()->SetRangeUser(0,MaxEvent-0.0001);

    h->Draw("p");
  }
  
  _Canvas = c;
  return 0;
}

//-----------------------------------------------------------------------------
// use default readout map
//-----------------------------------------------------------------------------
int plot_eflg_vs_evt(int RunNumber, const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure),Form("c_%05i",_Figure),1700,950);

  gPad->SetLogy(0);

  TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/evt_0/eflg_vs_evt"));
                                        // dont' need the error bars
  int nbx = h->GetNbinsX();
  for (int ib=0; ib<nbx; ib++) {
    h->SetBinError(ib+1,0);
  }
    
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);

  if (MaxEvent > 0) h->GetXaxis()->SetRangeUser(0,MaxEvent-0.0001);

  h->Draw("p");
  _Canvas = c;
  return 0;
}



//-----------------------------------------------------------------------------
int plot_error_code(int RunNumber, int Station, int Dtc, const char* Fn = nullptr, int Print = 0) {

  TFile* f = open_file(Fn,RunNumber);

  TCanvas* c = new TCanvas(Form("c_%05i",_Figure+Dtc*100),Form("c_%05i",_Figure+Dtc*100),1700,950);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogx(1);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/stn_%02i/dtc%i/roc%i/errcode",Station,Dtc,i));
                                        // dont' need the error bars
    h->Draw();
  }

  _Canvas = c;
  return 0;
}


//-----------------------------------------------------------------------------
// if Fn = nullptr, use defaults
//-----------------------------------------------------------------------------
int plot(int RunNumber, const char* Fn, int Figure, int Station = 0, int Dtc = 0, int Print = 0) {
  int rc(0);
  
  _Figure = Figure;
  
  if      (Figure == 1) plot_totals              (RunNumber, Fn, Print);
  else if (Figure == 2) plot_nhits_per_panel     (RunNumber, Station, Dtc, Fn, Print);
  else if (Figure == 3) plot_nh_vs_ch_per_panel  (RunNumber, Station, Dtc, MaxChannel,MaxNHits,Fn, Print);
  else if (Figure == 4) plot_nh_vs_adc_per_panel (RunNumber, Station, Dtc, MaxChannel,MaxNHits,Fn, Print);
  else if (Figure == 5) plot_eflg_vs_evt_panels  (RunNumber, Station, Dtc, MaxEvent  ,Fn, Print);
  else if (Figure == 6) plot_eflg_vs_evt         (RunNumber, Fn, Print);
  else if (Figure == 7) plot_error_code          (RunNumber, Station, Dtc, Fn, Print);
  else if (Figure == 8) plot_nh_vs_ich_per_panel (RunNumber, Station, Dtc, MaxChannel,MaxNHits,Fn, Print);
  else if (Figure == 9) plot_ph_per_panel        (RunNumber, Station, Dtc, MaxChannel,MaxNHits,Fn, Print);
  else if (Figure ==10) plot_nh100_per_panel     (RunNumber, Station, Dtc, MaxChannel,MaxNHits,Fn, Print);
  else {
    printf("ERROR: undefined figure %i\n",Figure);
    return -1;
  }

  if ((rc == 0) and (Print > 0)) {
    // everything looks ok, print
    char fn[1000];
    sprintf(fn,"trk_fragment_ana_run_%06i.figure_%05i.pdf",RunNumber,Figure+100*Dtc);
    _Canvas->Print(fn);
  }

  return 0;
}
