//


namespace {
  //  const char* HistDir   = "/exp/mu2e/data/projects/tracker/vst/hist" ;
  const char* HistDir   = "/data/tracker/vst/hist" ;
  int         _Figure   = 0;
  TCanvas*    _Canvas(nullptr);
};

//-----------------------------------------------------------------------------
TFile* open_file(const char* Fn, int RunNumber) {
  TFile* f(nullptr);
  
  if (Fn == nullptr) f = TFile::Open(Form("%s/hst.mu2e.vst00s000r000n000.trk_fragment_ana.%06i_000001.root",
                                          HistDir,RunNumber));
  else               f = TFile::Open(Fn);

  return f;
}

//-----------------------------------------------------------------------------
int fit_bl_width_roc(int RunNumber, const char* Fn, int Dtc, int Link) {

  TFile* f = open_file(Fn,RunNumber);

  TH1F* h_sig_vs_ch = new TH1F(Form("h_sig_vs_ch_%i_%i",Dtc,Link),
                               Form("run %06i: DTC%i:LINK%i sigma(BL) vs ch",RunNumber,Dtc,Link),100,0,100);

  for (int i=0; i<96; ++i) {

    std::string hist_name=Form("/TrkFragmentAna/stn_00/dtc%i/roc%i/ch_%02i/ch_%i_%02i_bl",Dtc,Link,i,Link,i);

    std::cout << "hist_name:" << hist_name << std::endl;
    
    TH1F* h = (TH1F*) f->Get(hist_name.data());

    // h->Draw();
                           
    TFitResultPtr res = h->Fit("gaus","S");

    TFitResult* pr = res.Get();

    std::cout << "pr:" << (void*) pr << std::endl;

    double x0 (-1.), sig(-1.), sig_err(0);
  
    if (pr and pr->IsValid()) {
      x0      = res.Get()->Parameter(1);
      sig     = res.Get()->Parameter(2);
      sig_err = res.Get()->ParError (2);
    }

    std::cout << "x0:" << x0 << " sig:" << sig << " sig_err:" << sig_err << std::endl;

    h_sig_vs_ch->SetBinContent(i+1,sig);
    h_sig_vs_ch->SetBinError  (i+1,sig_err);
  }

  h_sig_vs_ch->SetMarkerStyle(20);
  h_sig_vs_ch->SetMarkerSize (0.8);
  h_sig_vs_ch->GetYaxis()->SetRangeUser(-5,20);
  h_sig_vs_ch->Draw();
  
  return 0;
}

//-----------------------------------------------------------------------------
int fit_bl_width_dtc(int RunNumber, const char* Fn, int Dtc = 0, int Print = 0) {

  TCanvas* c = new TCanvas(Form("c_dtc%i",Dtc),Form("c_dtc%i",Dtc),1800,900);
  c->Divide(3,2);

  for (int ilink=0; ilink<6; ++ilink) {
    c->cd(ilink+1);
    fit_bl_width_roc(RunNumber, Fn, Dtc, ilink);
  }

  return 0;
}

