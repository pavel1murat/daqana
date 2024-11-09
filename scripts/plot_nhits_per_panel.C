//
namespace {
  const char* HistDir = "/exp/mu2e/data/projects/tracker/vst/hist" ;
};

//-----------------------------------------------------------------------------
int plot_totals(int RunNumber) {

  TFile* f = TFile::Open(Form("%s/trk_fragment_ana.%06i.hist",HistDir,RunNumber));

  TCanvas* c = new TCanvas("c","c",1200,800);
  c->Divide(3,2);

  c->cd(1);
  gPad->SetLogy(1);
  TH1* h1 = (TH1*) f->Get("//TrkFragmentAna/trk/nhits");
  h1->Draw();
  
  c->cd(2);
  gPad->SetLogy(1);
  TH1* h2 = (TH1*) f->Get("//TrkFragmentAna/trk/nbtot");
  h2->Draw();
  
  c->cd(3);
  gPad->SetLogy(1);
  TH1* h3 = (TH1*) f->Get("//TrkFragmentAna/trk/nfrag");
  h3->GetXaxis()->SetRangeUser(0,19.9);
  h3->Draw();
  
  c->cd(4);
  gPad->SetLogy(1);
  TH1* h4 = (TH1*) f->Get("//TrkFragmentAna/trk/fsize");
  h4->Draw();
  
  c->cd(5);
  gPad->SetLogy(1);
  TH1* h5 = (TH1*) f->Get("//TrkFragmentAna/trk/error");
  h5->Draw();
  
  c->cd(6);
  gPad->SetLogy(1);
  TH1* h6 = (TH1*) f->Get("//TrkFragmentAna/trk/valid");
  h6->Draw();
  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_nhits_per_panel(int RunNumber, int Dtc) {

  TFile* f = TFile::Open(Form("%s/trk_fragment_ana.%06i.hist",HistDir,RunNumber));

  TCanvas* c = new TCanvas(Form("c_nh_per_panel_dtc_%i",Dtc),Form("c_nh_per_panel_dtc_%i",Dtc),1200,800);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/trk/dtc_%03i/roc_%i/nhits",Dtc,i));
    h->Draw();
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int plot_nh_vs_ch_per_panel(int RunNumber, int Dtc, int MaxChannel = -1, int MaxNHits = -1) {

  TFile* f = TFile::Open(Form("%s/trk_fragment_ana.%06i.hist",HistDir,RunNumber));

  TCanvas* c = new TCanvas(Form("c_nh_vs_ch_dtc_%i",Dtc),Form("c_nh_vs_ch_dtc_%i",Dtc),1400,900);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH2* h = (TH2*) f->Get(Form("//TrkFragmentAna/trk/dtc_%03i/roc_%i/nh_vs_ch",Dtc,i));

    if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    h->SetFillColor(kBlue+2);
    h->Draw("box");
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// use default readout map
//-----------------------------------------------------------------------------
int plot_nh_vs_adc_per_panel(int RunNumber, int Dtc, int MaxChannel=-1, int MaxNHits = -1) {

  TFile* f = TFile::Open(Form("%s/trk_fragment_ana.%06i.hist",HistDir,RunNumber));

  TCanvas* c = new TCanvas(Form("c_nh_vs_adc_0_dtc_%i",Dtc),Form("c_nh_vs_adc_0_dtc_%i",Dtc),1400,900);
  c->Divide(3,2);

  for (int i=0; i<6; i++) {
    c->cd(i+1);
    gPad->SetLogy(1);

    TH1* h = (TH1*) f->Get(Form("//TrkFragmentAna/trk/dtc_%03i/roc_%i/nh_vs_adc_0",Dtc,i));

    if (MaxChannel > 0) h->GetXaxis()->SetRangeUser(0,MaxChannel-0.0001);
    if (MaxNHits   > 0) h->GetYaxis()->SetRangeUser(0,MaxNHits-0.0001);

    // h->SetFillColor(kBlue+2);
    // h->Draw("box");
    h->Draw("");
  }
  
  return 0;
}
