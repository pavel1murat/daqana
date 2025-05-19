//


TTree* get_tree(int Run) {
  TString datadir = "/data/tracker/vst/mu2etrk_daquser_001/vst00s000r000n002";
  
  TString fn = Form("%s/nts.mu2e.trk.vst00s000r000n002.%06i_000001.root",datadir.Data(),Run);
  TFile*  f  = (TFile*) gROOT->GetListOfFiles()->FindObject(fn.Data());
  if (!f || !f->IsOpen()) f = new TFile(fn.Data());

  TString dirname = Form("%s:/MakeDigiNtuple",fn.Data());
  TDirectory * dir = (TDirectory*)f->Get(dirname.Data());
  TTree* t;
  dir->GetObject("digis",t);

  return t;
}

//-----------------------------------------------------------------------------
// get_hist(106236,"plnset_00/MN224","edep")
//-----------------------------------------------------------------------------
TH1F* get_hist(int Run, const char* Module, const char* HistName) {
  TString datadir = "/data/tracker/vst/hist";
  
  TString fn = Form("%s/hst.mu2e.vst00s000r000n000.make_station_hist.%06i_000001.root",datadir.Data(),Run);
  TFile*  f  = (TFile*) gROOT->GetListOfFiles()->FindObject(fn.Data());
  if (!f || !f->IsOpen()) f = new TFile(fn.Data());

  TString dirname = Form("%s:/StationAna/%s",fn.Data(),Module);
  printf("dirname:%s\n",dirname.Data());
  TDirectory * dir = (TDirectory*)f->Get(dirname.Data());
  TH1F* hist;
  dir->GetObject(HistName,hist);

  return hist;
}


//-----------------------------------------------------------------------------
void plot_edep_sel(int Run, int MnID) {
  auto h1 = get_hist(Run,Form("pnlset_00/MN%03i",MnID),"edep");
  auto h2 = get_hist(Run,Form("pnlset_00/MN%03i",MnID),"edepg");

  h1->Draw();

  printf("h1 above 0.0006:%10.2lf\n",h1->Integral(27,100));

  h2->SetLineColor(kRed+2);
  h2->SetFillColor(kRed+2);
  h2->SetFillStyle(3005);
  h2->Draw("sames");
}

//-----------------------------------------------------------------------------
float livetime(int Run) {
  float lt = -1;
  if      (Run == 107176) lt =   2.85;
  else if (Run == 107235) lt =  82.65;
  else if (Run == 107236) lt = 103.50;
  else {
    printf("ERROR: no livetime for run %06i\n",Run);
  }
  return lt;
}

//-----------------------------------------------------------------------------
// how to normalize
//-----------------------------------------------------------------------------
void plot_occupg(int Run1, int Run2, int MnID) {
  
  auto h1 = (TH1F*) get_hist(Run1,Form("pnlset_00/MN%03i",MnID),"occupg")->Clone(Form("h_run_%06i",Run1));
  auto h2 = (TH1F*) get_hist(Run2,Form("pnlset_00/MN%03i",MnID),"occupg")->Clone(Form("h_run_%06i",Run2));

  TCanvas* c = new TCanvas(Form("c_MN%03i",MnID),Form("c_MN%03i",MnID),1000,800);
  
  float lt1 = livetime(Run1);
  h1->Scale(1./lt1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(kBlue+2);
  h1->SetLineColor(kBlue+2);
  h1->SetMarkerSize (0.8);
  h1->Draw();

  //  printf("h1 above 0.0006:%10.2lf\n",h1->Integral(27,100));

  float lt2 = livetime(Run2);
  h2->SetLineColor(kRed+1);
  // h2->SetFillColor(kRed+2);
  // h2->SetFillStyle(3005);
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(kRed+1);
  h2->SetMarkerSize (0.8);
  h2->Scale(1./lt2);
  h2->Draw("sames");

  c->Modified();
  c->Update();

  auto ps1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");
  ps1->SetTextColor(kBlue+2);

  auto ps2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
  ps2->SetTextColor(kRed+1);

  c->Modified();
  c->Update();

}


//-----------------------------------------------------------------------------
// how to normalize
//-----------------------------------------------------------------------------
void plot_occup_and_occupg(int Run, int MnID, int MakeNewCanvas = 0) {
  
  auto h1 = (TH1F*) get_hist(Run,Form("pnlset_00/MN%03i",MnID),"occup")->Clone(Form("h_occup_%06i",Run));
  auto h2 = (TH1F*) get_hist(Run,Form("pnlset_00/MN%03i",MnID),"occupg")->Clone(Form("h_occupg_%06i",Run));

  if (MakeNewCanvas) {
    TCanvas* c = new TCanvas(Form("c_%06i_MN%03i",Run,MnID),Form("c_%06i_MN%03i",Run,MnID),1200,800);
  }
  
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(kBlue+2);
  h1->SetLineColor(kBlue+2);
  h1->GetYaxis()->SetRangeUser(0.5,1.e8);
  h1->Draw();

  h2->SetLineColor(kRed+1);
  h2->SetFillStyle(3005);
  h2->SetFillColor(kRed+1);
  h2->SetMarkerStyle(20);
  h2->Draw("sames");

  gPad->SetLogy(1);
  gPad->Modified();
  
  gPad->Update();

  auto ps1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");
  ps1->SetTextColor(kBlue+2);

  auto ps2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
  ps2->SetTextColor(kRed+1);

  gPad->Modified();
  gPad->Update();

}


//-----------------------------------------------------------------------------
// how to normalize
//-----------------------------------------------------------------------------
void plot_nhits(int Run, int MnID, int MakeNewCanvas = 0) {
  
  auto h1 = (TH1F*) get_hist(Run,Form("pnlset_00/MN%03i",MnID),"nsht")->Clone(Form("h_nsht_%06i",Run));
  auto h2 = (TH1F*) get_hist(Run,Form("pnlset_00/MN%03i",MnID),"nshg")->Clone(Form("h_nshg_%06i",Run));

  if (MakeNewCanvas) {
    TCanvas* c = new TCanvas(Form("c_%06i_MN%03i",Run,MnID),Form("c_%06i_MN%03i",Run,MnID),1200,800);
  }
  
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(kBlue+2);
  h1->SetLineColor(kBlue+2);
  h1->GetXaxis()->SetRangeUser(-0.5,99.5);
  h1->GetYaxis()->SetRangeUser(0.5,1.e6);
  h1->Draw();

  h2->SetLineColor(kRed+1);
  h2->SetFillStyle(3005);
  h2->SetFillColor(kRed+1);
  h2->SetMarkerStyle(20);
  h2->Draw("sames");

  gPad->SetLogy(1);
  gPad->Modified();
  
  gPad->Update();

  auto ps1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");
  ps1->SetTextColor(kBlue+2);

  auto ps2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
  ps2->SetTextColor(kRed+1);

  gPad->Modified();
  gPad->Update();

}


//-----------------------------------------------------------------------------
int plot_hist(int Run1, int Run2, const char* Var, const char* Sel) {
  TTree *t1, *t2;

  t1 = get_tree(Run1);
  t2 = get_tree(Run2);

  std::cout << "t1:" << t1 << endl;

  t1->Draw(Var,Sel);
  t2->Draw(Var,Sel,"sames");
  
  return 0;
}
