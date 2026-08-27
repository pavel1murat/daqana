//-----------------------------------------------------------------------------
// T1 and T2 are the two clones of the same tree
//-----------------------------------------------------------------------------
#include <iostream>
#include "TH1D.h"
#include "daqana/obj/DaqEvent.hh"
#include "ana/plot_crv_tc_dt.hh"
#include "ana/check_split.hh"

//-----------------------------------------------------------------------------
int check_split(int Shift = 0) {

  DaqEvent*       e1(nullptr);               // #include "daqana_nt_format.hh"
  DaqEvent*       e2(nullptr);               // #include "daqana_nt_format.hh"

  TFile* f1 = new TFile("/data/mu2e/mu2etrk/datasets/vst00s000r000n104/nts.mu2e.trk.vst00s000r000n104.123680_000001.root");
  TFile* f2 = new TFile("/data/mu2e/mu2etrk/datasets/vst00s000r000n104/nts.mu2e.trk.vst00s000r000n104.123680_000001.root");

  TTree* t1 = (TTree*) f1->Get("/MakeDigiNtuple/digis");
  TTree* t2 = (TTree*) f2->Get("/MakeDigiNtuple/digis");
  
  int rc1 = t1->SetBranchAddress("evt",&e1);
  int rc2 = t2->SetBranchAddress("evt",&e2);
  // std::cout << rc1 << " " << rc2 << std::endl;

  //  t1->Print();
  
  TH1D* h1 = new TH1D("h1",Form("crvc.time-tc.t0, ROC=1, Shift:%d",Shift),400,-1000,1000);
  TH1D* h2 = new TH1D("h2",Form("crvc.time-tc.t0, ROC=2, Shift:%d",Shift),400,-1000,1000);
  
  Long64_t n = t1->GetEntries();
  std::cout << "0015 emoe n:" << n << std::endl;
  
  for (Long64_t i=2; i<n-2; ++i) {
    // std::cout << "001 emoe" << std::endl;
    
    int nb1 = t1->GetEntry(i  );
    // std::cout << "0011 evn1:" << e1->evn << " nb1: " << nb1 << std::endl;
    
    // std::cout << "002 emoe" << std::endl;
    // time clusters come from the first tree
    int ntc = e1->tc->GetEntriesFast();

    // std::cout << "003 emoe ntc:" << ntc << std::endl;
    
    t2->GetEntry(i+Shift);

    for (int i1=0; i1<ntc; i1++) {
      DaqTimeCluster* tc = e1->Tc(i1);
      
      for (int i2=0; i2<e2->ncrvp; i2++) {
        DaqCrvRecoPulse* crvp = e2->Crvp(i2);

        // for (int i2=0; i2<e1->ncrvp; i2++) {
      //   DaqCrvRecoPulse* crvp = e1->Crvp(i2);

        if (crvp->roc == 1) h1->Fill(crvp->time-tc->t0);
        if (crvp->roc == 2) h2->Fill(crvp->time-tc->t0);
      }
    }
  }

  auto c = new TCanvas("c","c",1200,700);
  c->Divide(2,1);
  c->cd(1);
  h1->Draw();
  c->cd(2);
  h2->Draw();

  return 0;
}
