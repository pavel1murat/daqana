///////////////////////////////////////////////////////////////////////////////
// pin = pulse injection
// do gSystem->Load("v001/.spack-env/view/lib/libdaqana_ana.so") first
// assuming ntuples for 8 pulse injection runs are produced and stored in
// /data/mu2e/mu2etrk/datasets/vst00s000r000n001/
// , make histogram file with the deltaT = T(cal)-T(hv) histograms 
///////////////////////////////////////////////////////////////////////////////
int make_pin_dt_hist() {
  int run1(122655), run2(122662);

  plot_n001_pin* x;
  
  for (int run=run1; run<run2+1; run++) {
    x = new plot_n001_pin("pin",run) ;
    x->Loop();
  }
  
  x->SaveHistograms(Form("pin_%06d_%06d.root",run1,run2));
  return 0;
}
