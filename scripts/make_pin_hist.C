///////////////////////////////////////////////////////////////////////////////
// pin = pulse injection
// do gSystem->Load("v001/.spack-env/view/lib/libdaqana_ana.so") first
// assuming ntuples for 8 pulse injection runs are produced and stored in
// /data/mu2e/mu2etrk/datasets/vst00s000r000n001/
// , make histogram file with the dt01 = T(cal)-T(hv) histograms 
///////////////////////////////////////////////////////////////////////////////
int make_pin_hist() {
  // int run1(122655), run2(122662);
  // run 122975 failed
  int rn[8] = {122968, 122969, 122970, 122971, 122972, 122973, 122974, 122976};

  plot_n001_pin* x;
  
  for (int i=0; i<8; i++) {
    int run = rn[i];
    x = new plot_n001_pin(run) ;
    x->Loop();
  }
  
  x->SaveHistograms(Form("hst.mu2e.trk.plot_n001_pin.%06d_%06d.root",rn[0],rn[7]));
  return 0;
}
