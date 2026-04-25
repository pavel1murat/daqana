#ifndef __daqana_nt_init_C__
#define __daqana_nt_init_C__
///////////////////////////////////////////////////////////////////////////////
// body of the init function - shared

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set array pointer
   for (int i=0; i<kMaxsd; ++i) sd_adc[i] = 0;

                                        // Set branch addresses and branch pointers
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_evt_run);
   fChain->SetBranchAddress("srn", &srn, &b_evt_srn);
   fChain->SetBranchAddress("evn", &evn, &b_evt_evn);
   fChain->SetBranchAddress("nsdtot", &nsdtot, &b_evt_nsdtot);
   fChain->SetBranchAddress("sd", &sd_, &b_evt_sd_);
   fChain->SetBranchAddress("sd.fUniqueID", sd_fUniqueID, &b_sd_fUniqueID);
   fChain->SetBranchAddress("sd.fBits", sd_fBits, &b_sd_fBits);
   fChain->SetBranchAddress("sd.sid", sd_sid, &b_sd_sid);
   fChain->SetBranchAddress("sd.mnid", sd_mnid, &b_sd_mnid);
   fChain->SetBranchAddress("sd.tdc0", sd_tdc0, &b_sd_tdc0);
   fChain->SetBranchAddress("sd.tdc1", sd_tdc1, &b_sd_tdc1);
   fChain->SetBranchAddress("sd.tot0", sd_tot0, &b_sd_tot0);
   fChain->SetBranchAddress("sd.tot1", sd_tot1, &b_sd_tot1);
   fChain->SetBranchAddress("sd.pmp", sd_pmp, &b_sd_pmp);
   fChain->SetBranchAddress("sd.flag", sd_flag, &b_sd_flag);
   fChain->SetBranchAddress("sd.fs", sd_fs, &b_sd_fs);
   fChain->SetBranchAddress("sd.bl", sd_bl, &b_sd_bl);
   fChain->SetBranchAddress("sd.ph", sd_ph, &b_sd_ph);
   fChain->SetBranchAddress("sd.ns", sd_ns, &b_sd_ns);
   fChain->SetBranchAddress("sd.adc", sd_adc, &b_sd_adc);
   fChain->SetBranchAddress("nshtot", &nshtot, &b_evt_nshtot);
   fChain->SetBranchAddress("nsh[36][6]", nsh, &b_evt_nsh);
   fChain->SetBranchAddress("pmp[36]", pmp, &b_evt_pmp);
   fChain->SetBranchAddress("sh", &sh_, &b_evt_sh_);
   fChain->SetBranchAddress("sh.fUniqueID", sh_fUniqueID, &b_sh_fUniqueID);
   fChain->SetBranchAddress("sh.fBits", sh_fBits, &b_sh_fBits);
   fChain->SetBranchAddress("sh.sid", sh_sid, &b_sh_sid);
   fChain->SetBranchAddress("sh.zface", sh_zface, &b_sh_zface);
   fChain->SetBranchAddress("sh.mnid", sh_mnid, &b_sh_mnid);
   fChain->SetBranchAddress("sh.time", sh_time, &b_sh_time);
   fChain->SetBranchAddress("sh.dt", sh_dt, &b_sh_dt);
   fChain->SetBranchAddress("sh.tot0", sh_tot0, &b_sh_tot0);
   fChain->SetBranchAddress("sh.tot1", sh_tot1, &b_sh_tot1);
   fChain->SetBranchAddress("sh.edep", sh_edep, &b_sh_edep);
   fChain->SetBranchAddress("maxEdep", &maxEdep, &b_evt_maxEdep);
   fChain->SetBranchAddress("nch", &nch, &b_evt_nch);
   fChain->SetBranchAddress("ch", &ch_, &b_evt_ch_);
   fChain->SetBranchAddress("ch.fUniqueID", &ch_fUniqueID, &b_ch_fUniqueID);
   fChain->SetBranchAddress("ch.fBits", &ch_fBits, &b_ch_fBits);
   fChain->SetBranchAddress("ch.sid", &ch_sid, &b_ch_sid);
   fChain->SetBranchAddress("ch.nsh", &ch_nsh, &b_ch_nsh);
   fChain->SetBranchAddress("ch.zface", &ch_zface, &b_ch_zface);
   fChain->SetBranchAddress("ch.mnid", &ch_mnid, &b_ch_mnid);
   fChain->SetBranchAddress("ch.time", &ch_time, &b_ch_time);
   fChain->SetBranchAddress("ch.dtime", &ch_dtime, &b_ch_dtime);
   fChain->SetBranchAddress("ch.x", &ch_x, &b_ch_x);
   fChain->SetBranchAddress("ch.y", &ch_y, &b_ch_y);
   fChain->SetBranchAddress("ch.z", &ch_z, &b_ch_z);
   fChain->SetBranchAddress("ch.ux", &ch_ux, &b_ch_ux);
   fChain->SetBranchAddress("ch.uy", &ch_uy, &b_ch_uy);
   fChain->SetBranchAddress("ch.ures", &ch_ures, &b_ch_ures);
   fChain->SetBranchAddress("ch.vres", &ch_vres, &b_ch_vres);
   fChain->SetBranchAddress("ch.edep", &ch_edep, &b_ch_edep);
   fChain->SetBranchAddress("ntc", &ntc, &b_evt_ntc);
   fChain->SetBranchAddress("tc", &tc_, &b_evt_tc_);
   fChain->SetBranchAddress("tc.fUniqueID", &tc_fUniqueID, &b_tc_fUniqueID);
   fChain->SetBranchAddress("tc.fBits", &tc_fBits, &b_tc_fBits);
   fChain->SetBranchAddress("tc.nsh", &tc_nsh, &b_tc_nsh);
   fChain->SetBranchAddress("tc.nch", &tc_nch, &b_tc_nch);
   fChain->SetBranchAddress("tc.t0", &tc_t0, &b_tc_t0);
   fChain->SetBranchAddress("tc.tmin", &tc_tmin, &b_tc_tmin);
   fChain->SetBranchAddress("tc.tmax", &tc_tmax, &b_tc_tmax);
   fChain->SetBranchAddress("tc.edep_max", &tc_edep_max, &b_tc_edep_max);
   fChain->SetBranchAddress("tc.y0", &tc_y0, &b_tc_y0);
   fChain->SetBranchAddress("tc.dydz", &tc_dydz, &b_tc_dydz);
   fChain->SetBranchAddress("tc.chi2yz", &tc_chi2yz, &b_tc_chi2yz);
   fChain->SetBranchAddress("tc.nplanes", &tc_nplanes, &b_tc_nplanes);
   fChain->SetBranchAddress("tc.nfaces", &tc_nfaces, &b_tc_nfaces);
   fChain->SetBranchAddress("tc.npanels", &tc_npanels, &b_tc_npanels);
   fChain->SetBranchAddress("tc._nhf[18][4]", &tc__nhf, &b_tc__nhf);
   fChain->SetBranchAddress("tc._timef[18][4]", &tc__timef, &b_tc__timef);
   fChain->SetBranchAddress("tc._nhp[18][2]", &tc__nhp, &b_tc__nhp);
   fChain->SetBranchAddress("tc._timep[18][2]", &tc__timep, &b_tc__timep);
   fChain->SetBranchAddress("tc._mnid[18][12]", &tc__mnid, &b_tc__mnid);
   fChain->SetBranchAddress("tc._nh_panel[18][12]", &tc__nh_panel, &b_tc__nh_panel);
   fChain->SetBranchAddress("tc._time_panel[18][12]", &tc__time_panel, &b_tc__time_panel);
   fChain->SetBranchAddress("tc._edep_panel[18][12]", &tc__edep_panel, &b_tc__edep_panel);
   fChain->SetBranchAddress("tc.max_nh_panel", &tc_max_nh_panel, &b_tc_max_nh_panel);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_evt_ntrk);
   fChain->SetBranchAddress("trk", &trk_, &b_evt_trk_);
   fChain->SetBranchAddress("trk.fUniqueID", &trk_fUniqueID, &b_trk_fUniqueID);
   fChain->SetBranchAddress("trk.fBits", &trk_fBits, &b_trk_fBits);
   fChain->SetBranchAddress("trk.nhits", &trk_nhits, &b_trk_nhits);
   fChain->SetBranchAddress("trk.t0", &trk_t0, &b_trk_t0);
   fChain->SetBranchAddress("trk.chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("ntrksh", &ntrksh, &b_evt_ntrksh);
   fChain->SetBranchAddress("trksh", &trksh_, &b_evt_trksh_);
   fChain->SetBranchAddress("trksh.fUniqueID", &trksh_fUniqueID, &b_trksh_fUniqueID);
   fChain->SetBranchAddress("trksh.fBits", &trksh_fBits, &b_trksh_fBits);
   fChain->SetBranchAddress("trksh.sid", &trksh_sid, &b_trksh_sid);
   fChain->SetBranchAddress("trksh.zface", &trksh_zface, &b_trksh_zface);
   fChain->SetBranchAddress("trksh.mnid", &trksh_mnid, &b_trksh_mnid);
   fChain->SetBranchAddress("trksh.time", &trksh_time, &b_trksh_time);
   fChain->SetBranchAddress("trksh.dt", &trksh_dt, &b_trksh_dt);
   fChain->SetBranchAddress("trksh.tot0", &trksh_tot0, &b_trksh_tot0);
   fChain->SetBranchAddress("trksh.tot1", &trksh_tot1, &b_trksh_tot1);
   fChain->SetBranchAddress("trksh.edep", &trksh_edep, &b_trksh_edep);
   fChain->SetBranchAddress("trksh.iseg", &trksh_iseg, &b_trksh_iseg);
   fChain->SetBranchAddress("trksh.itrk", &trksh_itrk, &b_trksh_itrk);
   fChain->SetBranchAddress("trksh.ihit", &trksh_ihit, &b_trksh_ihit);
   fChain->SetBranchAddress("trksh.rdrift", &trksh_rdrift, &b_trksh_rdrift);
   fChain->SetBranchAddress("trksh.doca", &trksh_doca, &b_trksh_doca);
   fChain->SetBranchAddress("trksh.dr", &trksh_dr, &b_trksh_dr);
   fChain->SetBranchAddress("trksh.drho", &trksh_drho, &b_trksh_drho);
   fChain->SetBranchAddress("nseg", &nseg, &b_evt_nseg);
   fChain->SetBranchAddress("seg", &seg_, &b_evt_seg_);
   fChain->SetBranchAddress("seg.fUniqueID", &seg_fUniqueID, &b_seg_fUniqueID);
   fChain->SetBranchAddress("seg.fBits", &seg_fBits, &b_seg_fBits);
   fChain->SetBranchAddress("seg.sid", &seg_sid, &b_seg_sid);
   fChain->SetBranchAddress("seg.nh", &seg_nh, &b_seg_nh);
   fChain->SetBranchAddress("seg.ngh", &seg_ngh, &b_seg_ngh);
   fChain->SetBranchAddress("seg.nghl[2]", &seg_nghl, &b_seg_nghl);
   fChain->SetBranchAddress("seg.nmhl[2]", &seg_nmhl, &b_seg_nmhl);
   fChain->SetBranchAddress("seg.ntrans", &seg_ntrans, &b_seg_ntrans);
   fChain->SetBranchAddress("seg.t0", &seg_t0, &b_seg_t0);
   fChain->SetBranchAddress("seg.chi2d", &seg_chi2d, &b_seg_chi2d);
   fChain->SetBranchAddress("seg.y0", &seg_y0, &b_seg_y0);
   fChain->SetBranchAddress("seg.z0", &seg_z0, &b_seg_z0);
   fChain->SetBranchAddress("seg.ymean", &seg_ymean, &b_seg_ymean);
   fChain->SetBranchAddress("seg.dzdy", &seg_dzdy, &b_seg_dzdy);
   fChain->SetBranchAddress("seg.y0t", &seg_y0t, &b_seg_y0t);
   fChain->SetBranchAddress("seg.dzdyt", &seg_dzdyt, &b_seg_dzdyt);
   fChain->SetBranchAddress("nsegsh", &nsegsh, &b_evt_nsegsh);
   fChain->SetBranchAddress("segsh", &segsh_, &b_evt_segsh_);
   fChain->SetBranchAddress("segsh.fUniqueID", &segsh_fUniqueID, &b_segsh_fUniqueID);
   fChain->SetBranchAddress("segsh.fBits", &segsh_fBits, &b_segsh_fBits);
   fChain->SetBranchAddress("segsh.sid", &segsh_sid, &b_segsh_sid);
   fChain->SetBranchAddress("segsh.zface", &segsh_zface, &b_segsh_zface);
   fChain->SetBranchAddress("segsh.mnid", &segsh_mnid, &b_segsh_mnid);
   fChain->SetBranchAddress("segsh.time", &segsh_time, &b_segsh_time);
   fChain->SetBranchAddress("segsh.dt", &segsh_dt, &b_segsh_dt);
   fChain->SetBranchAddress("segsh.tot0", &segsh_tot0, &b_segsh_tot0);
   fChain->SetBranchAddress("segsh.tot1", &segsh_tot1, &b_segsh_tot1);
   fChain->SetBranchAddress("segsh.edep", &segsh_edep, &b_segsh_edep);
   fChain->SetBranchAddress("segsh.iseg", &segsh_iseg, &b_segsh_iseg);
   fChain->SetBranchAddress("segsh.itrk", &segsh_itrk, &b_segsh_itrk);
   fChain->SetBranchAddress("segsh.ihit", &segsh_ihit, &b_segsh_ihit);
   fChain->SetBranchAddress("segsh.rdrift", &segsh_rdrift, &b_segsh_rdrift);
   fChain->SetBranchAddress("segsh.doca", &segsh_doca, &b_segsh_doca);
   fChain->SetBranchAddress("segsh.dr", &segsh_dr, &b_segsh_dr);
   fChain->SetBranchAddress("segsh.drho", &segsh_drho, &b_segsh_drho);
#endif
