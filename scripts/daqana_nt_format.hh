#ifndef __daqana_nt_format_hh__
#define __daqana_nt_format_hh__

  DaqEvent*   fEvent(nullptr);

 //   static constexpr Int_t kMaxsd = 1;
 //   static constexpr Int_t kMaxsh = 4090;
 //   static constexpr Int_t kMaxch = 1;
 //   static constexpr Int_t kMaxtc = 72;
 //   static constexpr Int_t kMaxtrk = 1;
 //   static constexpr Int_t kMaxtrksh = 1;
 //   static constexpr Int_t kMaxseg = 128;
 //   static constexpr Int_t kMaxsegsh = 3579;
 //   static constexpr Int_t kMaxcrvd = 1;
 //   static constexpr Int_t kMaxcrvp = 1;
 //   static constexpr Int_t kMaxcrvc = 1;

 //   // Declaration of leaf types
 // //DaqEvent        *evt;
 //   Int_t           run;
 //   Int_t           srn;
 //   Int_t           evn;
 //   Int_t           nsdtot;
 //   Int_t           sd_;
 //   UInt_t          sd_fUniqueID[kMaxsd];   //[sd_]
 //   UInt_t          sd_fBits[kMaxsd];   //[sd_]
 //   Int_t           sd__ns[kMaxsd];   //[sd_]
 //   Int_t           sd_sid[kMaxsd];   //[sd_]
 //   Int_t           sd_mnid[kMaxsd];   //[sd_]
 //   Int_t           sd_tdc0[kMaxsd];   //[sd_]
 //   Int_t           sd_tdc1[kMaxsd];   //[sd_]
 //   Int_t           sd_tot0[kMaxsd];   //[sd_]
 //   Int_t           sd_tot1[kMaxsd];   //[sd_]
 //   Int_t           sd_pmp[kMaxsd];   //[sd_]
 //   Int_t           sd_flag[kMaxsd];   //[sd_]
 //   Int_t           sd_fs[kMaxsd];   //[sd_]
 //   Float_t         sd_bl[kMaxsd];   //[sd_]
 //   Float_t         sd_ph[kMaxsd];   //[sd_]
 //   vector<unsigned short> sd_adc[kMaxsd];
 //   Int_t           nshtot;
 //   Int_t           nsh[36][6];
 //   Int_t           pmp[36];
 //   Int_t           sh_;
 //   UInt_t          sh_fUniqueID[kMaxsh];   //[sh_]
 //   UInt_t          sh_fBits[kMaxsh];   //[sh_]
 //   Int_t           sh_sid[kMaxsh];   //[sh_]
 //   Int_t           sh_zface[kMaxsh];   //[sh_]
 //   Int_t           sh_mnid[kMaxsh];   //[sh_]
 //   Float_t         sh_time[kMaxsh];   //[sh_]
 //   Float_t         sh_dt[kMaxsh];   //[sh_]
 //   Float_t         sh_tot0[kMaxsh];   //[sh_]
 //   Float_t         sh_tot1[kMaxsh];   //[sh_]
 //   Float_t         sh_edep[kMaxsh];   //[sh_]
 //   Float_t         sh_dped[kMaxsh];   //[sh_]
 //   Float_t         sh_dpmp[kMaxsh];   //[sh_]
 //   Float_t         maxEdep;
 //   Int_t           nch;
 //   Int_t           ch_;
 //   UInt_t          ch_fUniqueID[kMaxch];   //[ch_]
 //   UInt_t          ch_fBits[kMaxch];   //[ch_]
 //   Int_t           ch_sid[kMaxch];   //[ch_]
 //   Int_t           ch_nsh[kMaxch];   //[ch_]
 //   Int_t           ch_zface[kMaxch];   //[ch_]
 //   Int_t           ch_mnid[kMaxch];   //[ch_]
 //   Float_t         ch_time[kMaxch];   //[ch_]
 //   Float_t         ch_dtime[kMaxch];   //[ch_]
 //   Float_t         ch_x[kMaxch];   //[ch_]
 //   Float_t         ch_y[kMaxch];   //[ch_]
 //   Float_t         ch_z[kMaxch];   //[ch_]
 //   Float_t         ch_ux[kMaxch];   //[ch_]
 //   Float_t         ch_uy[kMaxch];   //[ch_]
 //   Float_t         ch_ures[kMaxch];   //[ch_]
 //   Float_t         ch_vres[kMaxch];   //[ch_]
 //   Float_t         ch_edep[kMaxch];   //[ch_]
 //   Int_t           ntc;
 //   Int_t           tc_;
 //   UInt_t          tc_fUniqueID[kMaxtc];   //[tc_]
 //   UInt_t          tc_fBits[kMaxtc];   //[tc_]
 //   Int_t           tc_nsh[kMaxtc];   //[tc_]
 //   Int_t           tc_nch[kMaxtc];   //[tc_]
 //   Float_t         tc_t0[kMaxtc];   //[tc_]
 //   Float_t         tc_tmin[kMaxtc];   //[tc_]
 //   Float_t         tc_tmax[kMaxtc];   //[tc_]
 //   Float_t         tc_edep_max[kMaxtc];   //[tc_]
 //   Float_t         tc_y0[kMaxtc];   //[tc_]
 //   Float_t         tc_dydz[kMaxtc];   //[tc_]
 //   Float_t         tc_chi2yz[kMaxtc];   //[tc_]
 //   Int_t           tc_nplanes[kMaxtc];   //[tc_]
 //   Int_t           tc_nfaces[kMaxtc];   //[tc_]
 //   Int_t           tc_npanels[kMaxtc];   //[tc_]
 //   Int_t           tc__nhf[kMaxtc][18][4];   //[tc_]
 //   Float_t         tc__timef[kMaxtc][18][4];   //[tc_]
 //   Int_t           tc__nhp[kMaxtc][18][2];   //[tc_]
 //   Float_t         tc__timep[kMaxtc][18][2];   //[tc_]
 //   Int_t           tc__mnid[kMaxtc][18][12];   //[tc_]
 //   Int_t           tc__nh_panel[kMaxtc][18][12];   //[tc_]
 //   Float_t         tc__time_panel[kMaxtc][18][12];   //[tc_]
 //   Float_t         tc__edep_panel[kMaxtc][18][12];   //[tc_]
 //   Int_t           tc_max_nh_panel[kMaxtc];   //[tc_]
 //   Int_t           tc_ngh[kMaxtc];   //[tc_]
 //   Int_t           ntrk;
 //   Int_t           trk_;
 //   UInt_t          trk_fUniqueID[kMaxtrk];   //[trk_]
 //   UInt_t          trk_fBits[kMaxtrk];   //[trk_]
 //   Int_t           trk_nhits[kMaxtrk];   //[trk_]
 //   Float_t         trk_t0[kMaxtrk];   //[trk_]
 //   Float_t         trk_chi2[kMaxtrk];   //[trk_]
 //   Int_t           ntrksh;
 //   Int_t           trksh_;
 //   UInt_t          trksh_fUniqueID[kMaxtrksh];   //[trksh_]
 //   UInt_t          trksh_fBits[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_sid[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_zface[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_mnid[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_time[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_dt[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_tot0[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_tot1[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_edep[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_iseg[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_itrk[kMaxtrksh];   //[trksh_]
 //   Int_t           trksh_ihit[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_rdrift[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_doca[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_dr[kMaxtrksh];   //[trksh_]
 //   Float_t         trksh_drho[kMaxtrksh];   //[trksh_]
 //   Int_t           nseg;
 //   Int_t           seg_;
 //   UInt_t          seg_fUniqueID[kMaxseg];   //[seg_]
 //   UInt_t          seg_fBits[kMaxseg];   //[seg_]
 //   Int_t           seg_sid[kMaxseg];   //[seg_]
 //   Int_t           seg_nh[kMaxseg];   //[seg_]
 //   Int_t           seg_ngh[kMaxseg];   //[seg_]
 //   Int_t           seg_nghl[kMaxseg][2];   //[seg_]
 //   Int_t           seg_nmhl[kMaxseg][2];   //[seg_]
 //   Int_t           seg_ntrans[kMaxseg];   //[seg_]
 //   Float_t         seg_t0[kMaxseg];   //[seg_]
 //   Float_t         seg_chi2d[kMaxseg];   //[seg_]
 //   Float_t         seg_y0[kMaxseg];   //[seg_]
 //   Float_t         seg_z0[kMaxseg];   //[seg_]
 //   Float_t         seg_ymean[kMaxseg];   //[seg_]
 //   Float_t         seg_dzdy[kMaxseg];   //[seg_]
 //   Float_t         seg_y0t[kMaxseg];   //[seg_]
 //   Float_t         seg_dzdyt[kMaxseg];   //[seg_]
 //   Int_t           nsegsh;
 //   Int_t           segsh_;
 //   UInt_t          segsh_fUniqueID[kMaxsegsh];   //[segsh_]
 //   UInt_t          segsh_fBits[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_sid[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_zface[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_mnid[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_time[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_dt[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_tot0[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_tot1[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_edep[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_iseg[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_itrk[kMaxsegsh];   //[segsh_]
 //   Int_t           segsh_ihit[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_rdrift[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_doca[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_dr[kMaxsegsh];   //[segsh_]
 //   Float_t         segsh_drho[kMaxsegsh];   //[segsh_]
 //   Int_t           ncrvd;
 //   Int_t           crvd_;
 //   UInt_t          crvd_fUniqueID[kMaxcrvd];   //[crvd_]
 //   UInt_t          crvd_fBits[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_sbid[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_tdc[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_nzs[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_odd_ts[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_sipm[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_roc[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_feb[kMaxcrvd];   //[crvd_]
 //   Int_t           crvd_ch[kMaxcrvd];   //[crvd_]
 //   Int_t           ncrvp;
 //   Int_t           crvp_;
 //   UInt_t          crvp_fUniqueID[kMaxcrvp];   //[crvp_]
 //   UInt_t          crvp_fBits[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_npes[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_pes_ph[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_time[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_ph[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_ped[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_beta[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_chi2[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_le_time[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_flags[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_npes_nofit[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_time_nofit[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_tstart[kMaxcrvp];   //[crvp_]
 //   Float_t         crvp_tend[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_sbid[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_sipm[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_roc[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_feb[kMaxcrvp];   //[crvp_]
 //   Int_t           crvp_ch[kMaxcrvp];   //[crvp_]
 //   Int_t           ncrvc;
 //   Int_t           crvc_;
 //   UInt_t          crvc_fUniqueID[kMaxcrvc];   //[crvc_]
 //   UInt_t          crvc_fBits[kMaxcrvc];   //[crvc_]
 //   Int_t           crvc_stype[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_tstart[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_tend[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_pes[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_time[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_x[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_y[kMaxcrvc];   //[crvc_]
 //   Float_t         crvc_z[kMaxcrvc];   //[crvc_]
 //   Int_t           crvc_nlayers[kMaxcrvc];   //[crvc_]
 //   Int_t           crvc_nsides[kMaxcrvc];   //[crvc_]

 //   // List of branches
 //   TBranch        *b_evt_run;   //!
 //   TBranch        *b_evt_srn;   //!
 //   TBranch        *b_evt_evn;   //!
 //   TBranch        *b_evt_nsdtot;   //!
 //   TBranch        *b_evt_sd_;   //!
 //   TBranch        *b_sd_fUniqueID;   //!
 //   TBranch        *b_sd_fBits;   //!
 //   TBranch        *b_sd__ns;   //!
 //   TBranch        *b_sd_sid;   //!
 //   TBranch        *b_sd_mnid;   //!
 //   TBranch        *b_sd_tdc0;   //!
 //   TBranch        *b_sd_tdc1;   //!
 //   TBranch        *b_sd_tot0;   //!
 //   TBranch        *b_sd_tot1;   //!
 //   TBranch        *b_sd_pmp;   //!
 //   TBranch        *b_sd_flag;   //!
 //   TBranch        *b_sd_fs;   //!
 //   TBranch        *b_sd_bl;   //!
 //   TBranch        *b_sd_ph;   //!
 //   TBranch        *b_sd_adc;   //!
 //   TBranch        *b_evt_nshtot;   //!
 //   TBranch        *b_evt_nsh;   //!
 //   TBranch        *b_evt_pmp;   //!
 //   TBranch        *b_evt_sh_;   //!
 //   TBranch        *b_sh_fUniqueID;   //!
 //   TBranch        *b_sh_fBits;   //!
 //   TBranch        *b_sh_sid;   //!
 //   TBranch        *b_sh_zface;   //!
 //   TBranch        *b_sh_mnid;   //!
 //   TBranch        *b_sh_time;   //!
 //   TBranch        *b_sh_dt;   //!
 //   TBranch        *b_sh_tot0;   //!
 //   TBranch        *b_sh_tot1;   //!
 //   TBranch        *b_sh_edep;   //!
 //   TBranch        *b_sh_dped;   //!
 //   TBranch        *b_sh_dpmp;   //!
 //   TBranch        *b_evt_maxEdep;   //!
 //   TBranch        *b_evt_nch;   //!
 //   TBranch        *b_evt_ch_;   //!
 //   TBranch        *b_ch_fUniqueID;   //!
 //   TBranch        *b_ch_fBits;   //!
 //   TBranch        *b_ch_sid;   //!
 //   TBranch        *b_ch_nsh;   //!
 //   TBranch        *b_ch_zface;   //!
 //   TBranch        *b_ch_mnid;   //!
 //   TBranch        *b_ch_time;   //!
 //   TBranch        *b_ch_dtime;   //!
 //   TBranch        *b_ch_x;   //!
 //   TBranch        *b_ch_y;   //!
 //   TBranch        *b_ch_z;   //!
 //   TBranch        *b_ch_ux;   //!
 //   TBranch        *b_ch_uy;   //!
 //   TBranch        *b_ch_ures;   //!
 //   TBranch        *b_ch_vres;   //!
 //   TBranch        *b_ch_edep;   //!
 //   TBranch        *b_evt_ntc;   //!
 //   TBranch        *b_evt_tc_;   //!
 //   TBranch        *b_tc_fUniqueID;   //!
 //   TBranch        *b_tc_fBits;   //!
 //   TBranch        *b_tc_nsh;   //!
 //   TBranch        *b_tc_nch;   //!
 //   TBranch        *b_tc_t0;   //!
 //   TBranch        *b_tc_tmin;   //!
 //   TBranch        *b_tc_tmax;   //!
 //   TBranch        *b_tc_edep_max;   //!
 //   TBranch        *b_tc_y0;   //!
 //   TBranch        *b_tc_dydz;   //!
 //   TBranch        *b_tc_chi2yz;   //!
 //   TBranch        *b_tc_nplanes;   //!
 //   TBranch        *b_tc_nfaces;   //!
 //   TBranch        *b_tc_npanels;   //!
 //   TBranch        *b_tc__nhf;   //!
 //   TBranch        *b_tc__timef;   //!
 //   TBranch        *b_tc__nhp;   //!
 //   TBranch        *b_tc__timep;   //!
 //   TBranch        *b_tc__mnid;   //!
 //   TBranch        *b_tc__nh_panel;   //!
 //   TBranch        *b_tc__time_panel;   //!
 //   TBranch        *b_tc__edep_panel;   //!
 //   TBranch        *b_tc_max_nh_panel;   //!
 //   TBranch        *b_tc_ngh;   //!
 //   TBranch        *b_evt_ntrk;   //!
 //   TBranch        *b_evt_trk_;   //!
 //   TBranch        *b_trk_fUniqueID;   //!
 //   TBranch        *b_trk_fBits;   //!
 //   TBranch        *b_trk_nhits;   //!
 //   TBranch        *b_trk_t0;   //!
 //   TBranch        *b_trk_chi2;   //!
 //   TBranch        *b_evt_ntrksh;   //!
 //   TBranch        *b_evt_trksh_;   //!
 //   TBranch        *b_trksh_fUniqueID;   //!
 //   TBranch        *b_trksh_fBits;   //!
 //   TBranch        *b_trksh_sid;   //!
 //   TBranch        *b_trksh_zface;   //!
 //   TBranch        *b_trksh_mnid;   //!
 //   TBranch        *b_trksh_time;   //!
 //   TBranch        *b_trksh_dt;   //!
 //   TBranch        *b_trksh_tot0;   //!
 //   TBranch        *b_trksh_tot1;   //!
 //   TBranch        *b_trksh_edep;   //!
 //   TBranch        *b_trksh_iseg;   //!
 //   TBranch        *b_trksh_itrk;   //!
 //   TBranch        *b_trksh_ihit;   //!
 //   TBranch        *b_trksh_rdrift;   //!
 //   TBranch        *b_trksh_doca;   //!
 //   TBranch        *b_trksh_dr;   //!
 //   TBranch        *b_trksh_drho;   //!
 //   TBranch        *b_evt_nseg;   //!
 //   TBranch        *b_evt_seg_;   //!
 //   TBranch        *b_seg_fUniqueID;   //!
 //   TBranch        *b_seg_fBits;   //!
 //   TBranch        *b_seg_sid;   //!
 //   TBranch        *b_seg_nh;   //!
 //   TBranch        *b_seg_ngh;   //!
 //   TBranch        *b_seg_nghl;   //!
 //   TBranch        *b_seg_nmhl;   //!
 //   TBranch        *b_seg_ntrans;   //!
 //   TBranch        *b_seg_t0;   //!
 //   TBranch        *b_seg_chi2d;   //!
 //   TBranch        *b_seg_y0;   //!
 //   TBranch        *b_seg_z0;   //!
 //   TBranch        *b_seg_ymean;   //!
 //   TBranch        *b_seg_dzdy;   //!
 //   TBranch        *b_seg_y0t;   //!
 //   TBranch        *b_seg_dzdyt;   //!
 //   TBranch        *b_evt_nsegsh;   //!
 //   TBranch        *b_evt_segsh_;   //!
 //   TBranch        *b_segsh_fUniqueID;   //!
 //   TBranch        *b_segsh_fBits;   //!
 //   TBranch        *b_segsh_sid;   //!
 //   TBranch        *b_segsh_zface;   //!
 //   TBranch        *b_segsh_mnid;   //!
 //   TBranch        *b_segsh_time;   //!
 //   TBranch        *b_segsh_dt;   //!
 //   TBranch        *b_segsh_tot0;   //!
 //   TBranch        *b_segsh_tot1;   //!
 //   TBranch        *b_segsh_edep;   //!
 //   TBranch        *b_segsh_iseg;   //!
 //   TBranch        *b_segsh_itrk;   //!
 //   TBranch        *b_segsh_ihit;   //!
 //   TBranch        *b_segsh_rdrift;   //!
 //   TBranch        *b_segsh_doca;   //!
 //   TBranch        *b_segsh_dr;   //!
 //   TBranch        *b_segsh_drho;   //!
 //   TBranch        *b_evt_ncrvd;   //!
 //   TBranch        *b_evt_crvd_;   //!
 //   TBranch        *b_crvd_fUniqueID;   //!
 //   TBranch        *b_crvd_fBits;   //!
 //   TBranch        *b_crvd_sbid;   //!
 //   TBranch        *b_crvd_tdc;   //!
 //   TBranch        *b_crvd_nzs;   //!
 //   TBranch        *b_crvd_odd_ts;   //!
 //   TBranch        *b_crvd_sipm;   //!
 //   TBranch        *b_crvd_roc;   //!
 //   TBranch        *b_crvd_feb;   //!
 //   TBranch        *b_crvd_ch;   //!
 //   TBranch        *b_evt_ncrvp;   //!
 //   TBranch        *b_evt_crvp_;   //!
 //   TBranch        *b_crvp_fUniqueID;   //!
 //   TBranch        *b_crvp_fBits;   //!
 //   TBranch        *b_crvp_npes;   //!
 //   TBranch        *b_crvp_pes_ph;   //!
 //   TBranch        *b_crvp_time;   //!
 //   TBranch        *b_crvp_ph;   //!
 //   TBranch        *b_crvp_ped;   //!
 //   TBranch        *b_crvp_beta;   //!
 //   TBranch        *b_crvp_chi2;   //!
 //   TBranch        *b_crvp_le_time;   //!
 //   TBranch        *b_crvp_flags;   //!
 //   TBranch        *b_crvp_npes_nofit;   //!
 //   TBranch        *b_crvp_time_nofit;   //!
 //   TBranch        *b_crvp_tstart;   //!
 //   TBranch        *b_crvp_tend;   //!
 //   TBranch        *b_crvp_sbid;   //!
 //   TBranch        *b_crvp_sipm;   //!
 //   TBranch        *b_crvp_roc;   //!
 //   TBranch        *b_crvp_feb;   //!
 //   TBranch        *b_crvp_ch;   //!
 //   TBranch        *b_evt_ncrvc;   //!
 //   TBranch        *b_evt_crvc_;   //!
 //   TBranch        *b_crvc_fUniqueID;   //!
 //   TBranch        *b_crvc_fBits;   //!
 //   TBranch        *b_crvc_stype;   //!
 //   TBranch        *b_crvc_tstart;   //!
 //   TBranch        *b_crvc_tend;   //!
 //   TBranch        *b_crvc_pes;   //!
 //   TBranch        *b_crvc_time;   //!
 //   TBranch        *b_crvc_x;   //!
 //   TBranch        *b_crvc_y;   //!
 //   TBranch        *b_crvc_z;   //!
 //   TBranch        *b_crvc_nlayers;   //!
 //   TBranch        *b_crvc_nsides;   //!

#endif
