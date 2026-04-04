//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 31 14:42:22 2026 by ROOT version 6.32.06
// from TTree digis/digis
// found on file: nts.mu2e.trk.vst00s000r000n001.120718_000001.root
//////////////////////////////////////////////////////////

#ifndef plot_n001_hist_h
#define plot_n001_hist_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.
#include "daqana/obj/DaqEvent.hh"
#include "TObject.h"
#include "daqana/obj/DaqStrawDigi.hh"
#include "daqana/obj/DaqStrawHit.hh"
#include "daqana/obj/DaqComboHit.hh"
#include "daqana/obj/DaqTimeCluster.hh"
#include "daqana/obj/DaqTrack.hh"
#include "daqana/obj/DaqTrkStrawHit.hh"
#include "daqana/obj/DaqSegment.hh"

class plot_n001_hist {
public :
  struct Hist_t {
    TH2F* fDtVsEvt;
    TH1F* fDt10k;
    TH1F* fNhCh1;
    TH1F* fNhCh2;
    TH1F* fGenDt1;
    TH1F* fGenDt2;
  };

  Hist_t         fHist;
  int            fRunNumber;
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxsd    = 1000;
   static constexpr Int_t kMaxsh    = 1000;
   static constexpr Int_t kMaxch    = 1;
   static constexpr Int_t kMaxtc    = 1;
   static constexpr Int_t kMaxtrk   = 1;
   static constexpr Int_t kMaxtrksh = 1;
   static constexpr Int_t kMaxseg   = 1;
   static constexpr Int_t kMaxsegsh = 1;

   // Declaration of leaf types
 //DaqEvent        *evt;
   Int_t           run;
   Int_t           srn;
   Int_t           evn;
   Int_t           nsdtot;
   Int_t           sd_;
   UInt_t          sd_fUniqueID[kMaxsd];   //[sd_]
   UInt_t          sd_fBits[kMaxsd];   //[sd_]
   Int_t           sd_sid[kMaxsd];   //[sd_]
   Int_t           sd_mnid[kMaxsd];   //[sd_]
   Int_t           sd_tdc0[kMaxsd];   //[sd_]
   Int_t           sd_tdc1[kMaxsd];   //[sd_]
   Int_t           sd_tot0[kMaxsd];   //[sd_]
   Int_t           sd_tot1[kMaxsd];   //[sd_]
   Int_t           sd_pmp[kMaxsd];   //[sd_]
   Int_t           sd_flag[kMaxsd];   //[sd_]
   Int_t           sd_fs[kMaxsd];   //[sd_]
   Float_t         sd_bl[kMaxsd];   //[sd_]
   Float_t         sd_ph[kMaxsd];   //[sd_]
   Int_t           sd_ns[kMaxsd];   //[sd_]
   Short_t        *sd_adc[kMaxsd];   //[sd_ns]
   Int_t           nshtot;
   Int_t           nsh[36][6];
   Int_t           pmp[36];
   Int_t           sh_;
   UInt_t          sh_fUniqueID[kMaxsh];   //[sh_]
   UInt_t          sh_fBits[kMaxsh];   //[sh_]
   Int_t           sh_sid[kMaxsh];   //[sh_]
   Int_t           sh_zface[kMaxsh];   //[sh_]
   Int_t           sh_mnid[kMaxsh];   //[sh_]
   Float_t         sh_time[kMaxsh];   //[sh_]
   Float_t         sh_dt[kMaxsh];   //[sh_]
   Float_t         sh_tot0[kMaxsh];   //[sh_]
   Float_t         sh_tot1[kMaxsh];   //[sh_]
   Float_t         sh_edep[kMaxsh];   //[sh_]
   Float_t         maxEdep;
   Int_t           nch;
   Int_t           ch_;
   UInt_t          ch_fUniqueID[kMaxch];   //[ch_]
   UInt_t          ch_fBits[kMaxch];   //[ch_]
   Int_t           ch_sid[kMaxch];   //[ch_]
   Int_t           ch_nsh[kMaxch];   //[ch_]
   Int_t           ch_zface[kMaxch];   //[ch_]
   Int_t           ch_mnid[kMaxch];   //[ch_]
   Float_t         ch_time[kMaxch];   //[ch_]
   Float_t         ch_dtime[kMaxch];   //[ch_]
   Float_t         ch_x[kMaxch];   //[ch_]
   Float_t         ch_y[kMaxch];   //[ch_]
   Float_t         ch_z[kMaxch];   //[ch_]
   Float_t         ch_ux[kMaxch];   //[ch_]
   Float_t         ch_uy[kMaxch];   //[ch_]
   Float_t         ch_ures[kMaxch];   //[ch_]
   Float_t         ch_vres[kMaxch];   //[ch_]
   Float_t         ch_edep[kMaxch];   //[ch_]
   Int_t           ntc;
   Int_t           tc_;
   UInt_t          tc_fUniqueID[kMaxtc];   //[tc_]
   UInt_t          tc_fBits[kMaxtc];   //[tc_]
   Int_t           tc_nsh[kMaxtc];   //[tc_]
   Int_t           tc_nch[kMaxtc];   //[tc_]
   Float_t         tc_t0[kMaxtc];   //[tc_]
   Float_t         tc_tmin[kMaxtc];   //[tc_]
   Float_t         tc_tmax[kMaxtc];   //[tc_]
   Float_t         tc_edep_max[kMaxtc];   //[tc_]
   Float_t         tc_y0[kMaxtc];   //[tc_]
   Float_t         tc_dydz[kMaxtc];   //[tc_]
   Float_t         tc_chi2yz[kMaxtc];   //[tc_]
   Int_t           tc_nplanes[kMaxtc];   //[tc_]
   Int_t           tc_nfaces[kMaxtc];   //[tc_]
   Int_t           tc_npanels[kMaxtc];   //[tc_]
   Int_t           tc__nhf[kMaxtc][18][4];   //[tc_]
   Float_t         tc__timef[kMaxtc][18][4];   //[tc_]
   Int_t           tc__nhp[kMaxtc][18][2];   //[tc_]
   Float_t         tc__timep[kMaxtc][18][2];   //[tc_]
   Int_t           tc__mnid[kMaxtc][18][12];   //[tc_]
   Int_t           tc__nh_panel[kMaxtc][18][12];   //[tc_]
   Float_t         tc__time_panel[kMaxtc][18][12];   //[tc_]
   Float_t         tc__edep_panel[kMaxtc][18][12];   //[tc_]
   Int_t           tc_max_nh_panel[kMaxtc];   //[tc_]
   Int_t           ntrk;
   Int_t           trk_;
   UInt_t          trk_fUniqueID[kMaxtrk];   //[trk_]
   UInt_t          trk_fBits[kMaxtrk];   //[trk_]
   Int_t           trk_nhits[kMaxtrk];   //[trk_]
   Float_t         trk_t0[kMaxtrk];   //[trk_]
   Float_t         trk_chi2[kMaxtrk];   //[trk_]
   Int_t           ntrksh;
   Int_t           trksh_;
   UInt_t          trksh_fUniqueID[kMaxtrksh];   //[trksh_]
   UInt_t          trksh_fBits[kMaxtrksh];   //[trksh_]
   Int_t           trksh_sid[kMaxtrksh];   //[trksh_]
   Int_t           trksh_zface[kMaxtrksh];   //[trksh_]
   Int_t           trksh_mnid[kMaxtrksh];   //[trksh_]
   Float_t         trksh_time[kMaxtrksh];   //[trksh_]
   Float_t         trksh_dt[kMaxtrksh];   //[trksh_]
   Float_t         trksh_tot0[kMaxtrksh];   //[trksh_]
   Float_t         trksh_tot1[kMaxtrksh];   //[trksh_]
   Float_t         trksh_edep[kMaxtrksh];   //[trksh_]
   Int_t           trksh_iseg[kMaxtrksh];   //[trksh_]
   Int_t           trksh_itrk[kMaxtrksh];   //[trksh_]
   Int_t           trksh_ihit[kMaxtrksh];   //[trksh_]
   Float_t         trksh_rdrift[kMaxtrksh];   //[trksh_]
   Float_t         trksh_doca[kMaxtrksh];   //[trksh_]
   Float_t         trksh_dr[kMaxtrksh];   //[trksh_]
   Float_t         trksh_drho[kMaxtrksh];   //[trksh_]
   Int_t           nseg;
   Int_t           seg_;
   UInt_t          seg_fUniqueID[kMaxseg];   //[seg_]
   UInt_t          seg_fBits[kMaxseg];   //[seg_]
   Int_t           seg_sid[kMaxseg];   //[seg_]
   Int_t           seg_nh[kMaxseg];   //[seg_]
   Int_t           seg_ngh[kMaxseg];   //[seg_]
   Int_t           seg_nghl[kMaxseg][2];   //[seg_]
   Int_t           seg_nmhl[kMaxseg][2];   //[seg_]
   Int_t           seg_ntrans[kMaxseg];   //[seg_]
   Float_t         seg_t0[kMaxseg];   //[seg_]
   Float_t         seg_chi2d[kMaxseg];   //[seg_]
   Float_t         seg_y0[kMaxseg];   //[seg_]
   Float_t         seg_z0[kMaxseg];   //[seg_]
   Float_t         seg_ymean[kMaxseg];   //[seg_]
   Float_t         seg_dzdy[kMaxseg];   //[seg_]
   Float_t         seg_y0t[kMaxseg];   //[seg_]
   Float_t         seg_dzdyt[kMaxseg];   //[seg_]
   Int_t           nsegsh;
   Int_t           segsh_;
   UInt_t          segsh_fUniqueID[kMaxsegsh];   //[segsh_]
   UInt_t          segsh_fBits[kMaxsegsh];   //[segsh_]
   Int_t           segsh_sid[kMaxsegsh];   //[segsh_]
   Int_t           segsh_zface[kMaxsegsh];   //[segsh_]
   Int_t           segsh_mnid[kMaxsegsh];   //[segsh_]
   Float_t         segsh_time[kMaxsegsh];   //[segsh_]
   Float_t         segsh_dt[kMaxsegsh];   //[segsh_]
   Float_t         segsh_tot0[kMaxsegsh];   //[segsh_]
   Float_t         segsh_tot1[kMaxsegsh];   //[segsh_]
   Float_t         segsh_edep[kMaxsegsh];   //[segsh_]
   Int_t           segsh_iseg[kMaxsegsh];   //[segsh_]
   Int_t           segsh_itrk[kMaxsegsh];   //[segsh_]
   Int_t           segsh_ihit[kMaxsegsh];   //[segsh_]
   Float_t         segsh_rdrift[kMaxsegsh];   //[segsh_]
   Float_t         segsh_doca[kMaxsegsh];   //[segsh_]
   Float_t         segsh_dr[kMaxsegsh];   //[segsh_]
   Float_t         segsh_drho[kMaxsegsh];   //[segsh_]

   // List of branches
   TBranch        *b_evt_run;   //!
   TBranch        *b_evt_srn;   //!
   TBranch        *b_evt_evn;   //!
   TBranch        *b_evt_nsdtot;   //!
   TBranch        *b_evt_sd_;   //!
   TBranch        *b_sd_fUniqueID;   //!
   TBranch        *b_sd_fBits;   //!
   TBranch        *b_sd_sid;   //!
   TBranch        *b_sd_mnid;   //!
   TBranch        *b_sd_tdc0;   //!
   TBranch        *b_sd_tdc1;   //!
   TBranch        *b_sd_tot0;   //!
   TBranch        *b_sd_tot1;   //!
   TBranch        *b_sd_pmp;   //!
   TBranch        *b_sd_flag;   //!
   TBranch        *b_sd_fs;   //!
   TBranch        *b_sd_bl;   //!
   TBranch        *b_sd_ph;   //!
   TBranch        *b_sd_ns;   //!
   TBranch        *b_sd_adc;   //!
   TBranch        *b_evt_nshtot;   //!
   TBranch        *b_evt_nsh;   //!
   TBranch        *b_evt_pmp;   //!
   TBranch        *b_evt_sh_;   //!
   TBranch        *b_sh_fUniqueID;   //!
   TBranch        *b_sh_fBits;   //!
   TBranch        *b_sh_sid;   //!
   TBranch        *b_sh_zface;   //!
   TBranch        *b_sh_mnid;   //!
   TBranch        *b_sh_time;   //!
   TBranch        *b_sh_dt;   //!
   TBranch        *b_sh_tot0;   //!
   TBranch        *b_sh_tot1;   //!
   TBranch        *b_sh_edep;   //!
   TBranch        *b_evt_maxEdep;   //!
   TBranch        *b_evt_nch;   //!
   TBranch        *b_evt_ch_;   //!
   TBranch        *b_ch_fUniqueID;   //!
   TBranch        *b_ch_fBits;   //!
   TBranch        *b_ch_sid;   //!
   TBranch        *b_ch_nsh;   //!
   TBranch        *b_ch_zface;   //!
   TBranch        *b_ch_mnid;   //!
   TBranch        *b_ch_time;   //!
   TBranch        *b_ch_dtime;   //!
   TBranch        *b_ch_x;   //!
   TBranch        *b_ch_y;   //!
   TBranch        *b_ch_z;   //!
   TBranch        *b_ch_ux;   //!
   TBranch        *b_ch_uy;   //!
   TBranch        *b_ch_ures;   //!
   TBranch        *b_ch_vres;   //!
   TBranch        *b_ch_edep;   //!
   TBranch        *b_evt_ntc;   //!
   TBranch        *b_evt_tc_;   //!
   TBranch        *b_tc_fUniqueID;   //!
   TBranch        *b_tc_fBits;   //!
   TBranch        *b_tc_nsh;   //!
   TBranch        *b_tc_nch;   //!
   TBranch        *b_tc_t0;   //!
   TBranch        *b_tc_tmin;   //!
   TBranch        *b_tc_tmax;   //!
   TBranch        *b_tc_edep_max;   //!
   TBranch        *b_tc_y0;   //!
   TBranch        *b_tc_dydz;   //!
   TBranch        *b_tc_chi2yz;   //!
   TBranch        *b_tc_nplanes;   //!
   TBranch        *b_tc_nfaces;   //!
   TBranch        *b_tc_npanels;   //!
   TBranch        *b_tc__nhf;   //!
   TBranch        *b_tc__timef;   //!
   TBranch        *b_tc__nhp;   //!
   TBranch        *b_tc__timep;   //!
   TBranch        *b_tc__mnid;   //!
   TBranch        *b_tc__nh_panel;   //!
   TBranch        *b_tc__time_panel;   //!
   TBranch        *b_tc__edep_panel;   //!
   TBranch        *b_tc_max_nh_panel;   //!
   TBranch        *b_evt_ntrk;   //!
   TBranch        *b_evt_trk_;   //!
   TBranch        *b_trk_fUniqueID;   //!
   TBranch        *b_trk_fBits;   //!
   TBranch        *b_trk_nhits;   //!
   TBranch        *b_trk_t0;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_evt_ntrksh;   //!
   TBranch        *b_evt_trksh_;   //!
   TBranch        *b_trksh_fUniqueID;   //!
   TBranch        *b_trksh_fBits;   //!
   TBranch        *b_trksh_sid;   //!
   TBranch        *b_trksh_zface;   //!
   TBranch        *b_trksh_mnid;   //!
   TBranch        *b_trksh_time;   //!
   TBranch        *b_trksh_dt;   //!
   TBranch        *b_trksh_tot0;   //!
   TBranch        *b_trksh_tot1;   //!
   TBranch        *b_trksh_edep;   //!
   TBranch        *b_trksh_iseg;   //!
   TBranch        *b_trksh_itrk;   //!
   TBranch        *b_trksh_ihit;   //!
   TBranch        *b_trksh_rdrift;   //!
   TBranch        *b_trksh_doca;   //!
   TBranch        *b_trksh_dr;   //!
   TBranch        *b_trksh_drho;   //!
   TBranch        *b_evt_nseg;   //!
   TBranch        *b_evt_seg_;   //!
   TBranch        *b_seg_fUniqueID;   //!
   TBranch        *b_seg_fBits;   //!
   TBranch        *b_seg_sid;   //!
   TBranch        *b_seg_nh;   //!
   TBranch        *b_seg_ngh;   //!
   TBranch        *b_seg_nghl;   //!
   TBranch        *b_seg_nmhl;   //!
   TBranch        *b_seg_ntrans;   //!
   TBranch        *b_seg_t0;   //!
   TBranch        *b_seg_chi2d;   //!
   TBranch        *b_seg_y0;   //!
   TBranch        *b_seg_z0;   //!
   TBranch        *b_seg_ymean;   //!
   TBranch        *b_seg_dzdy;   //!
   TBranch        *b_seg_y0t;   //!
   TBranch        *b_seg_dzdyt;   //!
   TBranch        *b_evt_nsegsh;   //!
   TBranch        *b_evt_segsh_;   //!
   TBranch        *b_segsh_fUniqueID;   //!
   TBranch        *b_segsh_fBits;   //!
   TBranch        *b_segsh_sid;   //!
   TBranch        *b_segsh_zface;   //!
   TBranch        *b_segsh_mnid;   //!
   TBranch        *b_segsh_time;   //!
   TBranch        *b_segsh_dt;   //!
   TBranch        *b_segsh_tot0;   //!
   TBranch        *b_segsh_tot1;   //!
   TBranch        *b_segsh_edep;   //!
   TBranch        *b_segsh_iseg;   //!
   TBranch        *b_segsh_itrk;   //!
   TBranch        *b_segsh_ihit;   //!
   TBranch        *b_segsh_rdrift;   //!
   TBranch        *b_segsh_doca;   //!
   TBranch        *b_segsh_dr;   //!
   TBranch        *b_segsh_drho;   //!

  plot_n001_hist(int RunNumber);
  plot_n001_hist(const char* Fn = "nts.mu2e.trk.vst00s000r000n001.120718_000001.root");
  virtual ~plot_n001_hist();
  
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual bool     Notify();
  virtual void     Show(Long64_t entry = -1);

  void             Loop          (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents, double DtMin = 0, double DtMax = -1);
  void             LoopScan      (int Mnid1, int Ch1, int Mnid2, int Ch2, int NEvents = -1);
  int              BookHistograms(Hist_t* Hist, int Mnid1, int Ch1, int Mnid2, int Ch2, float DtMin, float DtMax);
};

#endif
