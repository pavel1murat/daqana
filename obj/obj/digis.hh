//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 21 16:42:59 2025 by ROOT version 6.32.06
// from TTree digis/digis
// found on file: nts.mu2e.trk.vst00s000r011n003.107995_000001.root
//////////////////////////////////////////////////////////

#ifndef digis_h
#define digis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFolder.h>

// Header file for the classes stored in the TTree if any.
#include "DaqEvent.hh"
#include "TObject.h"
#include "DaqStrawDigi.hh"
#include "DaqStrawHit.hh"
#include "DaqComboHit.hh"
#include "DaqTimeCluster.hh"
#include "DaqTrack.hh"
#include "DaqTrkStrawHit.hh"
#include "DaqSegment.hh"

class digis {
public :

  enum {
    kNEventHistSets   = 100,
    kNSegmentHistSets = 100,
    kNTwoSegHistSets  =  10,
    kNSshtHistSets    = 100,
  };

  struct TwoSegPar_t {
    int   hid;
    float dt;
    int   diff_planes;
  };

  struct EventHist_t {
    TH1F*     fNSeg[2];
  };

  struct SegmentHist_t {
    TH1F* fNHits;
    TH1F* fNTrans;
    TH1F* fChi2d;
    TH1F* fDzDy;
    TH1F* fY0;
    TH2F* fNghl;
    TH2F* fNmhl;
  };

  struct TwoSegHist_t {
    TH1F* fDt[1100];
  };

  struct SshtHist_t {
    TH1F* fRDrift;
    TH1F* fDr;
    TH1F* fDrho;
    TH2F* fDrVsRDrift;
    TH2F* fDrhoVsRDrift;
    TH2F* fDrhoVsStraw;
    TH2F* fDrVsStraw;
  };

  struct Hist_t {
    EventHist_t*    fEvent  [kNEventHistSets  ];
    SegmentHist_t*  fSegment[kNSegmentHistSets];
    SshtHist_t*     fSsht   [kNSshtHistSets   ];
    TwoSegHist_t*   fTwoSeg [kNTwoSegHistSets ];
  };

  Hist_t          fHist;

  int             fNSeg6;
  int             fISeg6[100];   // index of a 6+ hit segments

                                        // to be initialized at start up
  int             fMinSegNHits;
  int             fMinSegNghl;
  float           fMaxSegChi2d;

  TFolder*        fFolder;
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any... mmm
// PM: this - max dimensions - is only true for a given file...

   static constexpr Int_t kMaxsd    = 1000;
   static constexpr Int_t kMaxsh    = 1000;
   static constexpr Int_t kMaxch    = 1000;
   static constexpr Int_t kMaxtc    = 100;
   static constexpr Int_t kMaxtrk   = 100;
   static constexpr Int_t kMaxtrksh = 1000;
   static constexpr Int_t kMaxseg   = 100;
   static constexpr Int_t kMaxsegsh = 1000;

   // Declaration of leaf types
 //DaqEvent        *evt;
   Int_t           run;
   Int_t           srn;
   Int_t           evn;
   Int_t           nsdtot;
   Short_t         nsd[2][6][96];
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
   Int_t           nsh[2][6];
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
   Int_t           tc__nhf[kMaxtc][4];   //[tc_]
   Float_t         tc__timef[kMaxtc][4];   //[tc_]
   Int_t           tc__nhp[kMaxtc][2];   //[tc_]
   Float_t         tc__timep[kMaxtc][2];   //[tc_]
   Int_t           tc__mnid[kMaxtc][12];   //[tc_]
   Int_t           tc__nh_panel[kMaxtc][12];   //[tc_]
   Float_t         tc__time_panel[kMaxtc][12];   //[tc_]
   Float_t         tc__edep_panel[kMaxtc][12];   //[tc_]
   Int_t           tc_max_nh_panel[kMaxtc];   //[tc_]
//-----------------------------------------------------------------------------
// tracks
//-----------------------------------------------------------------------------
   Int_t           ntrk;
   Int_t           trk_;
   UInt_t          trk_fUniqueID[kMaxtrk];   //[trk_]
   UInt_t          trk_fBits[kMaxtrk];   //[trk_]
   Int_t           trk_nhits[kMaxtrk];   //[trk_]
   Float_t         trk_t0[kMaxtrk];   //[trk_]
   Float_t         trk_chi2[kMaxtrk];   //[trk_]
   Int_t           ntrksh;
   Int_t           trksh_;
   DaqTrkStrawHit  *trksh_DaqStrawHit;
   Float_t         trksh_rdrift[kMaxtrksh];   //[trksh_]
   Float_t         trksh_doca[kMaxtrksh];   //[trksh_]
   Int_t           trksh_drho[kMaxtrksh];   //[trksh_]
   Int_t           trksh_iseg[kMaxtrksh];   //[trksh_]
   Int_t           trksh_itrk[kMaxtrksh];   //[trksh_]
//-----------------------------------------------------------------------------
// segments
//-----------------------------------------------------------------------------
   Int_t           nseg;
   Int_t           seg_;
   UInt_t          seg_fUniqueID[kMaxseg];   //[seg_]
   UInt_t          seg_fBits[kMaxseg];   //[seg_]
   Int_t           seg_sid[kMaxseg];   //[seg_]
   Int_t           seg_nh[kMaxseg];   //[seg_]
   Int_t           seg_ntrans[kMaxseg];   //[seg_]
   Int_t           seg_ngh[kMaxseg];   //[seg_]
   Int_t           seg_nghl[kMaxseg][2];   //[seg_]
   Int_t           seg_nmhl[kMaxseg][2];   //[seg_]
   Float_t         seg_t0[kMaxseg];   //[seg_]
   Float_t         seg_chi2d[kMaxseg];   //[seg_]
   Float_t         seg_y0[kMaxseg];   //[seg_]
   Float_t         seg_z0[kMaxseg];   //[seg_]
   Float_t         seg_ymean[kMaxseg];   //[seg_]
   Float_t         seg_dzdy[kMaxseg];   //[seg_]
   Float_t         seg_y0t[kMaxseg];   //[seg_]
   Float_t         seg_dzdyt[kMaxseg];   //[seg_]
//-----------------------------------------------------------------------------
// segment straw hits (SSHT)
//-----------------------------------------------------------------------------
   Int_t           nsegsh;
   Int_t           segsh_;
   UInt_t          segsh_fUniqueID[kMaxsh];   //[sh_]
   UInt_t          segsh_fBits[kMaxsh];   //[sh_]
   Int_t           segsh_sid[kMaxsh];   //[sh_]
   Int_t           segsh_zface[kMaxsh];   //[sh_]
   Int_t           segsh_mnid[kMaxsh];   //[sh_]
   Float_t         segsh_time[kMaxsh];   //[sh_]
   Float_t         segsh_dt[kMaxsh];   //[sh_]
   Float_t         segsh_tot0[kMaxsh];   //[sh_]
   Float_t         segsh_tot1[kMaxsh];   //[sh_]
   Float_t         segsh_edep[kMaxsh];   //[sh_]
   Float_t         segsh_rdrift[kMaxsegsh];   //[segsh_]
   Float_t         segsh_doca[kMaxsegsh];   //[segsh_]
   Float_t         segsh_drho[kMaxsegsh];   //[segsh_]
   Int_t           segsh_iseg[kMaxsegsh];   //[segsh_]
   Int_t           segsh_itrk[kMaxsegsh];   //[segsh_]

   // List of branches
   TBranch        *b_evt_run;   //!
   TBranch        *b_evt_srn;   //!
   TBranch        *b_evt_evn;   //!
   TBranch        *b_evt_nsdtot;   //!
   TBranch        *b_evt_nsd;   //!
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
   TBranch        *b_trksh_DaqStrawHit;   //!
   TBranch        *b_trksh_rdrift;   //!
   TBranch        *b_trksh_doca;   //!
   TBranch        *b_trksh_drho;   //!
   TBranch        *b_trksh_iseg;   //!
   TBranch        *b_trksh_itrk;   //!
   TBranch        *b_evt_nseg;   //!
   TBranch        *b_evt_seg_;   //!
   TBranch        *b_seg_fUniqueID;   //!
   TBranch        *b_seg_fBits;   //!
   TBranch        *b_seg_sid;   //!
   TBranch        *b_seg_nh;   //!
   TBranch        *b_seg_ntrans;   //!
   TBranch        *b_seg_ngh;   //!
   TBranch        *b_seg_nghl;   //!
   TBranch        *b_seg_nmhl;   //!
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
  //   TBranch        *b_segsh_DaqStrawHit;   //!
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
   TBranch        *b_segsh_rdrift;   //!
   TBranch        *b_segsh_doca;   //!
   TBranch        *b_segsh_drho;   //!
   TBranch        *b_segsh_iseg;   //!
   TBranch        *b_segsh_itrk;   //!

   digis(const char* Fn=0);
   virtual ~digis();

  int CalculateMissingParameters();

  int FillTwoSegHistograms (TwoSegHist_t*  Hist, TwoSegPar_t* Par);
  int FillSegmentHistograms(SegmentHist_t* Hist, int I);
  int FillSshtHistograms   (SshtHist_t*    Hist, int IHit);
  int FillEventHistograms  (EventHist_t*   Hist);
  int FillHistograms       (Hist_t*        Hist);

  int BookTwoSegHistograms (TwoSegHist_t*  Hist, const char* Folder);
  int BookSshtHistograms   (SshtHist_t*    Hist, const char* Folder);
  int BookSegmentHistograms(SegmentHist_t* Hist, const char* Folder);
  int BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
  int BookHistograms       (Hist_t*        Hist);

  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual bool     Notify();
  virtual void     Show(Long64_t entry = -1);

  int panel(int Sid) { return ((Sid >>  7) & 0x07) ; } 
  int plane(int Sid) { return ((Sid >> 10) & 0x3f) ; } 
//-----------------------------------------------------------------------------
// the following helper methods allow to save 1 line per request, which in 
// case of 100's histograms booked is a non-negligible number
//-----------------------------------------------------------------------------
  void  DeleteHistograms(TFolder* Folder = (TFolder*) -1);

  void  AddHistogram(TObject* hist, const char* FolderName = "Hist");

  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, const float* LowEdge, 
		const char* FolderName = "Hist");

  void  HBook1D(TH1D*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook1D(TH1D*& Hist, const char* Name, const char* Title,
		Int_t Nx, const double* LowEdge, 
		const char* FolderName = "Hist");

  void  HBook2F(TH2F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		Int_t Ny, Double_t YMin, Double_t YMax,
		const char* FolderName = "Hist");

  void  HProf (TProfile*& Hist, const char* Name, const char* Title,
	       Int_t Nx, Double_t XMin, Double_t XMax,
	       Double_t YMin, Double_t YMax,
	       const char* FolderName = "Hist");

  int  SaveFolder(TFolder* Folder, TDirectory* Dir);
  void SaveHist  (const char* Filename);
};

#endif

#ifdef digis_cxx
#endif // #ifdef digis_cxx
