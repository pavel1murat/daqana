#define digis_cxx
#include "daqana/obj/digis.hh"
#include <TH2.h>
// #include <TStyle.h>
// #include <TCanvas.h>
#include <TROOT.h>
#include <iostream>
#include <format>

//-----------------------------------------------------------------------------
void digis::Loop() {
//   In a ROOT session, you can do:
//      root> .L digis.C
//      root> digis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     CalculateMissingParameters();
     FillHistograms(&fHist);
   }
}




//-----------------------------------------------------------------------------
int digis::CalculateMissingParameters() {

  fNSeg6 = 0;
                                        // seelct "good" segments
  for (int i=0; i<nseg; i++) {
    if ((seg_ngh[i] >= fMinSegNHits) and (seg_chi2d[i] < fMaxSegChi2d) and
        (seg_nghl[i][0] >= fMinSegNghl) and  (seg_nghl[i][1] >= fMinSegNghl)) {
      fISeg6[fNSeg6] = i;
      fNSeg6        += 1;
    }
  }

  return 0;
};

//-----------------------------------------------------------------------------
int digis::FillTwoSegHistograms(TwoSegHist_t* Hist, TwoSegPar_t* Par) {

  printf("Par->hid:%i Par->dt: %10.3f hist:0x%8p\n",Par->hid,Par->dt,(void*) Hist->fDt[Par->hid]);

  Hist->fDt[Par->hid]->Fill(Par->dt);

  if (Par->diff_planes == 1) {
    Hist->fDt[0]->Fill(Par->dt);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int digis::FillSegmentHistograms(SegmentHist_t* Hist, int I) {
  Hist->fNHits->Fill(seg_nh   [I]);
  Hist->fNTrans->Fill(seg_ntrans   [I]);
  Hist->fChi2d->Fill(seg_chi2d[I]);
  Hist->fDzDy->Fill(seg_dzdy[I]);
  Hist->fY0->Fill(seg_y0[I]);
  Hist->fNghl->Fill (seg_nghl[I][0],seg_nghl[I][1]);
  Hist->fNmhl->Fill (seg_nmhl[I][0],seg_nmhl[I][1]);
  return 0;
}

//-----------------------------------------------------------------------------
int digis::FillEventHistograms(EventHist_t* Hist) {
  Hist->fNSeg[0]->Fill(nseg  );
  Hist->fNSeg[1]->Fill(fNSeg6);

  return 0;
}

//-----------------------------------------------------------------------------
int digis::FillSshtHistograms(SshtHist_t* Hist, int I) {
  Hist->fRDrift->Fill(segsh_rdrift[I]);
  Hist->fDr->Fill(segsh_doca[I]);

  Hist->fDrho->Fill(segsh_drho[I]);
  Hist->fDrVsRDrift->Fill(segsh_rdrift[I],segsh_doca[I]);
  Hist->fDrhoVsRDrift->Fill(segsh_rdrift[I],segsh_drho[I]);

  int straw = (segsh_sid[I] &0x7f);
  Hist->fDrVsStraw->Fill(straw,segsh_doca[I]);
  Hist->fDrhoVsStraw->Fill(straw,segsh_drho[I]);

  return 0;
}

//-----------------------------------------------------------------------------
int digis::FillHistograms(Hist_t* Hist) {
                                        // fill event histograms

  FillEventHistograms(Hist->fEvent[0]);
  if (fNSeg6 >= 2) FillEventHistograms(Hist->fEvent[1]);

                                        // segment histograms
  for (int i=0; i<nseg; i++) {
    int panel = (seg_sid[i] >>  7) & 0x7;
    int plane = (seg_sid[i] >> 10) & 0x3f;
    int unique_panel = 6*plane+panel;

    FillSegmentHistograms(Hist->fSegment[0],i);
    if (seg_nh[i] >= fMinSegNHits) {
      FillSegmentHistograms(Hist->fSegment[1],i);
      FillSegmentHistograms(Hist->fSegment[unique_panel+10],i);
      if ((seg_nghl[i][0] > 0) and (seg_nghl[i][1] > 0)) {
        FillSegmentHistograms(Hist->fSegment[2],i);
        if (seg_chi2d[i] <fMaxSegChi2d) {
          FillSegmentHistograms(Hist->fSegment[3],i);
          if (fNSeg6 >= 2) {
            FillSegmentHistograms(Hist->fSegment[4],i);
          }
        }
      }
    }
  }
//-----------------------------------------------------------------------------  
// dt between two 6+ chi2d<10 segments
//-----------------------------------------------------------------------------
  TwoSegPar_t par_ss;

  for (int i1=0; i1<fNSeg6-1; i1++) {
    int iseg1 = fISeg6[i1];
    int plane1 = plane(seg_sid[iseg1]);
    int panel1 = panel(seg_sid[iseg1]);
    int ip1    = 6*plane1+panel1;
    for (int i2=i1+1; i2<fNSeg6; i2++) {
      int iseg2  = fISeg6[i2];
      int plane2 = plane(seg_sid[iseg2]);
      int panel2 = panel(seg_sid[iseg2]);
      int ip2    = 6*plane2+panel2;
//-----------------------------------------------------------------------------
// two "good" segments, figure the histogram id
//-----------------------------------------------------------------------------
      // int   hid;
      // float dt;
      if (ip1 < ip2) {
        par_ss.hid = 100*ip1+ip2;
        par_ss.dt  = seg_t0[iseg1]-seg_t0[iseg2];
      }
      else {
        par_ss.hid = 100*ip2+ip1;
        par_ss.dt = seg_t0[iseg2]-seg_t0[iseg1];
      }
      if (plane1 == plane2) par_ss.diff_planes = 0;
      else                  par_ss.diff_planes = 1;

      FillTwoSegHistograms(Hist->fTwoSeg[0],&par_ss);
      if (seg_dzdy[i1] > 0) FillTwoSegHistograms(Hist->fTwoSeg[1],&par_ss);
      else                  FillTwoSegHistograms(Hist->fTwoSeg[2],&par_ss);
    }
  }
//-----------------------------------------------------------------------------
// segment straw hit histograms (SSHT) ... do that only for good segments
//-----------------------------------------------------------------------------
  for (int is=0; is<fNSeg6; is++) {
    int iseg = fISeg6[is];
    for (int ish=0; ish<nsegsh; ish++) {
      if (segsh_iseg[ish] != iseg)         continue;
      int plane = (segsh_sid[ish] >> 10) & 0x3f;
      int panel = (segsh_sid[ish] >>  7) & 0x7;
      // int straw = (segsh_sid[ish] >>  0) & 0x7f;
//-----------------------------------------------------------------------------
// all together
//-----------------------------------------------------------------------------
      FillSshtHistograms(Hist->fSsht[0]   , ish);
//-----------------------------------------------------------------------------
// each panel separately
//-----------------------------------------------------------------------------
      int iset = 10+plane*6+panel;
      FillSshtHistograms(Hist->fSsht[iset], ish);
    }
  }
  return 0;
};

//-----------------------------------------------------------------------------
int digis::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNSeg[0],"nseg_0","nseg_0",20,0,20, Folder);
  HBook1F(Hist->fNSeg[1],"nseg_1","nseg_1",20,0,20, Folder);
  return 0;
}

//-----------------------------------------------------------------------------
int digis::BookSegmentHistograms(SegmentHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNHits  ,"nhits"  ,"nhits", 100,    0, 100, Folder);
  HBook1F(Hist->fNTrans ,"ntrans"  ,"ntrans", 10,    0, 10, Folder);
  HBook2F(Hist->fNghl   ,"nghl"   ,"nghl" ,  50,    0,  50, 50,0,50,Folder);
  HBook2F(Hist->fNmhl   ,"nmhl"   ,"nmhl" ,  50,    0,  50, 50,0,50,Folder);
  HBook1F(Hist->fChi2d  ,"chi2d"  ,"chi2d", 100,    0,  50, Folder);
  HBook1F(Hist->fDzDy   ,"dzdy"   ,"dzdy" , 200,   -1,   1, Folder);
  HBook1F(Hist->fY0     ,"y0"     ,"y0"   , 100,    0, 100, Folder);
  return 0;
}
//-----------------------------------------------------------------------------
int digis::BookTwoSegHistograms(TwoSegHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fDt[0],Form("dt_%03i",0),Form("dt_%03i",0), 100,-50,50,Folder);
  
  for (int i1=0; i1<11; i1++) {
    for (int i2=i1+1; i2<12; i2++) {
      int hid = 100*i1+i2;
      HBook1F(Hist->fDt[hid],Form("dt_%03i",hid),Form("dt_%03i",hid), 100,-50,50,Folder);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// add more descriptive titles later
//-----------------------------------------------------------------------------
int digis::BookSshtHistograms(SshtHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fRDrift  ,"rdrift"  ,"rdrift, mm", 100,    0, 5, Folder);
  HBook1F(Hist->fDr ,"dr"  ,"dr, mm", 500,    -2.5, 2.5, Folder);
  HBook1F(Hist->fDrho ,"drho"  ,"drho, mm", 500,    -2.5, 2.5, Folder);
  HBook2F(Hist->fDrVsRDrift   ,"dr_vs_rdrift"   ,"dr_vs_rdrift" ,  100,    0,  5, 250,-2.5,2.5,Folder);
  HBook2F(Hist->fDrhoVsRDrift   ,"drho_vs_rdrift"   ,"drho_vs_rdrift" ,  100,    0,  5, 250,-2.5,2.5,Folder);
  HBook2F(Hist->fDrVsStraw   ,"dr_vs_straw"   ,"dr_vs_straw" ,  100,    0,  100, 250,-2.5,2.5,Folder);
  HBook2F(Hist->fDrhoVsStraw ,"drho_vs_straw"   ,"drho_vs_straw" ,  100,    0,  100, 250,-2.5,2.5,Folder);
  return 0;
}


//-----------------------------------------------------------------------------
int digis::BookHistograms(Hist_t* Hist) {
  TH1::AddDirectory(0);
  
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with 2 6+ segments

  TFolder* fol;
  char folder_name[64];

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_segment_histset[kNSegmentHistSets];
  for (int i=0; i<kNSegmentHistSets; i++) book_segment_histset[i] = 0;

  book_segment_histset[ 0] = 1;		// all segments
  book_segment_histset[ 1] = 1;	        // segments with 4+ good hits
  book_segment_histset[ 2] = 1;	        // segments with 4+ good hits and chi2d<10
  book_segment_histset[ 3] = 1;	        // segments with 4+ good hits and chi2d<10 and nghl[i] > 0
  book_segment_histset[ 4] = 1;	        // events with 2 such segments

  book_segment_histset[10] = 1;	        // events with segmentsin panel 0 
  book_segment_histset[11] = 1;	        // events with segmentsin panel 1
  book_segment_histset[12] = 1;	        // events with segmentsin panel 2
  book_segment_histset[13] = 1;	        // events with segmentsin panel 3
  book_segment_histset[14] = 1;	        // events with segmentsin panel 4
  book_segment_histset[15] = 1;	        // events with segmentsin panel 5
  book_segment_histset[16] = 1;	        // events with segmentsin panel 6
  book_segment_histset[17] = 1;	        // events with segmentsin panel 7
  book_segment_histset[18] = 1;	        // events with segmentsin panel 8
  book_segment_histset[19] = 1;	        // events with segmentsin panel 9
  book_segment_histset[20] = 1;	        // events with segmentsin panel 10
  book_segment_histset[21] = 1;	        // events with segmentsin panel 11

  for (int i=0; i<kNSegmentHistSets; i++) {
    if (book_segment_histset[i] != 0) {
      sprintf(folder_name,"seg_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fSegment[i] = new SegmentHist_t;
      BookSegmentHistograms(fHist.fSegment[i],Form("%s",folder_name));
    }
  }

//-----------------------------------------------------------------------------
// book 2-seg histograms
//-----------------------------------------------------------------------------
  int book_twoseg_histset[kNTwoSegHistSets];
  for (int i=0; i<kNTwoSegHistSets; i++) book_twoseg_histset[i] = 0;

  book_twoseg_histset[ 0] = 1;		// all twosegs
  book_twoseg_histset[ 1] = 1;		// dzdy(p1) > 0
  book_twoseg_histset[ 2] = 1;		// dzdy(p1) < 0

  for (int i=0; i<kNTwoSegHistSets; i++) {
    if (book_twoseg_histset[i] != 0) {
      sprintf(folder_name,"s2_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fTwoSeg[i] = new TwoSegHist_t;
      BookTwoSegHistograms(fHist.fTwoSeg[i],Form("%s",folder_name));
    }
  }

//-----------------------------------------------------------------------------
// book segment hit (ssht) histograms
//-----------------------------------------------------------------------------
  int book_ssht_histset[kNSshtHistSets];
  for (int i=0; i<kNSshtHistSets; i++) book_ssht_histset[i] = 0;

  book_ssht_histset[ 0] = 1;		// all segments

  book_ssht_histset[10] = 1;	        // segments with 4+ good hits panel#0  (+10)
  book_ssht_histset[11] = 1;	        // segments with 4+ good hits panel#1
  book_ssht_histset[12] = 1;	        // segments with 4+ good hits panel#2
  book_ssht_histset[13] = 1;	        // segments with 4+ good hits panel#3
  book_ssht_histset[14] = 1;	        // segments with 4+ good hits panel#4
  book_ssht_histset[15] = 1;	        // segments with 4+ good hits panel#5
  book_ssht_histset[16] = 1;	        // segments with 4+ good hits panel#6
  book_ssht_histset[17] = 1;	        // segments with 4+ good hits panel#7
  book_ssht_histset[18] = 1;	        // segments with 4+ good hits panel#8
  book_ssht_histset[19] = 1;	        // segments with 4+ good hits panel#9
  book_ssht_histset[20] = 1;	        // segments with 4+ good hits panel#10
  book_ssht_histset[21] = 1;	        // segments with 4+ good hits panel#10

  for (int i=0; i<kNSshtHistSets; i++) {
    if (book_ssht_histset[i] != 0) {
      sprintf(folder_name,"ssht_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fSsht[i] = new SshtHist_t;
      BookSshtHistograms(fHist.fSsht[i],Form("%s",folder_name));
    }
  }

  return 0;
};


//-----------------------------------------------------------------------------
digis::digis(const char* Fn) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  std::string filename;
  if (Fn) filename = Fn;
  else    filename = "nts/nts.mu2e.trk.vst00s000r011n003.107995_000001.root";

  TFile* f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
  if (!f || !f->IsOpen()) {
    f = new TFile(filename.data());
  }

  TTree* tree = (TTree*) f->Get("/MakeDigiNtuple/digis");

  fMinSegNHits =  4 ;
  fMinSegNghl  =  1;
  fMaxSegChi2d = 10.;

  fFolder = gROOT->GetRootFolder()->AddFolder("digis","digis");

  Init(tree);

  BookHistograms(&fHist);
}

//-----------------------------------------------------------------------------
digis::~digis() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
Int_t digis::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

//-----------------------------------------------------------------------------
Long64_t digis::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

//-----------------------------------------------------------------------------
void digis::Init(TTree *tree) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   // trksh_DaqStrawHit = 0;
   // segsh_DaqStrawHit = 0;
   // Set array pointer
   for(int i=0; i<kMaxsd; ++i) sd_adc[i] = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_evt_run);
   fChain->SetBranchAddress("srn", &srn, &b_evt_srn);
   fChain->SetBranchAddress("evn", &evn, &b_evt_evn);
   fChain->SetBranchAddress("nsdtot", &nsdtot, &b_evt_nsdtot);
   fChain->SetBranchAddress("nsd[2][6][96]", nsd, &b_evt_nsd);
   fChain->SetBranchAddress("sd", &sd_, &b_evt_sd_);
   fChain->SetBranchAddress("sd.fUniqueID", &sd_fUniqueID, &b_sd_fUniqueID);
   fChain->SetBranchAddress("sd.fBits", &sd_fBits, &b_sd_fBits);
   fChain->SetBranchAddress("sd.sid", &sd_sid, &b_sd_sid);
   fChain->SetBranchAddress("sd.mnid", &sd_mnid, &b_sd_mnid);
   fChain->SetBranchAddress("sd.tdc0", &sd_tdc0, &b_sd_tdc0);
   fChain->SetBranchAddress("sd.tdc1", &sd_tdc1, &b_sd_tdc1);
   fChain->SetBranchAddress("sd.tot0", &sd_tot0, &b_sd_tot0);
   fChain->SetBranchAddress("sd.tot1", &sd_tot1, &b_sd_tot1);
   fChain->SetBranchAddress("sd.pmp", &sd_pmp, &b_sd_pmp);
   fChain->SetBranchAddress("sd.flag", &sd_flag, &b_sd_flag);
   fChain->SetBranchAddress("sd.fs", &sd_fs, &b_sd_fs);
   fChain->SetBranchAddress("sd.bl", &sd_bl, &b_sd_bl);
   fChain->SetBranchAddress("sd.ph", &sd_ph, &b_sd_ph);
   fChain->SetBranchAddress("sd.ns", &sd_ns, &b_sd_ns);
   fChain->SetBranchAddress("sd.adc", &sd_adc, &b_sd_adc);
   fChain->SetBranchAddress("nshtot", &nshtot, &b_evt_nshtot);
   fChain->SetBranchAddress("nsh[2][6]", nsh, &b_evt_nsh);
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
   fChain->SetBranchAddress("nch", &nch, &b_evt_nch);
   fChain->SetBranchAddress("ch", &ch_, &b_evt_ch_);
   fChain->SetBranchAddress("ch.fUniqueID", ch_fUniqueID, &b_ch_fUniqueID);
   fChain->SetBranchAddress("ch.fBits", ch_fBits, &b_ch_fBits);
   fChain->SetBranchAddress("ch.sid", ch_sid, &b_ch_sid);
   fChain->SetBranchAddress("ch.nsh", ch_nsh, &b_ch_nsh);
   fChain->SetBranchAddress("ch.zface", ch_zface, &b_ch_zface);
   fChain->SetBranchAddress("ch.mnid", ch_mnid, &b_ch_mnid);
   fChain->SetBranchAddress("ch.time", ch_time, &b_ch_time);
   fChain->SetBranchAddress("ch.dtime", ch_dtime, &b_ch_dtime);
   fChain->SetBranchAddress("ch.x", ch_x, &b_ch_x);
   fChain->SetBranchAddress("ch.y", ch_y, &b_ch_y);
   fChain->SetBranchAddress("ch.z", ch_z, &b_ch_z);
   fChain->SetBranchAddress("ch.ux", ch_ux, &b_ch_ux);
   fChain->SetBranchAddress("ch.uy", ch_uy, &b_ch_uy);
   fChain->SetBranchAddress("ch.ures", ch_ures, &b_ch_ures);
   fChain->SetBranchAddress("ch.vres", ch_vres, &b_ch_vres);
   fChain->SetBranchAddress("ch.edep", ch_edep, &b_ch_edep);
   fChain->SetBranchAddress("ntc", &ntc, &b_evt_ntc);
   fChain->SetBranchAddress("tc", &tc_, &b_evt_tc_);
   fChain->SetBranchAddress("tc.fUniqueID", tc_fUniqueID, &b_tc_fUniqueID);
   fChain->SetBranchAddress("tc.fBits", tc_fBits, &b_tc_fBits);
   fChain->SetBranchAddress("tc.nsh", tc_nsh, &b_tc_nsh);
   fChain->SetBranchAddress("tc.nch", tc_nch, &b_tc_nch);
   fChain->SetBranchAddress("tc.t0", tc_t0, &b_tc_t0);
   fChain->SetBranchAddress("tc.tmin", tc_tmin, &b_tc_tmin);
   fChain->SetBranchAddress("tc.tmax", tc_tmax, &b_tc_tmax);
   fChain->SetBranchAddress("tc.edep_max", tc_edep_max, &b_tc_edep_max);
   fChain->SetBranchAddress("tc.y0", tc_y0, &b_tc_y0);
   fChain->SetBranchAddress("tc.dydz", tc_dydz, &b_tc_dydz);
   fChain->SetBranchAddress("tc.chi2yz", tc_chi2yz, &b_tc_chi2yz);
   fChain->SetBranchAddress("tc.nplanes", tc_nplanes, &b_tc_nplanes);
   fChain->SetBranchAddress("tc.nfaces", tc_nfaces, &b_tc_nfaces);
   fChain->SetBranchAddress("tc.npanels", tc_npanels, &b_tc_npanels);
   fChain->SetBranchAddress("tc._nhf[4]", tc__nhf, &b_tc__nhf);
   fChain->SetBranchAddress("tc._timef[4]", tc__timef, &b_tc__timef);
   fChain->SetBranchAddress("tc._nhp[2]", tc__nhp, &b_tc__nhp);
   fChain->SetBranchAddress("tc._timep[2]", tc__timep, &b_tc__timep);
   fChain->SetBranchAddress("tc._mnid[12]", tc__mnid, &b_tc__mnid);
   fChain->SetBranchAddress("tc._nh_panel[12]", tc__nh_panel, &b_tc__nh_panel);
   fChain->SetBranchAddress("tc._time_panel[12]", tc__time_panel, &b_tc__time_panel);
   fChain->SetBranchAddress("tc._edep_panel[12]", tc__edep_panel, &b_tc__edep_panel);
   fChain->SetBranchAddress("tc.max_nh_panel", tc_max_nh_panel, &b_tc_max_nh_panel);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_evt_ntrk);
   fChain->SetBranchAddress("trk", &trk_, &b_evt_trk_);
   fChain->SetBranchAddress("trk.fUniqueID", trk_fUniqueID, &b_trk_fUniqueID);
   fChain->SetBranchAddress("trk.fBits", trk_fBits, &b_trk_fBits);
   fChain->SetBranchAddress("trk.nhits", trk_nhits, &b_trk_nhits);
   fChain->SetBranchAddress("trk.t0", trk_t0, &b_trk_t0);
   fChain->SetBranchAddress("trk.chi2", trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("ntrksh", &ntrksh, &b_evt_ntrksh);
   fChain->SetBranchAddress("trksh", &trksh_, &b_evt_trksh_);
   //   fChain->SetBranchAddress("trksh.DaqStrawHit", &trksh_DaqStrawHit, &b_trksh_DaqStrawHit);
   fChain->SetBranchAddress("trksh.rdrift", &trksh_rdrift, &b_trksh_rdrift);
   fChain->SetBranchAddress("trksh.doca", &trksh_doca, &b_trksh_doca);
   fChain->SetBranchAddress("trksh.drho", &trksh_drho, &b_trksh_drho);
   fChain->SetBranchAddress("trksh.iseg", &trksh_iseg, &b_trksh_iseg);
   fChain->SetBranchAddress("trksh.itrk", &trksh_itrk, &b_trksh_itrk);
   fChain->SetBranchAddress("nseg", &nseg, &b_evt_nseg);
   fChain->SetBranchAddress("seg", &seg_, &b_evt_seg_);
   fChain->SetBranchAddress("seg.fUniqueID", seg_fUniqueID, &b_seg_fUniqueID);
   fChain->SetBranchAddress("seg.fBits", seg_fBits, &b_seg_fBits);
   fChain->SetBranchAddress("seg.sid", seg_sid, &b_seg_sid);
   fChain->SetBranchAddress("seg.nh", seg_nh, &b_seg_nh);
   fChain->SetBranchAddress("seg.ntrans", seg_ntrans, &b_seg_ntrans);
   fChain->SetBranchAddress("seg.ngh", seg_ngh, &b_seg_ngh);
   fChain->SetBranchAddress("seg.nghl[2]", seg_nghl, &b_seg_nghl);
   fChain->SetBranchAddress("seg.nmhl[2]", seg_nmhl, &b_seg_nmhl);
   fChain->SetBranchAddress("seg.t0", seg_t0, &b_seg_t0);
   fChain->SetBranchAddress("seg.chi2d", seg_chi2d, &b_seg_chi2d);
   fChain->SetBranchAddress("seg.y0", seg_y0, &b_seg_y0);
   fChain->SetBranchAddress("seg.z0", seg_z0, &b_seg_z0);
   fChain->SetBranchAddress("seg.ymean", seg_ymean, &b_seg_ymean);
   fChain->SetBranchAddress("seg.dzdy", seg_dzdy, &b_seg_dzdy);
   fChain->SetBranchAddress("seg.y0t", seg_y0t, &b_seg_y0t);
   fChain->SetBranchAddress("seg.dzdyt", seg_dzdyt, &b_seg_dzdyt);
   fChain->SetBranchAddress("nsegsh", &nsegsh, &b_evt_nsegsh);
   fChain->SetBranchAddress("segsh", &segsh_, &b_evt_segsh_);
   fChain->SetBranchAddress("segsh.fUniqueID", sh_fUniqueID, &b_sh_fUniqueID);
   fChain->SetBranchAddress("segsh.fBits", segsh_fBits, &b_sh_fBits);
   fChain->SetBranchAddress("segsh.sid", segsh_sid, &b_sh_sid);
   fChain->SetBranchAddress("segsh.zface", segsh_zface, &b_segsh_zface);
   fChain->SetBranchAddress("segsh.mnid", segsh_mnid, &b_segsh_mnid);
   fChain->SetBranchAddress("segsh.time", segsh_time, &b_segsh_time);
   fChain->SetBranchAddress("segsh.dt"  , segsh_dt, &b_segsh_dt);
   fChain->SetBranchAddress("segsh.tot0", segsh_tot0, &b_segsh_tot0);
   fChain->SetBranchAddress("segsh.tot1", segsh_tot1, &b_segsh_tot1);
   fChain->SetBranchAddress("segsh.edep", segsh_edep, &b_segsh_edep);
   //   fChain->SetBranchAddress("segsh.DaqStrawHit", segsh_DaqStrawHit, &b_segsh_DaqStrawHit);
   fChain->SetBranchAddress("segsh.rdrift", segsh_rdrift, &b_segsh_rdrift);
   fChain->SetBranchAddress("segsh.doca", segsh_doca, &b_segsh_doca);
   fChain->SetBranchAddress("segsh.drho", segsh_drho, &b_segsh_drho);
   fChain->SetBranchAddress("segsh.iseg", segsh_iseg, &b_segsh_iseg);
   fChain->SetBranchAddress("segsh.itrk", segsh_itrk, &b_segsh_itrk);
   Notify();
}

bool digis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

//-----------------------------------------------------------------------------
void digis::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t digis::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


//_____________________________________________________________________________
void     digis::AddHistogram(TObject* hist, const char* FolderName) {
  TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
  fol->Add(hist); 
}

//_____________________________________________________________________________
void digis::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void digis::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const float* LowEdge,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,LowEdge);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void digis::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // this is introduced specifically for weighting the DIO Mu2e spectrum
  // book 1D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void digis::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const double* LowEdge,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,LowEdge);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void digis::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 Int_t Ny, Double_t YMin, Double_t YMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void digis::HProf(TProfile*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       Double_t YMin, Double_t YMax,
		       const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
  AddHistogram(Hist,FolderName);
}




//_____________________________________________________________________________
void digis::DeleteHistograms(TFolder* Folder) {
  // internal method...

  if (((long int) Folder) == -1) Folder = fFolder;

  TObject  *o1;

  TIter    it1(Folder->GetListOfFolders());

  while ((o1 = it1.Next())) {
    if (o1->InheritsFrom("TFolder")) {
      DeleteHistograms((TFolder*) o1);
    }
    else if (o1->InheritsFrom("TH1")) {
      Folder->Remove(o1);
      delete o1;
    }
  }
}

//_____________________________________________________________________________
int  digis::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

  TDirectory*  dir;
  TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
  Dir->cd();
  //  dir = new TDirectory(Folder->GetName(),Folder->GetName(),"");
  dir = Dir->mkdir(Folder->GetName(),Folder->GetTitle());
  dir->cd();

//   printf(" ------------------- Dir: %s, new dir: %s\n",
// 	 Dir->GetName(),dir->GetName());


  TIter  it(Folder->GetListOfFolders());
  while ((o = it.Next())) {
//     printf(" o->GetName, o->ClassName : %-20s %-20s\n",
// 	   o->GetName(),
// 	   o->ClassName());

    if (strcmp(o->ClassName(),"TFolder") == 0) {
      SaveFolder((TFolder*) o, dir);
      //      dir->cd();
    }
    else if (! o->InheritsFrom("TStnModule")) {
      //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
      o->Write();
      //      gDirectory->GetListOfKeys()->Print();
    }
  }

  Dir->cd();
  return 0;
}

//_____________________________________________________________________________
void digis::SaveHist(const char* Filename) {
  // save histograms booked by all the modules into a file with the given name
  // Mode = 1: save folders

  TFile* f = new TFile(Filename,"recreate");

  SaveFolder(fFolder,f);

  f->Close();
  delete f;
}


