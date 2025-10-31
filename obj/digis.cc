#define digis_cxx
#include "daqana/obj/digis.hh"
#include <TH2.h>
// #include <TStyle.h>
// #include <TCanvas.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <format>

//-----------------------------------------------------------------------------
void digis::Loop(long int NEvents) {
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

   long int nevents = NEvents;
   if (nevents < 0) nevents = nentries;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nevents; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     CalculateMissingParameters();
     FillHistograms(&fHist);
   }
}




//-----------------------------------------------------------------------------
int digis::CalculateMissingParameters() {

  fNSeg4 = 0;
  fNSeg6 = 0;
  fNSeg8 = 0;
                                        // seelct "good" segments
  for (int i=0; i<nseg; i++) {
    if ((seg_ngh[i] >= fMinSegNHits[0]) and (seg_chi2d[i] < fMaxSegChi2d) and
        (seg_nghl[i][0] >= fMinSegNghl) and  (seg_nghl[i][1] >= fMinSegNghl)) {
      fISeg4[fNSeg4] = i;
      fNSeg4        += 1;
    }
    if ((seg_ngh[i] >= fMinSegNHits[1]) and (seg_chi2d[i] < fMaxSegChi2d) and
        (seg_nghl[i][0] >= fMinSegNghl) and  (seg_nghl[i][1] >= fMinSegNghl)) {
      fISeg6[fNSeg6] = i;
      fNSeg6        += 1;
    }
    if ((seg_ngh[i] >= fMinSegNHits[2]) and (seg_chi2d[i] < fMaxSegChi2d) and
        (seg_nghl[i][0] >= fMinSegNghl) and  (seg_nghl[i][1] >= fMinSegNghl)) {
      fISeg8[fNSeg8] = i;
      fNSeg8        += 1;
    }
  }

  return 0;
};

//-----------------------------------------------------------------------------
int digis::FillTwoSegHistograms(TwoSegHist_t* Hist, TwoSegPar_t* Par) {

  //  printf("Par->hid:%i Par->dt: %10.3f hist:0x%8p\n",Par->hid,Par->dt,(void*) Hist->fDt[Par->hid]);

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
int digis::FillPanelShHistograms(PanelShHist_t* Hist, int I) {
  Hist->fRdrift->Fill(segsh_rdrift[I]);
  Hist->fDrhoVsRdrift->Fill(segsh_rdrift[I],segsh_drho[I]);

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
//-----------------------------------------------------------------------------
// fMinSegNHits[0]: 4-hit segments
//-----------------------------------------------------------------------------
    if (seg_nh[i] >= fMinSegNHits[0]) {
      FillSegmentHistograms(Hist->fSegment[401],i);
      FillSegmentHistograms(Hist->fSegment[unique_panel+410],i);
      if ((seg_nghl[i][0] > 0) and (seg_nghl[i][1] > 0)) {
        FillSegmentHistograms(Hist->fSegment[402],i);
        if (seg_chi2d[i] <fMaxSegChi2d) {
          FillSegmentHistograms(Hist->fSegment[403],i);
          if (fNSeg6 >= 2) {
            FillSegmentHistograms(Hist->fSegment[404],i);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// fMinSegNHits[1]: 6-hit segments
//-----------------------------------------------------------------------------
    if (seg_nh[i] >= fMinSegNHits[1]) {
      FillSegmentHistograms(Hist->fSegment[601],i);
      FillSegmentHistograms(Hist->fSegment[unique_panel+610],i);
      if ((seg_nghl[i][0] > 0) and (seg_nghl[i][1] > 0)) {
        FillSegmentHistograms(Hist->fSegment[602],i);
        if (seg_chi2d[i] <fMaxSegChi2d) {
          FillSegmentHistograms(Hist->fSegment[603],i);
          if (fNSeg6 >= 2) {
            FillSegmentHistograms(Hist->fSegment[604],i);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// fMinSegNHits[2]: 8-hit segments
//-----------------------------------------------------------------------------
    if (seg_nh[i] >= fMinSegNHits[2]) {
      FillSegmentHistograms(Hist->fSegment[801],i);
      FillSegmentHistograms(Hist->fSegment[unique_panel+810],i);
      if ((seg_nghl[i][0] > 0) and (seg_nghl[i][1] > 0)) {
        FillSegmentHistograms(Hist->fSegment[802],i);
        if (seg_chi2d[i] <fMaxSegChi2d) {
          FillSegmentHistograms(Hist->fSegment[803],i);
          if (fNSeg8 >= 2) {
            FillSegmentHistograms(Hist->fSegment[804],i);
          }
        }
      }
    }
  }
//-----------------------------------------------------------------------------  
// dt between two 4+ chi2d<10 segments
//-----------------------------------------------------------------------------
  TwoSegPar_t par_ss;

  for (int i1=0; i1<fNSeg4-1; i1++) {
    int iseg1 = fISeg4[i1];
    int plane1 = plane(seg_sid[iseg1]);
    int panel1 = panel(seg_sid[iseg1]);
    int ip1    = 6*plane1+panel1;
    for (int i2=i1+1; i2<fNSeg4; i2++) {
      int iseg2  = fISeg4[i2];
      int plane2 = plane(seg_sid[iseg2]);
      int panel2 = panel(seg_sid[iseg2]);
      int ip2    = 6*plane2+panel2;
//-----------------------------------------------------------------------------
// two "good" 4-hit segments, figure the histogram id
//-----------------------------------------------------------------------------
      // int   hid;
      // float dt;
      if (ip1 < ip2) {
        par_ss.hid = 100*ip1+ip2;     // 4-hit ones start from hist_id=40
        par_ss.dt  = seg_t0[iseg1]-seg_t0[iseg2];
      }
      else {
        par_ss.hid = 100*ip2+ip1;
        par_ss.dt  = seg_t0[iseg2]-seg_t0[iseg1];
      }
      if (plane1 == plane2) par_ss.diff_planes = 0;
      else                  par_ss.diff_planes = 1;

      FillTwoSegHistograms(Hist->fTwoSeg[0],&par_ss);
      if (seg_dzdy[i1] > 0) FillTwoSegHistograms(Hist->fTwoSeg[1],&par_ss);
      else                  FillTwoSegHistograms(Hist->fTwoSeg[2],&par_ss);
    }
  }
//-----------------------------------------------------------------------------  
// dt between two 6+ chi2d<10 segments
//-----------------------------------------------------------------------------
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
// two "good" 6-hit segments, figure the histogram id
//-----------------------------------------------------------------------------
      // int   hid;
      // float dt;
      if (ip1 < ip2) {
        par_ss.hid = 100*ip1+ip2;     // 4-hit ones start from hist_id=40
        par_ss.dt  = seg_t0[iseg1]-seg_t0[iseg2];
      }
      else {
        par_ss.hid = 100*ip2+ip1;
        par_ss.dt  = seg_t0[iseg2]-seg_t0[iseg1];
      }
      if (plane1 == plane2) par_ss.diff_planes = 0;
      else                  par_ss.diff_planes = 1;

      FillTwoSegHistograms(Hist->fTwoSeg[10],&par_ss);
      if (seg_dzdy[i1] > 0) FillTwoSegHistograms(Hist->fTwoSeg[11],&par_ss);
      else                  FillTwoSegHistograms(Hist->fTwoSeg[12],&par_ss);
    }
  }
//-----------------------------------------------------------------------------
// 4-hit segment straw hit histograms (SSHT) ... do that only for good segments
// and use only hits with 0.5 < R < 2.0 mm
//-----------------------------------------------------------------------------
  for (int is=0; is<fNSeg4; is++) {
    int iseg = fISeg4[is];
    for (int ish=0; ish<nsegsh; ish++) {
      if (segsh_iseg[ish] != iseg)         continue;
      int plane = (segsh_sid[ish] >> 10) & 0x3f;
      int panel = (segsh_sid[ish] >>  7) & 0x7;
//-----------------------------------------------------------------------------
// each panel separately
//-----------------------------------------------------------------------------
      float rdrift = segsh_rdrift[ish];
      if ((rdrift > 0.5) and (rdrift < 2.0)) {
        int iset = 400+plane*6+panel;
        FillSshtHistograms(Hist->fSsht[iset], ish);
      }
    }
  }
//-----------------------------------------------------------------------------
// 6-hit segment straw hit histograms (SSHT) ... do that only for good segments
// and use only hits with 0.5 < R < 2.0 mm
//-----------------------------------------------------------------------------
  for (int is=0; is<fNSeg6; is++) {
    int iseg = fISeg6[is];
    for (int ish=0; ish<nsegsh; ish++) {
      if (segsh_iseg[ish] != iseg)         continue;
      int plane = (segsh_sid[ish] >> 10) & 0x3f;
      int panel = (segsh_sid[ish] >>  7) & 0x7;
//-----------------------------------------------------------------------------
// each panel separately
//-----------------------------------------------------------------------------
      float rdrift = segsh_rdrift[ish];
      if ((rdrift > 0.5) and (rdrift < 2.0)) {
        int iset = 600+plane*6+panel;
        FillSshtHistograms(Hist->fSsht[iset], ish);
      }
    }
  }
//-----------------------------------------------------------------------------
// 8-hit segment straw hit histograms (SSHT) ... do that only for good segments
// and use only hits with 0.5 < R < 2.0 mm
//-----------------------------------------------------------------------------
  for (int is=0; is<fNSeg8; is++) {
    int iseg = fISeg8[is];
    for (int ish=0; ish<nsegsh; ish++) {
      if (segsh_iseg[ish] != iseg)         continue;
      int plane = (segsh_sid[ish] >> 10) & 0x3f;
      int panel = (segsh_sid[ish] >>  7) & 0x7;
      int unique_panel = 6*plane+panel;
       FillPanelShHistograms(Hist->fPanelSh[unique_panel],ish);
//-----------------------------------------------------------------------------
// each panel separately
//-----------------------------------------------------------------------------
      float rdrift = segsh_rdrift[ish];
      if ((rdrift > 0.5) and (rdrift < 2.0)) {
        int iset = 800+plane*6+panel;
        FillSshtHistograms(Hist->fSsht[iset], ish);
      }
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
// add more descriptive titles later
//-----------------------------------------------------------------------------
int digis::BookPanelShHistograms(PanelShHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fRdrift      ,"r"        ,"rdrift, mm"   , 100,    0, 5, Folder);
  HBook2F(Hist->fDrhoVsRdrift,"drho_vs_r","drho_vs_r, mm",  100,   0, 5, 200,-1,1,Folder);
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

  book_segment_histset[  0] = 1;		// all segments

  book_segment_histset[401] = 1;	        // segments with 4+ good hits
  book_segment_histset[402] = 1;	        // segments with 4+ good hits and chi2d<10
  book_segment_histset[403] = 1;	        // segments with 4+ good hits and chi2d<10 and nghl[i] > 0
  book_segment_histset[404] = 1;	        // events with 2 such segments

  book_segment_histset[410] = 1;	        // events with segments in panel 0 
  book_segment_histset[411] = 1;	        // events with segments in panel 1
  book_segment_histset[412] = 1;	        // events with segments in panel 2
  book_segment_histset[413] = 1;	        // events with segments in panel 3
  book_segment_histset[414] = 1;	        // events with segments in panel 4
  book_segment_histset[415] = 1;	        // events with segments in panel 5
  book_segment_histset[416] = 1;	        // events with segments in panel 6
  book_segment_histset[417] = 1;	        // events with segments in panel 7
  book_segment_histset[418] = 1;	        // events with segments in panel 8
  book_segment_histset[419] = 1;	        // events with segments in panel 9
  book_segment_histset[420] = 1;	        // events with segments in panel 10
  book_segment_histset[421] = 1;	        // events with segments in panel 11

  book_segment_histset[601] = 1;	        // segments with 6+ good hits
  book_segment_histset[602] = 1;	        // segments with 6+ good hits and chi2d<10
  book_segment_histset[603] = 1;	        // segments with 6+ good hits and chi2d<10 and nghl[i] > 0
  book_segment_histset[604] = 1;	        // events with 2 such segments

  book_segment_histset[610] = 1;	        // events with segments in panel 0 
  book_segment_histset[611] = 1;	        // events with segments in panel 1
  book_segment_histset[612] = 1;	        // events with segments in panel 2
  book_segment_histset[613] = 1;	        // events with segments in panel 3
  book_segment_histset[614] = 1;	        // events with segments in panel 4
  book_segment_histset[615] = 1;	        // events with segments in panel 5
  book_segment_histset[616] = 1;	        // events with segments in panel 6
  book_segment_histset[617] = 1;	        // events with segments in panel 7
  book_segment_histset[618] = 1;	        // events with segments in panel 8
  book_segment_histset[619] = 1;	        // events with segments in panel 9
  book_segment_histset[620] = 1;	        // events with segments in panel 10
  book_segment_histset[621] = 1;	        // events with segments in panel 11

  book_segment_histset[801] = 1;	        // segments with 6+ good hits
  book_segment_histset[802] = 1;	        // segments with 6+ good hits and chi2d<10
  book_segment_histset[803] = 1;	        // segments with 6+ good hits and chi2d<10 and nghl[i] > 0
  book_segment_histset[804] = 1;	        // events with 2 such segments

  book_segment_histset[810] = 1;	        // events with good 8+ hit segments in panel 0 
  book_segment_histset[811] = 1;	        // events with good 8+ hit segments in panel 1
  book_segment_histset[812] = 1;	        // events with good 8+ hit segments in panel 2
  book_segment_histset[813] = 1;	        // events with good 8+ hit segments in panel 3
  book_segment_histset[814] = 1;	        // events with good 8+ hit segments in panel 4
  book_segment_histset[815] = 1;	        // events with good 8+ hit segments in panel 5
  book_segment_histset[816] = 1;	        // events with good 8+ hit segments in panel 6
  book_segment_histset[817] = 1;	        // events with good 8+ hit segments in panel 7
  book_segment_histset[818] = 1;	        // events with good 8+ hit segments in panel 8
  book_segment_histset[819] = 1;	        // events with good 8+ hit segments in panel 9
  book_segment_histset[820] = 1;	        // events with good 8+ hit segments in panel 10
  book_segment_histset[821] = 1;	        // events with good 8+ hit segments in panel 11

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

  book_twoseg_histset[ 0] = 1;		// all twosegs nhits >= 4
  book_twoseg_histset[ 1] = 1;		// all twosegs nhits >= 4 dzdy(p1) > 0
  book_twoseg_histset[ 2] = 1;		// all twosegs nhits >= 4 dzdy(p1) < 0
  book_twoseg_histset[10] = 1;		// all twosegs nhits >= 6
  book_twoseg_histset[11] = 1;		// all twosegs nhits >= 6 dzdy(p1) > 0
  book_twoseg_histset[12] = 1;		// all twosegs nhits >= 6 dzdy(p1) < 0

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

  book_ssht_histset[400] = 1;	        // good segments with 4+ hits panel#0  (+10)
  book_ssht_histset[401] = 1;	        // good segments with 4+ hits panel#1
  book_ssht_histset[402] = 1;	        // good segments with 4+ hits panel#2
  book_ssht_histset[403] = 1;	        // good segments with 4+ hits panel#3
  book_ssht_histset[404] = 1;	        // good segments with 4+ hits panel#4
  book_ssht_histset[405] = 1;	        // good segments with 4+ hits panel#5
  book_ssht_histset[406] = 1;	        // good segments with 4+ hits panel#6
  book_ssht_histset[407] = 1;	        // good segments with 4+ hits panel#7
  book_ssht_histset[408] = 1;	        // good segments with 4+ hits panel#8
  book_ssht_histset[409] = 1;	        // good segments with 4+ hits panel#9
  book_ssht_histset[410] = 1;	        // good segments with 4+ hits panel#10
  book_ssht_histset[411] = 1;	        // good segments with 4+ hits panel#11

  book_ssht_histset[600] = 1;	        // good segments with 6+ hits panel#0  (+10)
  book_ssht_histset[601] = 1;	        // good segments with 6+ hits panel#1
  book_ssht_histset[602] = 1;	        // good segments with 6+ hits panel#2
  book_ssht_histset[603] = 1;	        // good segments with 6+ hits panel#3
  book_ssht_histset[604] = 1;	        // good segments with 6+ hits panel#4
  book_ssht_histset[605] = 1;	        // good segments with 6+ hits panel#5
  book_ssht_histset[606] = 1;	        // good segments with 6+ hits panel#6
  book_ssht_histset[607] = 1;	        // good segments with 6+ hits panel#7
  book_ssht_histset[608] = 1;	        // good segments with 6+ hits panel#8
  book_ssht_histset[609] = 1;	        // good segments with 6+ hits panel#9
  book_ssht_histset[610] = 1;	        // good segments with 6+ hits panel#10
  book_ssht_histset[611] = 1;	        // good segments with 6+ hits panel#11

  book_ssht_histset[800] = 1;	        // good segments with 8+ hits panel#0  (+10)
  book_ssht_histset[801] = 1;	        // good segments with 8+ hits panel#1
  book_ssht_histset[802] = 1;	        // good segments with 8+ hits panel#2
  book_ssht_histset[803] = 1;	        // good segments with 8+ hits panel#3
  book_ssht_histset[804] = 1;	        // good segments with 8+ hits panel#4
  book_ssht_histset[805] = 1;	        // good segments with 8+ hits panel#5
  book_ssht_histset[806] = 1;	        // good segments with 8+ hits panel#6
  book_ssht_histset[807] = 1;	        // good segments with 8+ hits panel#7
  book_ssht_histset[808] = 1;	        // good segments with 8+ hits panel#8
  book_ssht_histset[809] = 1;	        // good segments with 8+ hits panel#9
  book_ssht_histset[810] = 1;	        // good segments with 8+ hits panel#10
  book_ssht_histset[811] = 1;	        // good segments with 8+ hits panel#11

  for (int i=0; i<kNSshtHistSets; i++) {
    if (book_ssht_histset[i] != 0) {
      sprintf(folder_name,"ssht_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fSsht[i] = new SshtHist_t;
      BookSshtHistograms(fHist.fSsht[i],Form("%s",folder_name));
    }
  }

//-----------------------------------------------------------------------------
// book panel straw hit histograms
//-----------------------------------------------------------------------------
  int book_panel_sh_histset[kNPanelShHistSets];
  for (int i=0; i<kNPanelShHistSets; i++) book_panel_sh_histset[i] = 0;

  book_panel_sh_histset[ 0] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 1] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 2] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 3] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 4] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 5] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 6] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 7] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 8] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[ 9] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[10] = 1;		// panel hist for segments with nhits >= 8
  book_panel_sh_histset[11] = 1;		// panel hist for segments with nhits >= 8

  for (int i=0; i<kNTwoSegHistSets; i++) {
    if (book_panel_sh_histset[i] != 0) {
      sprintf(folder_name,"psh_%i",i);
      fol = (TFolder*) fFolder->FindObject(folder_name);
      if (! fol) fol = fFolder->AddFolder(folder_name,folder_name);
      fHist.fPanelSh[i] = new PanelShHist_t;
      BookPanelShHistograms(fHist.fPanelSh[i],Form("%s",folder_name));
    }
  }

  return 0;
};


//-----------------------------------------------------------------------------
// if Dsid < 0, a single file
// otherwise, a specific dataset, ignore Fn
//-----------------------------------------------------------------------------
digis::digis(const std::string& Fn, const std::string& Fileset) : fChain(0) {
  struct NestedFunctor {
    std::string operator()(std::string& s) {
      std::string x = s;
      size_t start = x.find_first_not_of(" \t\r\n");
      size_t end   = x.find_last_not_of (" \t\r\n");
      if (start == std::string::npos) return "";
      x = x.substr(start, end - start + 1);
      return x;
    }
  };

  fMinSegNHits[0] =  4;
  fMinSegNHits[1] =  6;
  fMinSegNHits[2] =  8;
  
  fMinSegNghl     =  2;
  fMaxSegChi2d    = 10.;

  fFolder  = gROOT->GetRootFolder()->AddFolder("digis","digis");

  fChain    = new TChain("/MakeDigiNtuple/digis");

  if (Fileset == "file") {
    fChain->AddFile(Fn.data(),TChain::kBigNumber);
  }
  else if (Fn == "vst00s0s10r0000") {
    fChain->AddFile("./nts/nts.murat.vst00s0s10r0000.daqana.107995_000001.root",TChain::kBigNumber);
    fChain->AddFile("./nts/nts.murat.vst00s0s10r0000.daqana.107995_000148.root",TChain::kBigNumber);
    fChain->AddFile("./nts/nts.murat.vst00s0s10r0000.daqana.107995_000288.root",TChain::kBigNumber);
    fChain->AddFile("./nts/nts.murat.vst00s0s10r0000.daqana.107995_000429.root",TChain::kBigNumber);
  }
  else if (Fn == "vst04s0s10r0000") {
                                        // read text file from a catalog
    std::string fn = "daqana/datasets/vst04s0/catalog/nts.murat."+Fn+".daqana.root.files";
    if (Fileset != "") fn = fn + "." + Fileset;

    // std::cout << "fn:" << fn << std::endl;
    std::ifstream input(fn);
    std::string line;
    NestedFunctor nested_trim;
    while (std::getline(input, line)) {
      std::string trimmed_line = nested_trim(line);
      //      std::cout << std::format("line:{} trimmed_line:{}\n",line,trimmed_line);
      if (trimmed_line.empty() or (trimmed_line[0] == '#'))                   continue;
//------------------------------------------------------------------------------
// expect the line to contain the filename
//-----------------------------------------------------------------------------
      fChain->AddFile(trimmed_line.data(),TChain::kBigNumber);
    }
    input.close();
  }

  fCurrent = -1;
  fChain->SetMakeClass(1);

  Init();

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
void digis::Init() {
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
      digis::SaveFolder((TFolder*) o, dir);
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

  digis::SaveFolder(fFolder,f);

  f->Close();
  delete f;
}


