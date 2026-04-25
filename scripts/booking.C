///////////////////////////////////////////////////////////////////////////////
// helper class , helps to store histograms
///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TFolder.h"

class Booking {
public:
  TFolder* fFolder;                     // not owned, just cached

  Booking(TFolder* Folder) {
    fFolder = Folder;
  }

  void     AddHistogram(TObject* hist, TFolder* Folder) {
    if (Folder) {
      Folder->Add(hist);
    }
    else {
      std::cout << std::format("ERROR: folder not found\n");
    }
  }

//_____________________________________________________________________________
  void HBook1F(TH1F*& Hist, const char* Name, const char* Title,
               Int_t Nx, Double_t XMin, Double_t XMax,
               TFolder* Folder)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TH1F(Name,Title,Nx,XMin,XMax);
    AddHistogram(Hist,Folder);
  }

//_____________________________________________________________________________
  void HBook2F(TH2F*& Hist, const char* Name, const char* Title,
               Int_t Nx, Double_t XMin, Double_t XMax,
               Int_t Ny, Double_t YMin, Double_t YMax,
               TFolder* Folder)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user
    
    Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
    AddHistogram(Hist,Folder);
  }


//_____________________________________________________________________________
  void HProf(TProfile*& Hist, const char* Name, const char* Title,
                         Int_t Nx, Double_t XMin, Double_t XMax,
                         Double_t YMin, Double_t YMax,
                         TFolder* Folder)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
    AddHistogram(Hist,Folder);
  }

  //-----------------------------------------------------------------------------
  int  SaveFolder(TFolder* Folder, TDirectory* Dir) {
    // save Folder into a subdirectory
    // do not write TStnModule's - for each TStnModule save contents of its
    // fFolder

    //  TFolder*     fol;
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
      // else if (! o->InheritsFrom("TStnModule")) {
      else {
        //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
        o->Write();
        //      gDirectory->GetListOfKeys()->Print();
      }
    }

    Dir->cd();
    return 0;
  }

//_____________________________________________________________________________
  void SaveHist(const char* Filename) 
  {

    TFile* f = new TFile(Filename,"recreate");

    SaveFolder(fFolder,f);
    f->Close();
    delete f;
  }
  
};
