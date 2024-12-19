//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 17 16:35:19 2024 by ROOT version 6.30/04
// from TTree CalibTree/CalibTree
// found on file: EGamma0_40to60_v4.root
//////////////////////////////////////////////////////////

#ifndef IsoTrackV2_h
#define IsoTrackV2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class IsoTrackV2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           t_Run;
   Int_t           t_Event;
   Int_t           t_DataType;
   Int_t           t_ieta;
   Int_t           t_iphi;
   Double_t        t_EventWeight;
   Int_t           t_nVtx;
   Int_t           t_nTrk;
   Int_t           t_goodPV;
   Double_t        t_l1pt;
   Double_t        t_l1eta;
   Double_t        t_l1phi;
   Double_t        t_l3pt;
   Double_t        t_l3eta;
   Double_t        t_l3phi;
   Double_t        t_p;
   Double_t        t_pt;
   Double_t        t_phi;
   Double_t        t_mindR1;
   Double_t        t_mindR2;
   Double_t        t_eMipDR;
  //Double_t        t_eMipDR2;
  //Double_t        t_eMipDR3;
  //Double_t        t_eMipDR4;
  //Double_t        t_eMipDR5;
   Double_t        t_eHcal;
   Double_t        t_eHcal10;
   Double_t        t_eHcal30;
   Double_t        t_hmaxNearP;
  //Double_t        t_emaxNearP;
  //Double_t        t_eAnnular;
  //Double_t        t_hAnnular;
   Double_t        t_rhoh;
   Bool_t          t_selectTk;
   Bool_t          t_qltyFlag;
   Bool_t          t_qltyMissFlag;
   Bool_t          t_qltyPVFlag;
   Double_t        t_gentrackP;
   vector<unsigned int> *t_DetIds;
   vector<double>  *t_HitEnergies;
   vector<bool>    *t_trgbits;
   vector<unsigned int> *t_DetIds1;
   vector<unsigned int> *t_DetIds3;
   vector<double>  *t_HitEnergies1;
   vector<double>  *t_HitEnergies3;

   // List of branches
   TBranch        *b_t_Run;   //!
   TBranch        *b_t_Event;   //!
   TBranch        *b_t_DataType;   //!
   TBranch        *b_t_ieta;   //!
   TBranch        *b_t_iphi;   //!
   TBranch        *b_t_EventWeight;   //!
   TBranch        *b_t_nVtx;   //!
   TBranch        *b_t_nTrk;   //!
   TBranch        *b_t_goodPV;   //!
   TBranch        *b_t_l1pt;   //!
   TBranch        *b_t_l1eta;   //!
   TBranch        *b_t_l1phi;   //!
   TBranch        *b_t_l3pt;   //!
   TBranch        *b_t_l3eta;   //!
   TBranch        *b_t_l3phi;   //!
   TBranch        *b_t_p;   //!
   TBranch        *b_t_pt;   //!
   TBranch        *b_t_phi;   //!
   TBranch        *b_t_mindR1;   //!
   TBranch        *b_t_mindR2;   //!
   TBranch        *b_t_eMipDR;   //!
  //TBranch        *b_t_eMipDR2;   //!
  //TBranch        *b_t_eMipDR3;   //!
  //TBranch        *b_t_eMipDR4;   //!
  //TBranch        *b_t_eMipDR5;   //!
   TBranch        *b_t_eHcal;   //!
   TBranch        *b_t_eHcal10;   //!
   TBranch        *b_t_eHcal30;   //!
   TBranch        *b_t_hmaxNearP;   //!
  //TBranch        *b_t_emaxNearP;   //!
  //TBranch        *b_t_eAnnular;   //!
  //TBranch        *b_t_hAnnular;   //!
   TBranch        *b_t_rhoh;   //!
   TBranch        *b_t_selectTk;   //!
   TBranch        *b_t_qltyFlag;   //!
   TBranch        *b_t_qltyMissFlag;   //!
   TBranch        *b_t_qltyPVFlag;   //!
   TBranch        *b_t_gentrackP;   //!
   TBranch        *b_t_DetIds;   //!
   TBranch        *b_t_HitEnergies;   //!
   TBranch        *b_t_trgbits;   //!
   TBranch        *b_t_DetIds1;   //!
   TBranch        *b_t_DetIds3;   //!
   TBranch        *b_t_HitEnergies1;   //!
   TBranch        *b_t_HitEnergies3;   //!

   IsoTrackV2(TTree *tree=0);
   virtual ~IsoTrackV2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef IsoTrackV2_cxx
IsoTrackV2::IsoTrackV2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("EGamma0_40to60_v4.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("EGamma0_40to60_v4.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("EGamma0_40to60_v4.root:/hcalIsoTrackAnalyzer");
      dir->GetObject("CalibTree",tree);

   }
   Init(tree);
}

IsoTrackV2::~IsoTrackV2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t IsoTrackV2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t IsoTrackV2::LoadTree(Long64_t entry)
{
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

void IsoTrackV2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   t_DetIds = 0;
   t_HitEnergies = 0;
   t_trgbits = 0;
   t_DetIds1 = 0;
   t_DetIds3 = 0;
   t_HitEnergies1 = 0;
   t_HitEnergies3 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("t_Run", &t_Run, &b_t_Run);
   fChain->SetBranchAddress("t_Event", &t_Event, &b_t_Event);
   fChain->SetBranchAddress("t_DataType", &t_DataType, &b_t_DataType);
   fChain->SetBranchAddress("t_ieta", &t_ieta, &b_t_ieta);
   fChain->SetBranchAddress("t_iphi", &t_iphi, &b_t_iphi);
   fChain->SetBranchAddress("t_EventWeight", &t_EventWeight, &b_t_EventWeight);
   fChain->SetBranchAddress("t_nVtx", &t_nVtx, &b_t_nVtx);
   fChain->SetBranchAddress("t_nTrk", &t_nTrk, &b_t_nTrk);
   fChain->SetBranchAddress("t_goodPV", &t_goodPV, &b_t_goodPV);
   fChain->SetBranchAddress("t_l1pt", &t_l1pt, &b_t_l1pt);
   fChain->SetBranchAddress("t_l1eta", &t_l1eta, &b_t_l1eta);
   fChain->SetBranchAddress("t_l1phi", &t_l1phi, &b_t_l1phi);
   fChain->SetBranchAddress("t_l3pt", &t_l3pt, &b_t_l3pt);
   fChain->SetBranchAddress("t_l3eta", &t_l3eta, &b_t_l3eta);
   fChain->SetBranchAddress("t_l3phi", &t_l3phi, &b_t_l3phi);
   fChain->SetBranchAddress("t_p", &t_p, &b_t_p);
   fChain->SetBranchAddress("t_pt", &t_pt, &b_t_pt);
   fChain->SetBranchAddress("t_phi", &t_phi, &b_t_phi);
   fChain->SetBranchAddress("t_mindR1", &t_mindR1, &b_t_mindR1);
   fChain->SetBranchAddress("t_mindR2", &t_mindR2, &b_t_mindR2);
   fChain->SetBranchAddress("t_eMipDR", &t_eMipDR, &b_t_eMipDR);
   //fChain->SetBranchAddress("t_eMipDR2", &t_eMipDR2, &b_t_eMipDR2);
   //fChain->SetBranchAddress("t_eMipDR3", &t_eMipDR3, &b_t_eMipDR3);
   //fChain->SetBranchAddress("t_eMipDR4", &t_eMipDR4, &b_t_eMipDR4);
   //fChain->SetBranchAddress("t_eMipDR5", &t_eMipDR5, &b_t_eMipDR5);
   fChain->SetBranchAddress("t_eHcal", &t_eHcal, &b_t_eHcal);
   fChain->SetBranchAddress("t_eHcal10", &t_eHcal10, &b_t_eHcal10);
   fChain->SetBranchAddress("t_eHcal30", &t_eHcal30, &b_t_eHcal30);
   fChain->SetBranchAddress("t_hmaxNearP", &t_hmaxNearP, &b_t_hmaxNearP);
   //fChain->SetBranchAddress("t_emaxNearP", &t_emaxNearP, &b_t_emaxNearP);
   //fChain->SetBranchAddress("t_eAnnular", &t_eAnnular, &b_t_eAnnular);
   //fChain->SetBranchAddress("t_hAnnular", &t_hAnnular, &b_t_hAnnular);
   fChain->SetBranchAddress("t_rhoh", &t_rhoh, &b_t_rhoh);
   fChain->SetBranchAddress("t_selectTk", &t_selectTk, &b_t_selectTk);
   fChain->SetBranchAddress("t_qltyFlag", &t_qltyFlag, &b_t_qltyFlag);
   fChain->SetBranchAddress("t_qltyMissFlag", &t_qltyMissFlag, &b_t_qltyMissFlag);
   fChain->SetBranchAddress("t_qltyPVFlag", &t_qltyPVFlag, &b_t_qltyPVFlag);
   fChain->SetBranchAddress("t_gentrackP", &t_gentrackP, &b_t_gentrackP);
   fChain->SetBranchAddress("t_DetIds", &t_DetIds, &b_t_DetIds);
   fChain->SetBranchAddress("t_HitEnergies", &t_HitEnergies, &b_t_HitEnergies);
   fChain->SetBranchAddress("t_trgbits", &t_trgbits, &b_t_trgbits);
   fChain->SetBranchAddress("t_DetIds1", &t_DetIds1, &b_t_DetIds1);
   fChain->SetBranchAddress("t_DetIds3", &t_DetIds3, &b_t_DetIds3);
   fChain->SetBranchAddress("t_HitEnergies1", &t_HitEnergies1, &b_t_HitEnergies1);
   fChain->SetBranchAddress("t_HitEnergies3", &t_HitEnergies3, &b_t_HitEnergies3);
   Notify();
}

Bool_t IsoTrackV2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void IsoTrackV2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t IsoTrackV2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef IsoTrackV2_cxx
