#define IsoTrack_cxx
#include "IsoTrack.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <iostream>

void unpackDetId(unsigned int detId, int& subdet, int& zside, int& ieta, int& iphi, int& depth) {
  // The maskings are defined in DataFormats/DetId/interface/DetId.h
  //                      and in DataFormats/HcalDetId/interface/HcalDetId.h
  // The macro does not invoke the classes there and use them
  subdet = ((detId >> 25) & (0x7));
  if ((detId & 0x1000000) == 0) {
    ieta = ((detId >> 7) & 0x3F);
    zside = (detId & 0x2000) ? (1) : (-1);
    depth = ((detId >> 14) & 0x1F);
    iphi = (detId & 0x3F);
  } else {
    ieta = ((detId >> 10) & 0x1FF);
    zside = (detId & 0x80000) ? (1) : (-1);
    depth = ((detId >> 20) & 0xF);
    iphi = (detId & 0x3FF);
  }
}

// From Long Wang, Re: ROOT file with HCalRespCorrs for comparison, 2024/2/1
double etaVal(int ieta) {
  //if (ieta>=0) ++ieta; // patch iEta!=0 (Mikko's addition)
  double etavl(0.);
  if (ieta <= -24)
    etavl = .1695*ieta + 1.9931;
  else if (ieta <= -1)
    etavl = .0875*ieta + .0489;
  else if  (ieta < 24)
    etavl = .0875*ieta - .0489;
  else
    etavl = .1695*ieta - 1.9931;
  return etavl;
} // etaVal


// Estimate effective area based on half granularity on iphi after |ieta|>=21
double areaScale3x5(int ieta) {

  double areaScale1 = 3.*5/(3.*9.-3.*5.); // = 1.25, |ieta|<21
  double areaScale2 = 3.*3./(3.*5-3.*3.); // = 1.50, |ieta|>=21
  if (abs(ieta)<20)  return areaScale1;
  if (abs(ieta)==20) return (2./3.*areaScale1+1./3.*areaScale2);
  if (abs(ieta)==21) return (1./3.*areaScale1+2./3.*areaScale2);
  if (abs(ieta)>21)  return areaScale2;
  return 1.5;
}

// Merge depths 1+2
bool doPerDepth = true;
bool enforce5x5 = false;
bool enforce3x3 = false;
bool enforce3x5 = true;
bool updateSingleDepth = true;
bool checkConsistency = true;
bool doSingleDepth = true;
bool mergeDepths1and2 = false;//true;

void IsoTrack::Loop()
{
//   In a ROOT session, you can do:
//      root> .L IsoTrack.C
//      root> IsoTrack t
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

   cout << "Calling IsoTrack::Loop()..." << endl << flush;
   
   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("t_p",1);
   fChain->SetBranchStatus("t_eHcal",1);
   fChain->SetBranchStatus("t_ieta",1);
   // MIP correction
   fChain->SetBranchStatus("t_eMipDR",1); // also track ID
   // PU
   fChain->SetBranchStatus("t_eHcal10",1);
   fChain->SetBranchStatus("t_eHcal30",1);
   // track ID
   fChain->SetBranchStatus("t_selectTk",1);
   fChain->SetBranchStatus("t_qltyMissFlag",1);
   fChain->SetBranchStatus("t_hmaxNearP",1);
   fChain->SetBranchStatus("t_mindR1",1);

   // Process depths (heavy!)
   if (doPerDepth) {
     // Sanity checks
     fChain->SetBranchStatus("t_iphi",1);
     // Core energy per depth
     fChain->SetBranchStatus("t_DetIds",1);
     fChain->SetBranchStatus("t_HitEnergies",1);
     // PU per depth (just first cone now)
     //fChain->SetBranchStatus("t_DetIds1",1);
     //fChain->SetBranchStatus("t_HitEnergies1",1);
     fChain->SetBranchStatus("t_DetIds3",1);
     fChain->SetBranchStatus("t_HitEnergies3",1);
   }
   //fChain->SetBranchStatus("t_DetIds",1);
   
   Long64_t nentries = fChain->GetEntriesFast();

   cout << "Opening output file IsoTrack.root" << endl << flush;
   cout << "Starting loop";
   TFile *fout = new TFile("IsoTrack.root","RECREATE");

   // Depth-independent results
   TH2D *h2raw, *h2pu1, *h2pu3, *h2mip, *h2c;
   h2raw = new TH2D("h2raw",";ieta;t_eHcal/(t_p-eMIP)", 61,-30.5,30.5, 400,0,4);
   h2pu1 = new TH2D("h2pu1",";ieta;(t_eHcal10-t_eHcal)/(t_p-eMIP)",
		    61,-30.5,30.5, 400,0,4);
   h2pu3 = new TH2D("h2pu3",";ieta;(t_eHcal30-t_eHcal10)/(t_p-eMIP)",
		    61,-30.5,30.5, 400,0,4);
   h2mip = new TH2D("h2mip",";ieta;t_eMipDR/t_p", 61,-30.5,30.5, 400,0,4);
   h2c = new TH2D("h2c",";ieta;(t_eHcal-ePU)/(t_p-eMIP)",
		  61,-30.5,30.5, 400,0,4);

   // Means
   TProfile *praw, *ppu1, *ppu3, *pmip, *pc, *pc_even, *pc_odd;
   praw = new TProfile("praw",";ieta;t_eHcal/(t_p-eMIP) [0.15,1.85]",
		       61,-30.5,30.5);
   ppu1 = new TProfile("ppu1",";ieta;(t_eHcal10-t_eHcal)/(t_p-eMIP)"
		       " [0.15,1.85]", 61,-30.5,30.5);
   ppu3 = new TProfile("ppu3",";ieta;(t_eHcal30-t_eHcal10)/(t_p-eMIP)"
		       " [0.15,1.85]", 61,-30.5,30.5);
   pmip = new TProfile("pmip",";ieta;t_eMipDR/(t_p-eMIP)"
		       " [0.15,1.85]", 61,-30.5,30.5);   
   pc = new TProfile("pc",";ieta;(t_eHcal-ePU)/(t_p-eMIP) [0.15,1.85]",
		     61,-30.5,30.5);
   pc_even = new TProfile("pc_even",";ieta;(t_eHcal-ePU)/(t_p-eMIP)"
			  " [0.15,1.85]", 61,-30.5,30.5);
   pc_odd = new TProfile("pc_odd",";ieta;(t_eHcal-ePU)/(t_p-eMIP)"
			 " [0.15,1.85]", 61,-30.5,30.5);

   // Same exercise now per depth
   // Depth=0 for ECAL (MIP or early shower)
   TH3D *h3raw, *h3pu1, *h3c;
   h3raw = new TH3D("h3raw",";ieta;depth; t_eHcal/(t_p-eMIP)",
		    61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);
   h3pu1 = new TH3D("h3pu1",";ieta;depth;(t_eHcal10-t_eHcal)/(t_p-eMIP)",
		    61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);
   h3c = new TH3D("h3c",";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)",
		  61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);

   // Means
   TProfile2D *p2raw, *p2pu1, *p2c, *p2c_even, *p2c_odd;
   p2raw = new TProfile2D("p2raw",";ieta;depth;t_eHcal/(t_p-eMIP) [0.15,1.85]",
			  61,-30.5,30.5, 10,-0.5,9.5);
   p2pu1 = new TProfile2D("p2pu1",";ieta;depth;(t_eHcal10-t_eHcal)/(t_p-eMIP)"
			  " [0.15,1.85]", 61,-30.5,30.5, 10,-0.5,9.5);
   p2c = new TProfile2D("p2c",";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)"
			" [0.15,1.85]", 61,-30.5,30.5, 10,-0.5,9.5);
   p2c_even = new TProfile2D("p2c_even",";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)"
			     " [0.15,1.85]", 61,-30.5,30.5, 10,-0.5,9.5);
   p2c_odd = new TProfile2D("p2c_odd",";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)"
			    " [0.15,1.85]", 61,-30.5,30.5, 10,-0.5,9.5);

   // Depth Covariances
   TProfile3D *p3c, *p3c_even, *p3c_odd;
   p3c = new TProfile3D("p3c",";ieta;depth;depth;(t_eHcal-ePU)/(t_p-eMIP) "
			"[0.15,1.85]", 61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);
   p3c_even = new TProfile3D("p3c_even",";ieta;depth;depth;"
			     "(t_eHcal-ePU)/(t_p-eMIP) [0.15,1.85]",
			     61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);
   p3c_odd = new TProfile3D("p3c_odd",";ieta;depth;depth;"
			    "(t_eHcal-ePU)/(t_p-eMIP) [0.15,1.85]",
			    61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);

   // Checking shape and energy distribution of isolation cone
   TH2D *h2draw_hb, *h2dpu1_hb;
   TH2D *h2draw_he1, *h2dpu1_he1, *h2draw_he2, *h2dpu1_he2;
   h2draw_hb = new TH2D("h2draw_hb",";#Delta|i#eta|;#Delta|i#phi|;",
			21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_hb = new TH2D("h2dpu1_hb",";#Delta|i#eta|;#Delta|i#phi|;",
			21,-10.5,10.5, 21,-10.5,10.5);
   h2draw_he1 = new TH2D("h2draw_he1",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_he1 = new TH2D("h2dpu1_he1",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2draw_he2 = new TH2D("h2draw_he2",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_he2 = new TH2D("h2dpu1_he2",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   
   TProfile2D *p2draw_hb, *p2dpu1_hb;
   TProfile2D *p2draw_he1, *p2dpu1_he1, *p2draw_he2, *p2dpu1_he2;
   p2draw_hb= new TProfile2D("p2draw_hb",";#Deltai#eta;#Deltai#phi;"
			     "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_hb = new TProfile2D("p2dpu1_hb",";#Deltai#eta;#Deltai#phi;"
			      "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2draw_he1 = new TProfile2D("p2draw_he1",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_he1 = new TProfile2D("p2dpu1_he1",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2draw_he2 = new TProfile2D("p2draw_he2",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_he2 = new TProfile2D("p2dpu1_he2",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%10000==0) cout << "." << flush;
      
      bool ok = ((t_selectTk) && (t_qltyMissFlag) && (t_hmaxNearP < 10.0) && (t_eMipDR < 1.0) && (t_mindR1 > 1.0));
      if (!ok) continue;
      
      //double eta = etaVal(t_ieta);
      if (doPerDepth) {
	double esum(0), esum1(0);
	double e[10] =  {0,0,0,0,0, 0,0,0,0,0};
	double e1[10] = {0,0,0,0,0, 0,0,0,0,0};
	e[0] = t_eMipDR;
	int subdet, zside, ieta, iphi, depth;

	// Core energies
	for (unsigned int idet = 0; idet != t_DetIds->size(); ++idet) {
	  unpackDetId((*t_DetIds)[idet], subdet, zside, ieta, iphi, depth);
	  ieta *= zside;
	  double edet = (*t_HitEnergies)[idet];

	  /*
	  if (fabs(ieta-t_ieta)>3) { // looking into 5x5? sometimes off-center 
	    cout << "jentry = " << jentry
		 << ", ieta = " << ieta << ", iphi = " << iphi
		 << ", t_ieta = " << t_ieta << ", edet = " << edet << endl;
	    continue;
	  }
	  if (depth<0 || depth>9) {
	    cout << "jentry = " << jentry
		 << ", ieta = " << ieta << ", iphi = " << iphi
		 << ", depth = " << depth << ", edet = " << edet << endl;
	    continue;
	  }
	  */
	  double dieta = (abs(ieta)-abs(t_ieta));
	  double diphi = (abs(iphi)-abs(t_iphi));
	  if (ieta*t_ieta<0) dieta += -1; // ieta=0 does not exist
	  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic

	  bool is5x5 = (fabs(dieta)<3 && fabs(diphi)<3);
	  bool is3x3 = (fabs(dieta)<2 && fabs(diphi)<2);
	  bool is3x5 = (fabs(dieta)<2 && fabs(diphi)<3);
	  bool noForce = (!enforce5x5 && !enforce3x3 && !enforce3x5);
	  
	  if ((is5x5 && enforce5x5) || (is3x3 && enforce3x3) ||
	      (is3x5 && enforce3x5) || noForce) {
	    
	    esum += edet;
	    e[depth] += edet;

	    if (abs(ieta)<=16) {
	      h2draw_hb->Fill(dieta, diphi);
	      p2draw_hb->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else if (abs(ieta)<=23) {
	      h2draw_he1->Fill(dieta, diphi);
	      p2draw_he1->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else {
	      h2draw_he2->Fill(dieta, diphi);
	      p2draw_he2->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	  } // is5x5
	} // for idet

	// Cone energies for PU estimates
	//for (unsigned int idet = 0; idet != t_DetIds1->size(); ++idet) {
	//unpackDetId((*t_DetIds1)[idet], subdet, zside, ieta, iphi, depth);
	//ieta *= zside;
	//double edet = (*t_HitEnergies1)[idet];

	for (unsigned int idet = 0; idet != t_DetIds3->size(); ++idet) {
	  unpackDetId((*t_DetIds3)[idet], subdet, zside, ieta, iphi, depth);
	  ieta *= zside;
	  double edet = (*t_HitEnergies3)[idet];

	  /*
	  if (fabs(ieta-t_ieta)>4) { // looking into 7x7?
	    cout << "jentry = " << jentry
		 << ", ieta1 = " << ieta << ", iphi = " << iphi
		 << ", t_ieta = " << t_ieta << ", edet = " << edet << endl;
	    continue;
	  }
	  if (depth<0 || depth>9) {
	    cout << "jentry = " << jentry
		 << ", ieta1 = " << ieta << ", iphi = " << iphi
		 << ", depth = " << depth << ", edet = " << edet << endl;
	    continue;
	  }
	  */

	  double dieta = (abs(ieta)-abs(t_ieta));
	  double diphi = (abs(iphi)-abs(t_iphi));
	  if (ieta*t_ieta<0) dieta += -1; // ieta=0 does not exist
	  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic
	  
	  //bool is5x9 = (fabs(dieta)<3 && fabs(diphi)<5);
	  //bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<5);
	  bool noForce = (!enforce5x5 && !enforce3x3 && !enforce3x5);
	  bool is5x9 = (fabs(dieta)<3 && fabs(diphi)<7 &&
			fabs(diphi)!=3 && fabs(diphi)!=4);
	  //bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<7 &&
	  //		fabs(diphi)!=2 && fabs(diphi)!=3);
	  bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<5 &&
	  		fabs(diphi)!=2);
	  bool is3x7 = (fabs(dieta)<2 && fabs(diphi)<5);
	  
	  if ((is5x9 && enforce5x5) || (is3x9 && enforce3x3) ||
	      (is3x7 && enforce3x5) || noForce) {
	    
	    esum1 += edet;
	    e1[depth] += edet;

	    // Control plots of pileup subtraction
	    if (abs(ieta)<=16) {
	      h2dpu1_hb->Fill(dieta, diphi);
	      p2dpu1_hb->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else if (abs(ieta)<=23) {
	      h2dpu1_he1->Fill(dieta, diphi);
	      p2dpu1_he1->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else {
	      h2dpu1_he2->Fill(dieta, diphi);
	      p2dpu1_he2->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	  } // is5x9
	} // for idet

	if (updateSingleDepth) {
	  t_eHcal = esum;
	  t_eHcal10 = esum1;
	}
	
	// Verify sums
	if (checkConsistency && fabs(esum-t_eHcal)>1) {
	  cout << "jentry = " << jentry
	       << ", t_ieta = " << t_ieta << ", t_iphi = " << t_iphi
	       << ", esum = " << esum << ", t_eHcal = " << t_eHcal << endl;
	  continue;
	}
	if (checkConsistency && fabs(esum1-t_eHcal10)>1) {
	  cout << "jentry = " << jentry
	       << ", t_ieta = " << t_ieta << ", t_iphi = " << t_iphi
	       << ", esum1 = " << esum1 << ", t_eHcal1 = " << t_eHcal10 << endl;
	  continue;
	}

	// Merge depths 1+2 by replacing them both with half-average
	// This will preserve the array sizes for later
	// => matrix was singular, try setting e[1] to eps and e[2] to sum
	if (mergeDepths1and2) {
	  //double e12_avg = 0.5*(e[1]+e[2]);
	  //double pu12_avg = 0.5*(e1[1]+e1[2]);
	  //e[1] = e[2] = e12_avg;
	  //e1[1] = e1[2] = pu12_avg;
	  double eps = 0.01; // small non-zero energy offset
	  double e12_sum = (e[1]+e[2]);
	  double pu12_sum = (e1[1]+e1[2]);
	  e[1] = eps; e[2] = e12_sum;
	  e1[1] = 0; e1[2] = pu12_sum;
	}
	
	// Fill histograms and profiles per depth
	const int ndepth = 10;
	double vrc[ndepth];
	//double rcsum = rc; // precalculated (or sum explicitly)
	double areaScale(0.5);
	if (enforce5x5) areaScale = 0.5;//(5.*5.)/(5.*9.-5.*5.);
	if (enforce3x3) areaScale = 0.5;//(3.*3.)/(3.*9.-3.*3.);
	if (enforce3x5) areaScale = areaScale3x5(t_ieta);
	double epu1sum = (esum-esum1)*areaScale;
	double rcsum = (esum - epu1sum) / (t_p - t_eMipDR);
	for (int i = 0; i != ndepth; ++i) {

	  double depth = i;
	  double epu1 = (i==0 ? 0 : e1[i]-e[i]) * areaScale;
	  //double epu3 = (t_eHcal30-t_eHcal10)*0.5;
	  double rc = (e[i]-epu1) / (t_p - t_eMipDR);

	  h3raw->Fill(t_ieta, depth, e[i] / (t_p - t_eMipDR));
	  h3pu1->Fill(t_ieta, depth, epu1 / (t_p - t_eMipDR));
	  h3c->Fill(t_ieta, depth, rc);
	  
	  if (rcsum>0.15 && rcsum<1.85) {
	    p2raw->Fill(t_ieta, depth, e[i] / (t_p - t_eMipDR));
	    p2pu1->Fill(t_ieta, depth, epu1 / (t_p - t_eMipDR));
	    p2c->Fill(t_ieta, depth, rc);
	    if (jentry%2==0) p2c_even->Fill(t_ieta, depth, rc);
	    if (jentry%2==1) p2c_odd->Fill(t_ieta, depth, rc);
	    vrc[i] = rc;
	  }
	  else {
	    vrc[i] = 0.;
	  }
	} // for i
	
	// Calculate covariance
	for (int i = 0; i != ndepth; ++i) {
	  for (int j = 0; j != ndepth; ++j) {
	    if (rcsum>0.15 && rcsum<1.85) {
	      p3c->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	      if (jentry%2==0) p3c_even->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	      if (jentry%2==1) p3c_odd->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	    }
	  } // for i
	} // for j

      } // doPerDepth

      
      if (doSingleDepth) {

	double areaScale(0.5);
	if (enforce5x5) areaScale = 0.5;//(5.*5.)/(5.*9.-5.*5.);
	if (enforce3x3) areaScale = 0.5;//(3.*3.)/(3.*9.-3.*3.);
	if (enforce3x5) areaScale = areaScale3x5(t_ieta);
	double epu1 = (t_eHcal10-t_eHcal) * areaScale;
	double epu3 = (t_eHcal30-t_eHcal10)*0.5;
	double rc = (t_eHcal-epu1) / (t_p - t_eMipDR);
	
	h2raw->Fill(t_ieta, t_eHcal / (t_p - t_eMipDR));
	h2pu1->Fill(t_ieta, epu1 / (t_p - t_eMipDR));
	h2pu3->Fill(t_ieta, epu3 / (t_p - t_eMipDR));
	h2mip->Fill(t_ieta, t_eMipDR / (t_p - t_eMipDR));
	h2c->Fill(t_ieta, rc);
	
	if (rc>0.15 && rc<1.85) {
	  praw->Fill(t_ieta, t_eHcal / (t_p - t_eMipDR));
	  ppu1->Fill(t_ieta, epu1 / (t_p - t_eMipDR));
	  ppu3->Fill(t_ieta, epu3 / (t_p - t_eMipDR));
	  pmip->Fill(t_ieta, t_eMipDR / (t_p - t_eMipDR));
	  pc->Fill(t_ieta, rc);
	  if (jentry%2==0) pc_even->Fill(t_ieta, rc);
	  if (jentry%2==1) pc_odd->Fill(t_ieta, rc);
	}
      }

   } // for jentry

   fout->Write();
   fout->Close();

   cout << "\nClosed output file IsoTrack.root" << endl << flush;
   cout << "IsoTrack::Loop() finished.\n" << endl << flush;
} // Loop
