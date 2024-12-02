// Purpose: Draw average pileup correction per depth per ieta
#ifndef __drawPileupScan_C__
#define __drawPileupScan_C__
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH2D.h"

#include "tdrstyle_mod22.C"

void drawPileupScan() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  gROOT->ProcessLine(".! mkdir pdf/pileup");
    
  //TFile *f = new TFile("IsoTrack_lxplus_v8_24CDEF.root","READ");
  TFile *f = new TFile("IsoTrack.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  TH1D *h = tdrHist("h","#Delta|i#phi|",-12.5,12.5,"#Delta|i#eta|",-12.5,12.5);
  lumi_136TeV = "2024F";
  extraText = "Private";
  //TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  for (int ieta = 1; ieta != 28; ++ieta) {
    //int ieta = 7;
    
    TCanvas *c1 = new TCanvas(Form("c1_%d",ieta),"c1",4*300,2*300);
    c1->Divide(4,2,0,0);
  

    TProfile3D *p3_0, *p3_1, *p3_3;
    p3_0 = (TProfile3D*)f->Get(Form("pileup/p3_0_ieta%d",ieta));  assert(p3_0);
    p3_1 = (TProfile3D*)f->Get(Form("pileup/p3_1_ieta%d",ieta));  assert(p3_1);
    p3_3 = (TProfile3D*)f->Get(Form("pileup/p3_3_ieta%d",ieta));  assert(p3_3);
    
    for (int idepth = 1; idepth != 9; ++idepth) {
      
      c1->cd(idepth);
      
      p3_0->GetZaxis()->SetRange(idepth,idepth);
      p3_1->GetZaxis()->SetRange(idepth,idepth);
      p3_3->GetZaxis()->SetRange(idepth,idepth);
      if (idepth==8) { // all depths
	p3_0->GetZaxis()->SetRange(1,8);
	p3_1->GetZaxis()->SetRange(1,8);
	p3_3->GetZaxis()->SetRange(1,8);
      }

      // NB: yx = "y vs x" or y vertical and x-horizontal
      TProfile2D *p2_0, *p2_1, *p2_3, *p2_3b;
      p2_0 = p3_0->Project3DProfile("yx");
      p2_0->SetName(Form("p2_0_ieta%d_depth%d",ieta,idepth));
      p2_1 = p3_1->Project3DProfile("yx");
      p2_1->SetName(Form("p2_1_ieta%d_depth%d",ieta,idepth));
      p2_3 = p3_3->Project3DProfile("yx");
      p2_3->SetName(Form("p2_3_ieta%d_depth%d",ieta,idepth));
      p2_3b = p3_3->Project3DProfile("yx");
      p2_3b->SetName(Form("p2_3b_ieta%d_depth%d",ieta,idepth));
      
      if (idepth%4==0) gPad->SetRightMargin(0.15);
      gPad->SetLogz();
      
      h->Draw("AXIS");
      //p2_3->Draw("SAME COLZ");
      p2_0->Draw("SAME COLZ");
      //p2_1->Draw("SAME BOX");
      p2_3b->Draw("SAME BOX");
      
      double pmin(0.8e-3), pmax(0.20);
      p2_0->GetZaxis()->SetRangeUser(pmin,pmax);
      p2_1->GetZaxis()->SetRangeUser(pmin,pmax);
      p2_3->GetZaxis()->SetRangeUser(pmin,pmax);
      p2_3b->GetZaxis()->SetRangeUser(pmin,pmax);
      
      TH2D *h2_0 = p2_0->ProjectionXY(Form("h2_0_ieta%d_depth%d",ieta,idepth));
      h2_0->SetContour(1);
      //h2_0->SetLineColor(kBlack);
      h2_0->SetFillColorAlpha(0, 0.5);
      //h2_0->Draw("SAME CONT3");
      
      gPad->RedrawAxis();
      gPad->Update();
    } // for idepth
    
    c1->SaveAs(Form("pdf/pileup/drawPileUpScan_ieta%d.pdf",ieta));
  } // for ieta
    
} // drawPileupScan
#endif 
