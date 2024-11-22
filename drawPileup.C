// Purpose: Draw average pileup correction 
#ifndef __drawPileup_C__
#define __drawPileup_C__
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"

void drawPileup() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  //TFile *f = new TFile("IsoTrack_lxplus_v8_24CDEF.root","READ");
  TFile *f = new TFile("IsoTrack.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  //////////////////////////////
  // First, 1D case as warmup //
  //////////////////////////////
  TProfile *praw = (TProfile*)f->Get("praw"); assert(praw);
  TProfile *ppu1 = (TProfile*)f->Get("ppu1"); assert(ppu1);
  TProfile *ppu3 = (TProfile*)f->Get("ppu3"); assert(ppu3);
  TProfile *pmip = (TProfile*)f->Get("pmip"); assert(pmip);
  TProfile *pc = (TProfile*)f->Get("pc"); assert(pc);

  // Adjust big ring to match small ring better
  // Areas presumably 3x3(=9) for core
  // 5x5(=25) for small
  TH1D *hpu3 = ppu3->ProjectionX("hpu3");
  //hpu3->Scale(
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  
  TH1D *h = tdrHist("h","E_{calo} / (p_{track} - E_{MIP})",0,2.8,
		    "i#eta",-27.5,27.5);
  //lumi_136TeV = "2024CDEF";
  lumi_136TeV = "2024F";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  h->GetYaxis()->SetTitleOffset(1.2);
  
  //TLegend *leg = tdrLeg(0.35,0.90-0.05*06,0.60,0.90);
  TLegend *leg = tdrLeg(0.35,0.90-0.05*5,0.60,0.90);
  leg->SetHeader("40 #leq p_{track} < 60 GeV");
  
  tdrDraw(praw,"HIST",kNone,kBlue,kSolid,-1,1001,kBlue-9);
  praw->SetFillColorAlpha(kBlue-9,0.7);
  l->DrawLine(-27.5,1,+27.5,1);
  
  tdrDraw(ppu1,"P",kFullSquare,kMagenta+2,kSolid,-1,kNone,0,0.6);
  //tdrDraw(ppu3,"P",kOpenSquare,kRed,kSolid,-1,kNone,0,0.6);
  //tdrDraw(hpu3,"P",kOpenSquare,kRed,kSolid,-1,kNone,0,0.6);
  tdrDraw(pmip,"HIST",kNone,kBlue,kSolid,-1,kNone,0,0.6);
  tdrDraw(pc,"P",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);

  leg->AddEntry(praw,"Raw energy","F");
  
  leg->AddEntry(pc,"Corrected energy","PLE");
  //leg->AddEntry(hpu3,"Offset (big cone)","PLE");
  //leg->AddEntry(ppu1,"Offset (small cone)","PLE");
  leg->AddEntry(ppu1,"Offset","PLE");
  leg->AddEntry(pmip,"MIP energy (ECAL)","L");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/IsoTrack_drawPileup.pdf");

  
  ////////////////////////////////////////////////////////////////////////////
  // Draw 2D distribution of RecHits around track to verify size of PU cone //
  ////////////////////////////////////////////////////////////////////////////

  TProfile2D *p2draw_hb = (TProfile2D*)f->Get("p2draw_hb");
  assert(p2draw_hb);
  TProfile2D *p2dpu1_hb = (TProfile2D*)f->Get("p2dpu1_hb");
  assert(p2dpu1_hb);
  TH2D *h2draw_hb = p2draw_hb->ProjectionXY("h2draw1_hb");
  TH2D *h2dpu1_hb = p2dpu1_hb->ProjectionXY("h2dpu1_hb");
  TH2D *h2dpu0_hb = p2dpu1_hb->ProjectionXY("h2dpu0_hb");
  h2dpu0_hb->Add(h2draw_hb,-1);

  TProfile2D *p2draw_he1 = (TProfile2D*)f->Get("p2draw_he1");
  assert(p2draw_he1);
  TProfile2D *p2dpu1_he1 = (TProfile2D*)f->Get("p2dpu1_he1");
  assert(p2dpu1_he1);
  TH2D *h2draw_he1 = p2draw_he1->ProjectionXY("h2draw1_he1");
  TH2D *h2dpu1_he1 = p2dpu1_he1->ProjectionXY("h2dpu1_he1");
  TH2D *h2dpu0_he1 = p2dpu1_he1->ProjectionXY("h2dpu0_he1");
  h2dpu0_he1->Add(h2draw_he1,-1);
  
  TProfile2D *p2draw_he2 = (TProfile2D*)f->Get("p2draw_he2");
  assert(p2draw_he2);
  TProfile2D *p2dpu1_he2 = (TProfile2D*)f->Get("p2dpu1_he2");
  assert(p2dpu1_he2);
  TH2D *h2draw_he2 = p2draw_he2->ProjectionXY("h2draw1_he2");
  TH2D *h2dpu1_he2 = p2dpu1_he2->ProjectionXY("h2dpu1_he2");
  TH2D *h2dpu0_he2 = p2dpu1_he2->ProjectionXY("h2dpu0_he2");
  h2dpu0_he2->Add(h2draw_he2,-1);
  
  TCanvas *c2 = new TCanvas("c2","c2",3*300,3*300);
  c2->Divide(3,3,0,0);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);

  //double pmin(0.20e-3), pmax(0.20);
  // pmin is greater than 0.2 GeV (calo) / 60 GeV (track) = 3.3e-3 
  double pmin(0.8e-3), pmax(0.20);
  c2->cd(1);
  gPad->SetLogz();
  h2dpu1_hb->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu1_hb->UseCurrentStyle();
  h2dpu1_hb->SetYTitle("#Delta|i#phi|");
  h2dpu1_hb->Draw("COL");
  tex->DrawLatex(0.60,0.92,"HB |i#eta|#leq16 all");

  c2->cd(2);
  gPad->SetLogz();
  h2dpu1_he1->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu1_he1->UseCurrentStyle();
  h2dpu1_he1->Draw("COL");
  tex->DrawLatex(0.40,0.92,"HE1 17#leq|i#eta|#leq24 all");

  c2->cd(3);
  gPad->SetLogz();
  gPad->SetRightMargin(0.20);
  h2dpu1_he2->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu1_he2->UseCurrentStyle();
  h2dpu1_he2->Draw("COLZ");
  tex->DrawLatex(0.20,0.92,"HE2 25#leq|i#eta|#leq27 all");

  
  c2->cd(4);
  gPad->SetLogz();
  h2draw_hb->GetZaxis()->SetRangeUser(pmin,pmax);
  h2draw_hb->UseCurrentStyle();
  h2draw_hb->SetYTitle("#Delta|i#phi|");
  h2draw_hb->Draw("COL");
  tex->DrawLatex(0.60,0.92,"HB |i#eta|#leq16 core");

  c2->cd(5);
  gPad->SetLogz();
  h2draw_he1->GetZaxis()->SetRangeUser(pmin,pmax);
  h2draw_he1->UseCurrentStyle();
  h2draw_he1->Draw("COL");
  tex->DrawLatex(0.40,0.92,"HE1 17#leq|i#eta|#leq24 core");

  c2->cd(6);
  gPad->SetLogz();
  gPad->SetRightMargin(0.20);
  h2draw_he2->GetZaxis()->SetRangeUser(pmin,pmax);
  h2draw_he2->UseCurrentStyle();
  h2draw_he2->Draw("COLZ");
  tex->DrawLatex(0.20,0.92,"HE2 25#leq|i#eta|#leq27 core");

  
  c2->cd(7);
  gPad->SetLogz();
  h2dpu0_hb->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu0_hb->UseCurrentStyle();
  h2dpu0_hb->SetXTitle("#Delta|i#eta|");
  h2dpu0_hb->SetYTitle("#Delta|i#phi|");
  h2dpu0_hb->Draw("COL");
  tex->DrawLatex(0.60,0.92,"HB |i#eta|#leq16 PU");

  c2->cd(8);
  gPad->SetLogz();
  h2dpu0_he1->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu0_he1->UseCurrentStyle();
  h2dpu0_he1->SetXTitle("#Delta|i#eta|");
  h2dpu0_he1->Draw("COL");
  tex->DrawLatex(0.40,0.92,"HE1 17#leq|i#eta|#leq24 PU");

  c2->cd(9);
  gPad->SetLogz();
  gPad->SetRightMargin(0.20);
  h2dpu0_he2->GetZaxis()->SetRangeUser(pmin,pmax);
  h2dpu0_he2->UseCurrentStyle();
  h2dpu0_he2->SetXTitle("#Delta|i#eta|");
  h2dpu0_he2->Draw("COLZ");
  tex->DrawLatex(0.20,0.92,"HE2 25#leq|i#eta|#leq27 PU");
  
  c2->SaveAs("pdf/IsoTrack_drawPileup_etaphi.pdf");
} // drawPileup
#endif
