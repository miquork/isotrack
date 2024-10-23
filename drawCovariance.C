//Purpose: draw covariance between energy measurements at different depths
#ifndef __drawCovariance_C__
#define __drawCovariance_C__
#include "TFile.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include <iostream>

#include "tdrstyle_mod22.C"

void drawCovariances(int ieta, string era, string version);

void drawCovariance(int ieta = 0, string era="24F", string version="vX") {

  cout << "Running drawCovariance("<<ieta<<",\""<<era
       << ",\""<<version<<"\")..."<<endl<<flush;
  
  if (ieta==0) {
    for (int i = -27; i != +27+1; ++i) {
      drawCovariances(i, era, version);
    }
  }
  else
    drawCovariances(ieta, era, version);

  cout << "drawCovariance("<<ieta<<") finished.\n" << endl << flush;
}

void drawCovariances(int ieta, string era, string version) {

  //cout << "drawCovariances("<<ieta<<")" << endl << flush;
  const char *ce = era.c_str();
  const char *cv = version.c_str();
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile(Form("IsoTrack_%s_%s.root",cv,ce),"READ");
  //TFile *f = new TFile("IsoTrack_lxplus_v4.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  if (ieta==-27) cout << "Input from " << f->GetName() << endl << flush;

  TProfile3D *p3 = (TProfile3D*)f->Get("p3c"); assert(p3);
  
  int ix = p3->GetXaxis()->FindBin(ieta);
  p3->GetXaxis()->SetRange(ix, ix);
  TProfile2D *p2 = (TProfile2D*)p3->Project3DProfile("yz");
  p2->SetName(Form("p2yz_cov_%d",ieta));
  TH2D *h2 = p2->ProjectionXY(); h2->Reset();
  
  for (int i = 1; i != p2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != p2->GetNbinsY()+1; ++j) {
      double xy = p2->GetBinContent(i, j);
      double xx = p2->GetBinContent(i, i);
      double yy = p2->GetBinContent(j, j);
      if (xx!=0 && yy!=0)
	h2->SetBinContent(i, j, xy/sqrt(xx*yy));
      else
	h2->SetBinContent(i, j, 0.);
    } // for j
  } // for i
  
  TH1D *h_1 = tdrHist(Form("h_1_%d",ieta),"depth",-0.5,7.5,"depth",-0.5,7.5);
  lumi_136TeV = Form("20%s EGamma",ce);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1_%d",ieta),h_1,8,11,kSquare);
  h_1->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.22);

  h2->Draw("SAME COLZ");
  h2->GetZaxis()->SetRangeUser(-1,1);
  h2->GetZaxis()->SetTitle("Correlation E[xy] / (E[x]#timesE[y])");
  h2->GetZaxis()->SetTitleOffset(1.4);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0.5,0.5,0.5,7.5);
  l->DrawLine(0.5,0.5,7.5,0.5);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.25,0.87,Form("i#eta = %d",ieta));//[%1.0f,%1.0f]",ix1,);
  
  gPad->RedrawAxis();
  gPad->Update();
  c1->SaveAs(Form("pdf/drawCovariance/%s_%s/IsoTrack_Covariance_"
		  "ieta%s%02d_%s_%s.pdf", cv, ce,
		  (ieta>0 ? "p" : "m"), abs(ieta), cv, ce));
} // drawCovariance
#endif
