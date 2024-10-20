// Purpose: compare even and odd corrfact.root
#include "TFile.h"
#include "tdrstyle_mod22.C"

#include <iostream>

void drawCorrFact(string mode = "") {

  cout << "Running drawCorrFact(\""<<mode<<"\")..." << endl << flush;
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();

  const char *cm = mode.c_str();
  TFile *fe = new TFile(Form("CorrFact_even%s.root",cm),"READ");
  //TFile *fe = new TFile("CorrFact_even.root","READ");
  //TFile *fe = new TFile("CorrFact_even_wide.root","READ");
  assert(fe && !fe->IsZombie());
  TFile *fo = new TFile(Form("CorrFact_odd%s.root",cm),"READ");
  //TFile *fo = new TFile("CorrFact_odd.root","READ");
  //TFile *fo = new TFile("CorrFact_odd_wide.root","READ");
  assert(fo && !fo->IsZombie());

  cout << "Input files " << fe->GetName() << " and "
       << fo->GetName() << endl << flush;

  // Sunanda's refered
  //TFile *fs = new TFile("../IsoTrackSunanda/uncX/MFitC.root","READ");
  //TFile *fs = new TFile("../IsoTrackSunanda/uncX/MFitF.root","READ");
  //TFile *fs = new TFile("../IsoTrackSunanda/uncX/SFitC.root","READ");
  //TFile *fs = new TFile("../IsoTrackSunanda/uncX/SFitF.root","READ");
  //assert(fs && !fs->IsZombie());

  curdir->cd();

  TH1D *h_1 = tdrHist(Form("h_1dcf%s",cm),"CF_{odd} / CF_{even}",0.8,1.2,
		      "i#eta",-29,29);
  lumi_136TeV = "2024F EGamma";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1_dcf%s",cm),h_1,8,11,kRectangular);

  TLegend *leg = tdrLeg(0.40,0.90-5*0.05,0.65,0.90);
  TLegend *leg2 = tdrLeg(0.40,0.15,0.65,0.15+3*0.05);//4*0.05);

  TH1D *he = (TH1D*)fe->Get("h_di"); assert(he);
  TH1D *ho = (TH1D*)fo->Get("h_di"); assert(ho);
  TH1 *hr = (TH1D*)ho->Clone(Form("hr%s",cm));
  hr->Divide(he);

  // Reference results from Sunanda
  //TH1D *hs = (TH1D*)fs->Get("24F12CM60Z3"); assert(hs);
  //TH1D *hs = (TH1D*)fs->Get("24F12EAM60Z3"); assert(hs);
  //TH1D *hs = (TH1D*)fs->Get("24F12CS60Z3"); assert(hs);
  //TH1D *hs = (TH1D*)fs->Get("24F12EAS60Z3"); assert(hs);
  TH1D *hb = (TH1D*)he->Clone(Form("hb%s",cm));
  hb->Add(ho);
  hb->Scale(0.5);
  //TH1 *hr2 = (TH1D*)hs->Clone("hr2");
  //hr2->Divide(hb);
  /*
  // Do this a bit harder way due to different number of bins
  for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
    int j = hb->FindBin(hs->GetBinCenter(i));
    if (hb->GetBinContent(j)!=0) {
      hr2->SetBinContent(i, hs->GetBinContent(i)/hb->GetBinContent(j));
      hr2->SetBinError(i, hs->GetBinError(i)/hb->GetBinContent(j));
    }
    else {
      hr2->SetBinContent(i, 0.);
      hr2->SetBinError(i, 0.);
    }
  }
  */

  leg->AddEntry(hr,"Depth independent","PLE");
  //leg2->AddEntry(hr2,"Sunanda's SFitF","PLE");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlue);
  l->DrawLine(-29,1,29,1);
  l->SetLineStyle(kDotted);
  l->SetLineColor(kGray+1);
  l->DrawLine(-15.5,0.8,-15.5,1.2);
  l->DrawLine(+15.5,0.8,+15.5,1.2);
  l->SetLineColor(kGray+1);
  l->DrawLine(-17.5,0.8,-17.5,1.2);
  l->DrawLine(+17.5,0.8,+17.5,1.2);
  
  // Calculate correction ratios per depth
  const int ndepth = 10;
  int color[ndepth] =
    {kBlack, kRed, kGreen+2, kBlue, kMagenta+1,
     kOrange+1,kCyan+2,kGray+2, 1,1};
  int marker[ndepth] =
    {1, kOpenTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp,
     kOpenStar, kOpenDiamond, kOpenCross, 1,1};

  cout << "Drawing plots per depth" << endl << flush;
  
  // Get even+odd, do ratio, plot
  for (int i = 1; i != 8; ++i) {
    /*
    TH1D *he = (TH1D*)fe->Get(Form("h_dd_%d",i)); assert(he);
    TH1D *ho = (TH1D*)fo->Get(Form("h_dd_%d",i)); assert(ho);
    TH1 *hr = (TH1D*)ho->Clone(Form("hr_dd_%d",i));
    */
    TH1D *he = (TH1D*)fe->Get(Form("hf_dd_%d",i)); assert(he);
    TH1D *ho = (TH1D*)fo->Get(Form("hf_dd_%d",i)); assert(ho);
    TH1 *hr = (TH1D*)ho->Clone(Form("hfr_dd_%d%s",i,cm));
    hr->Divide(he);

    tdrDraw(hr,"Pz",marker[i],color[i], kSolid,-1,kNone,0, 0.8);
    if (i<5) leg->AddEntry(hr,Form("DD depth %d",i),"PLE");
    else leg2->AddEntry(hr,Form("DD depth %d (HB)",i),"PLE");
  }

  // Reference results from depth-independent case
  tdrDraw(hr,"Pz",kFullCircle,kBlack, kSolid,-1,kNone,0, 0.8);
  hr->SetLineWidth(2);

  // Reference result from Sunanda
  //tdrDraw(hr2,"Pz",kOpenCircle,kBlack, kSolid,-1,kNone,0, 1.2);
  //hr2->SetLineWidth(2);

  c1->SaveAs(Form("pdf/IsoTrack_CFRatio%s.pdf",cm));

  cout << "Done with drawCorrFact()." << endl << flush;
  
} // drawCorrFact
