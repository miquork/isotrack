// Purpose: compare two sets of CorrFact.root
#include "TFile.h"
#include "TH1D.h"

#include <fstream>
#include <map>
#include <string>

#include "tdrstyle_mod22.C"

void parseSunanda();

void compareCorrFacts(string file1, string file2,
		      string name1="file1", string name2="filel2",
		      string title = "IsoTrack", string name="",
		      double ymin=0.8, double ymax=1.2, bool addDI=true) {

  cout << "Running compareCorrFact(\""<<file1<<",\""<<file2<<"\")..."
       << endl << flush;
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *fe = new TFile(Form("%s",file1.c_str()),"READ");
  assert(fe && !fe->IsZombie());
  TFile *fo = new TFile(Form("%s",file2.c_str()),"READ");
  assert(fo && !fo->IsZombie());

  cout << "Input files " << fe->GetName() << " and "
       << fo->GetName() << endl << flush;

  curdir->cd();

  const char *cf1 = name1.c_str();
  const char *cf2 = name2.c_str();
  const char *cn = name.c_str();
  TH1D *h_1 = tdrHist(Form("h_1dcf%s%s",cf1,cf2),
		      //Form("CF(%s) / CF(%s)",cf2,cf1), 0.8,1.2, "i#eta",-29,29);
		      Form("CF(%s) / CF(%s)",cf2,cf1),ymin,ymax,"i#eta",-29,29);
  lumi_136TeV = title.c_str();
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1_dcf%s%s",cf1,cf2),h_1,8,11,kRectangular);

  TLegend *leg = tdrLeg(0.40,0.90-(addDI ? 5 : 4)*0.05,0.65,0.90);
  TLegend *leg2 = tdrLeg(0.40,0.15,0.65,0.15+3*0.05);

  TH1D *he = (TH1D*)fe->Get("h_di"); assert(he);
  TH1D *ho = (TH1D*)fo->Get("h_di"); assert(ho);
  TH1 *hr = (TH1D*)ho->Clone(Form("hr%s%s",cf1,cf2));
  hr->Divide(he);

  // Patch the missing negative side by reflecting denominator to -|i#eta|
  if (name=="abs") {
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
      double ix = hr->GetBinCenter(i);
      if (ix<0 && hr->GetBinContent(i)==0) {
	int i2 = hr->FindBin(-ix);
	if (he->GetBinContent(i2)!=0) {
	  hr->SetBinContent(i, ho->GetBinContent(i) / he->GetBinContent(i2));
	  hr->SetBinError(i, hr->GetBinError(i2));
	}
      }
    }
  }

  if (addDI) leg->AddEntry(hr,"Depth independent","PLE");

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
    TH1D *he = (TH1D*)fe->Get(Form("hf_dd_%d",i)); assert(he);
    TH1D *ho = (TH1D*)fo->Get(Form("hf_dd_%d",i)); assert(ho);
    TH1 *hr = (TH1D*)ho->Clone(Form("hfr_dd_%d%s%s",i,cf1,cf2));
    hr->Divide(he);

    // Patch the missing negative side by reflecting denominator to -|i#eta|
    if (name=="abs") {
      for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
	double ix = hr->GetBinCenter(i);
	if (ix<0 && hr->GetBinContent(i)==0) {
	  int i2 = he->FindBin(-ix);
	  if (he->GetBinContent(i2)!=0) {
	    hr->SetBinContent(i, ho->GetBinContent(i) / he->GetBinContent(i2));
	    hr->SetBinError(i, hr->GetBinError(i2));
	  }
	}
      }
    }

    tdrDraw(hr,"Pz",marker[i],color[i], kSolid,-1,kNone,0, 0.8);
    hr->SetLineWidth(2);
    if (i<5) leg->AddEntry(hr,Form("DD depth %d",i),"PLE");
    else leg2->AddEntry(hr,Form("DD depth %d (HE)",i),"PLE");
  }

  // Reference results from depth-independent case
  if (addDI) {
    tdrDraw(hr,"Pz",kFullCircle,kBlack, kSolid,-1,kNone,0, 0.8);
    hr->SetLineWidth(2);
  }
    
  c1->SaveAs(Form("pdf/compareCorrFact_%s_vs_%s%s.pdf",cf1,cf2,cn));

  cout << "Done with compareCorrFact(\""<<file1<<",\""<<file2<<"\")."
       << endl << flush;  
} // compareCorrFacts

void compareCorrFact() {

  /*
  // First and second half of early 2024
  compareCorrFacts("rootfiles/CorrFact_lxplus_v8_24CDE.root",
		   "rootfiles/CorrFact_lxplus_v8_24F.root",
		   "24CDE","24F","IsoTrack");

  // First and second half of early 2024 with wider bins
  compareCorrFacts("rootfiles/CorrFact_wide_lxplus_v8_24CDE.root",
		   "rootfiles/CorrFact_wide_lxplus_v8_24F.root",
		   "24CDE","24F","IsoTrack wide bins","wide");

  // Plus vs minus (average of plus and minus) for early 2024 
  compareCorrFacts("rootfiles/CorrFact_abs_lxplus_v8_24CDEF.root",
		   "rootfiles/CorrFact_lxplus_v8_24CDEF.root",
		   "Both","Signed","IsoTrack |i#eta| bins","abs");
  */

  parseSunanda();

  // v11: useClassic
  // v17: useSunanda (w/ fixed esum3)
  // v18: useSunanda w/ CalibCorr.C gains + phi
  // v19: useSunanda w/ CalibCorr.C gains + phi + cuts (only worked for HE)
  // v20: useSunanda w/ CalibCorr.C gains + phi + cuts (year 3; fixed for HB)
  // v21: useSunanda w/ CalibCorr.C gains + phi + cuts (year 3->2)
  // v22: useSunanda w/ CalibCorr.C gains + phi + cuts=2 + puFactor(-8)
  // v23: useSunanda w/ CalibCorr.C gains + phi + cuts=3 + puFactor(-8)
  // v24: useSunanda w/ CalibCorr.C gains + phi + cuts=3 + puFactor(+8) [== v20]
  // v25: useSunanda w/ CalibCorr.C cuts=3 + puFactor(+8) only
  compareCorrFacts("rootfiles/CorrFact_Sunanda_24CDEFGHI.root",
		   "rootfiles/CorrFact_lxplus_v25_24CDEFGHI.root",
		   "Sunanda","Mikko","IsoTrack: teams","_zoomout",0.,2.,false);
  compareCorrFacts("rootfiles/CorrFact_Sunanda_24CDEFGHI.root",
		   "rootfiles/CorrFact_lxplus_v25_24CDEFGHI.root",
		   "Sunanda","Mikko","IsoTrack: teams","_zoomin",0.8,1.3,false);
  /*
  compareCorrFacts("rootfiles/CorrFact_lxplus_v13_24CDEFGHI.root",
		   "rootfiles/CorrFact_lxplus_v14_24CDEFGHI.root",
		   "v13_regular","v14_3x5","IsoTrack: containment 3#times5",
		   "",0.8,1.6);

  compareCorrFacts("rootfiles/CorrFact_lxplus_v14_24CDEFGHI.root",
		   "rootfiles/CorrFact_hybrid_lxplus_v14_24CDEFGHI.root",
		   "signed","abs","IsoTrack: hybrid correction",
		   "",0.9,1.15);

  compareCorrFacts("rootfiles/CorrFact_odd_hybrid_lxplus_v14_24CDEFGHI.root",
		   "rootfiles/CorrFact_even_hybrid_lxplus_v14_24CDEFGHI.root",
		   "odd","even","IsoTrack: even vs odd",
		   "",0.9,1.15);
  */

  /*
  compareCorrFacts("rootfiles/CorrFact_lxplus_v15_24CDEFGHI.root",
		   "rootfiles/CorrFact_lxplus_v14_24CDEFGHI.root",
		   "v15","v14","IsoTrack: v15 with gain corrections",
		   "",0.8,1.3);
  */
}

// Parse Sunanda's text files into .root 
void parseSunanda() {

  // For example histogram
  TFile *f = new TFile("rootfiles/CorrFact_lxplus_v11_24CDEFGHI.root","READ");
  assert(f && !f->IsZombie());

  // For output histograms
  TFile *fout = new TFile("rootfiles/CorrFact_Sunanda_24CDEFGHI.root","RECREATE");
  assert(fout && !fout->IsZombie());
  
  TH1D *h = (TH1D*)f->Get("h_dd_1"); assert(h);
  TH1D *h0 = (TH1D*)h->Clone("hf_dd_0"); h0->Reset();
  
  // For Sunanda's results in text file format
  ifstream fin("textfiles/24CDEFGHIEBS00corr.txt");

  // Remove header
  char cline[1024];
  fin.getline(cline,1024);
  cout << cline << endl << flush;

  map<int, TH1D*> mh;
  mh[0] = h0;
  
  while (fin.getline(cline,1024)) {

    char detid[512];
    int ieta, depth;
    double corr, err;
    if (sscanf(cline,"%s %d %d %lf %lf",detid,&ieta,&depth,&corr,&err)==5) {
      TH1D *h = mh[depth];
      if (!h) {
	h = (TH1D*)h0->Clone(Form("hf_dd_%d",depth));
	mh[depth] = h;
      }
      int i = h->FindBin(ieta);
      h->SetBinContent(i, corr);
      h->SetBinError(i, err);
    }
    else {
      cout << "Error reading line:"<<endl<<cline<<endl<<flush;
    }
    
  } // while fin

  fout->cd();
  if (mh[0]) mh[0]->Write("h_di",TObject::kOverwrite);
  for (unsigned int i = 0; i != mh.size(); ++i) {
    TH1D *h = mh[i];
    if (!h) continue;
    h->Write(Form("hf_dd_%d",i),TObject::kOverwrite);
  }
  fout->Close();
  
} // parseSunanda
