#include "IsoTrack.h"
#include "rebinProfiles.C"
#include "drawIsoTrack.C"
#include "drawCorrFact.C"
#include "drawCovariance.C"
// Forward declarations (dont' work here?)
//void rebinProfiles();
//void drawIsoTrack();
//void drawCorrFact();
//void drawCovariance();

//R__LOAD_LIBRARY(IsoTrack.C+g)
R__LOAD_LIBRARY(IsoTrack_C.so)


void mk_IsoTrack() {

  string path = gSystem->pwd();
  cout << "Current path: " << path << endl << flush;
  //bool runLXPLUS = (path=="/afs/cern.ch/user/v/voutila/scratch0/IsoTrack");
  bool runLXPLUS = TString(path.c_str()).Contains("/afs/cern.ch/user/");
  
  TChain *c = new TChain("hcalIsoTrackAnalyzer/CalibTree");

  if (runLXPLUS) {
    // Copy files to lxplus and back:
    // rsync -rutP *.C *.h 24FEB.txt voutila@lxplus.cern.ch:~/scratch0/IsoTrack/
    // rsync -rutP voutila@lxplus.cern.ch:~/scratch0/IsoTrack/IsoTrack.root ./IsoTrack_lxplus.root
    cout << "Running on lxplus..." << endl << flush;
    ifstream fin("24FEB.txt");
    string s;
    while (fin >> s) c->AddFile(s.c_str());
  }
  else {
    cout << "Running on a local test file..." << endl << flush;
    c->AddFile("../data/IsoTrack/EGamma0_40to60_v4.root");
  }

  // Create IsoTrack.root output file for interactive analysis
  //IsoTrack it(c);
  //it.Loop();

  // Rebin ieta for more stable depths
  gROOT->ProcessLine(".L rebinProfiles.C+g");
  rebinProfiles();
  
  // Solve corrections, IsoTrack.root -> CorrFact*.root
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".L drawIsoTrack.C+g");
  drawIsoTrack();

  // Even/odd tests of uncertainties for CorrFact*.root
  gROOT->ProcessLine(".L drawCorrFact.C+g");
  drawCorrFact();
  drawCorrFact("_wide");

  // Control plots of depth covariance in IsoTrack.root
  gROOT->ProcessLine(".! mkdir pdf/drawCovariance");
  gROOT->ProcessLine(".L drawCovariance.C+g");
  drawCovariance();
} // mk_IsoTrack
