#include "IsoTrackV2.h"

// Forward declarations (dont' work here?)
//void rebinProfiles();
//void drawIsoTrack();
//void drawCorrFact();
//void drawCovariance();
//void hybridCorrFact();
//void drawPileup();

//R__LOAD_LIBRARY(IsoTrackV2.C+g)
R__LOAD_LIBRARY(IsoTrackV2_C.so)
R__LOAD_LIBRARY(drawIsoTrack_C.so)
R__LOAD_LIBRARY(compareCorrFact_C.so)
R__LOAD_LIBRARY(drawPileUp_C.so)
R__LOAD_LIBRARY(rebinProfiles_C.so)
R__LOAD_LIBRARY(drawCovariance_C.so)
R__LOAD_LIBRARY(hybridCorrFact_C.so)

void mk_IsoTrackV2(string era = "24F", string version = "local") {

  string path = gSystem->pwd();
  cout << "Current path: " << path << endl << flush;
  //bool runLXPLUS = (path=="/afs/cern.ch/user/v/voutila/scratch0/IsoTrack");
  bool runLXPLUS = TString(path.c_str()).Contains("/afs/cern.ch/");

  cout << "Era "<<era<<", version "<<version<<endl<<flush;

  TString t(era.c_str());
  TChain *c(0);
  if (t.Contains("MC"))
    c = new TChain("hcalIsoTrkAnalyzer/CalibTree");
  else
    c = new TChain("hcalIsoTrackAnalyzer/CalibTree");

  const char *ce = era.c_str();
  const char *cv = version.c_str();
  
  if (runLXPLUS) {
    // Copy files to lxplus and back:
    // rsync -rutP *.C *.h 24FEB.txt voutila@lxplus.cern.ch:~/scratch0/IsoTrack/
    // rsync -rutP voutila@lxplus.cern.ch:~/scratch0/IsoTrack/IsoTrack.root ./IsoTrack_lxplus.root
    cout << "Running on lxplus..." << endl << flush;
    ifstream fin(Form("%sE%s.txt",ce,era=="24F"?"B":"A"));
    string s;
    while (fin >> s) c->AddFile(s.c_str());
  }
  else {
    cout << "Running on a local test file..." << endl << flush;
    c->AddFile("../data/IsoTrack/EGamma0_40to60_v4.root");
  }

  // Create IsoTrack.root output file for interactive analysis (lxplus only)

  // Uncomment on lxplus

  IsoTrackV2 it(c);
  //if (runLXPLUS) // don't redo histograms if running locally
  it.Loop(); // uncomment on lxplus
  //exit(0);

} // mk_IsoTrackV2
