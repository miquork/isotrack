#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cctype>

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TSystem.h"

// Include your gain corrections code
#include "gainCorrections.C"

// Use two-space indentation throughout

static std::vector<int> ieta_values() {
  std::vector<int> etas;
  //for (int e=-29; e<=-16; e++) etas.push_back(e);
  //for (int e=16; e<=29; e++) etas.push_back(e);
  for (int e=-29; e<=-1; e++) etas.push_back(e);
  for (int e=1; e<=29; e++) etas.push_back(e);
  return etas;
}

static std::vector<int> depth_values() {
  std::vector<int> depths;
  for (int d=1; d<=7; d++) depths.push_back(d);
  return depths;
}

// Simple parser for the golden JSON file
std::vector<int> loadGoldenRuns(const std::string &goldenJsonFile) {
  std::vector<int> runs;
  std::ifstream jf(goldenJsonFile);
  if (!jf.is_open()) {
    std::cerr << "Could not open " << goldenJsonFile << "\n";
    return runs;
  }

  std::string line;
  while (std::getline(jf, line)) {
    int run;
    // Try to parse line like: "355374": ...
    // Pattern: optional whitespace, double-quote, digits, double-quote, colon
    if (std::sscanf(line.c_str(), " \"%d\":", &run) == 1) {
      runs.push_back(run);
    }
  }

  std::sort(runs.begin(), runs.end());
  return runs;
}

// Function to load luminosity per run from the file
std::map<int,double> loadLuminosity(const std::string &lumiFile) {
  std::map<int,double> runToLumi;
  std::ifstream lf(lumiFile);
  if(!lf.is_open()) {
    std::cerr << "Could not open " << lumiFile << ", proceeding without luminosities.\n";
    return runToLumi;
  }
  std::string line;
  while (std::getline(lf, line)) {
    int run;
    double lumi;
    // Skip lines that don't match the data pattern
    // Data lines start with "| " and have the run number in the format "| run:fill  |"
    // Let's try to parse the run number and the recorded luminosity
    // Sample line:
    // | 378985:9474  | 04/06/24 00:24:34 | 768  | 768  | 0.000251636    | 0.000248158   |
    // We need to extract the run number and the last number (recorded luminosity)
    if (std::sscanf(line.c_str(), "| %d:%*d  |%*[^|]|%*[^|]|%*[^|]|%*[^|]| %lf   |", &run, &lumi) == 2) {
      runToLumi[run] = lumi;
    }
  }
  return runToLumi;
}

// Dummy energy fraction function
double getEnergyFraction(int /*eta*/, int depth) {
  return depth / 7.0;
}

void unittest_gainCorrections() {
  gSystem->Load("libTree"); // Ensure ROOT libs are loaded if needed

  // Initialize corrections
  gainCorrections();

  // Load the golden runs from JSON
  std::string goldenJsonFile = "textfiles/Cert_Collisions2024_378981_386951_Golden.json";
  std::vector<int> runs = loadGoldenRuns(goldenJsonFile);
  if (runs.empty()) {
    std::cerr << "No runs found in " << goldenJsonFile << ". Exiting.\n";
    return;
  }
  std::cout << "Read in " << runs.size() << " golden runs from file '" << goldenJsonFile << "'" << std::endl;

  // Load luminosity per run
  std::string lumiFile = "textfiles/lumi_2024.txt";
  std::map<int,double> runToLumi = loadLuminosity(lumiFile);
  std::cout << "Read in luminosity for " << runs.size() << " runs from file '" << lumiFile << "'" << std::endl;

  auto etas = ieta_values();
  auto depths = depth_values();

  // Create the 1D and 2D profiles
  TProfile *prof_corr = new TProfile("prof_corr","Average Correction vs ieta; ieta; Average Correction", 59, -29.5, 29.5);
  TProfile2D *prof2_corr = new TProfile2D("prof2_corr","Average Correction; ieta; depth", 59, -29.5, 29.5, 7, 0.5, 7.5);

  // Create output file and histograms
  TFile *fout = new TFile("gainCorrections_unitTest_output.root","RECREATE");

  for (auto run : runs) {
    if (run < 379415) continue; // Exclude runs before 379415 (2024B excluded)
    double runLumi = 0.0;
    auto lumiIt = runToLumi.find(run);
    if (lumiIt != runToLumi.end()) {
      runLumi = lumiIt->second;
    } else {
      // If no luminosity info, skip
      //continue;
      // or set default and print error
      runLumi = 1.0;
      std::cerr << "Missing luminosity for run " << runLumi << endl;
    }

    for (auto eta : etas) {
      for (auto depth : depths) {
        double corr = _gainCorrectionRetriever->getCorrection(run, eta, depth);
        if (corr < 0.0) corr = 1.0; // Default if missing

        double eFrac = getEnergyFraction(eta, depth);
        double weight = runLumi * eFrac;

        // Fill the 2D profile
        prof2_corr->Fill((double)eta, (double)depth, corr, weight);

        // For the 1D profile, accumulate over depth
        prof_corr->Fill((double)eta, corr, weight);
      }
    }
  }

  fout->cd();
  prof_corr->Write();
  prof2_corr->Write();
  fout->Close();

  std::cout << "Unit test completed. Results written to gainCorrections_unitTest_output.root\n";

  // Draw the histograms on canvases
  // Since the file is closed, histograms remain in memory and can be drawn

  // Draw the 2D profile
  TCanvas *c1 = new TCanvas("c1", "2D Profile", 800, 600);
  prof2_corr->SetStats(0);
  prof2_corr->GetZaxis()->SetRangeUser(0.70, 1.1);
  prof2_corr->Draw("COLZ");

  // Draw the 1D profile
  TCanvas *c2 = new TCanvas("c2", "1D Profile", 800, 600);
  prof_corr->SetStats(0);
  prof_corr->GetYaxis()->SetRangeUser(0.95, 1.05);
  prof_corr->Draw();

  // Optionally, save the canvases
  c1->SaveAs("pdf/unittest/gainCorrections_2DProfile.pdf");
  c2->SaveAs("pdf/unittest/gainCorrections_1DProfile.pdf");
}
