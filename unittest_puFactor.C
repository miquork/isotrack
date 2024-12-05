// Purpose: test puFactor correction from Sunanda

double puFactor(int ieta, double pmom, double eHcal, double ediff, bool debug = false) {
  
  double fac(1.0);
  if (debug)
    std::cout << "Input Type " << 8 << " ieta " << ieta << " pmon " << pmom << " E " << eHcal << ":" << ediff;
  
  int jeta = std::abs(ieta);
  double d2p = (ediff / pmom);
  const double DELTA_CUT = 0.03;
  
  // type == 8: Run3 MAHI. Deleted all the rest
  //} else {  // Mahi 22pu (Jan, 2022)
  const double CONST_COR_COEF[6] = {0.995902, 0.991240, 0.981019, 0.788052, 0.597956, 0.538731};
  const double LINEAR_COR_COEF[6] = {-0.0540563, -0.104361, -0.215936, -0.147801, -0.160845, -0.154359};
  const double SQUARE_COR_COEF[6] = {0, 0, 0.0365911, 0.0161266, 0.0180053, 0.0184295};
  const int PU_IETA_1 = 7;
  const int PU_IETA_2 = 16;
  const int PU_IETA_3 = 25;
  const int PU_IETA_4 = 26;
  const int PU_IETA_5 = 27;
  unsigned icor = (unsigned(jeta >= PU_IETA_1) +
		   unsigned(jeta >= PU_IETA_2) +
		   unsigned(jeta >= PU_IETA_3) +
		   unsigned(jeta >= PU_IETA_4) +
		   unsigned(jeta >= PU_IETA_5));
  double deltaCut = (icor > 2) ? 1.0 : DELTA_CUT;
  if (d2p > deltaCut)
    fac = (CONST_COR_COEF[icor] + LINEAR_COR_COEF[icor] * d2p + SQUARE_COR_COEF[icor] * d2p * d2p);
  if (debug)
    std::cout << " d2p " << d2p << ":" << DELTA_CUT << " coeff " << icor << ":" << CONST_COR_COEF[icor] << ":"
	      << LINEAR_COR_COEF[icor] << ":" << SQUARE_COR_COEF[icor] << " Fac " << fac;
//}
//}
  if (fac < 0 || fac > 1)
    fac = 0;
  if (debug)
    std::cout << " Final factor " << fac << std::endl;
  return fac;
} // puFactor

void unittest_puFactor() {

  double pmom = 50.;
  double etot = 50; // not used

  // Average epu3 / pmom ~ 0.05 in barrel, where
  // double epu3 = (t_eHcal30-t_eHcal10)*0.5;
  // We want ediff = etot3 - etot1
  // double ediff = (0.05*pmom)*2.; // rough barrel
  // Better estimate from ppu3
  TFile *f = new TFile("rootfiles/IsoTrack_lxplus_v11_24CDEFGHI_orig.root","READ");
  assert(f && !f->IsZombie());
  TProfile *ppu3 = (TProfile*)f->Get("ppu3");
  TH1D *hpu3 = ppu3->ProjectionX("hpu3");
  hpu3->Scale(2.);
  
  TH1D *hkpu = new TH1D("hkpu",";i#eta;Offset scale factor",59,-29.5,29.5);
  for (int i = 1; i != hkpu->GetNbinsX()+1; ++i) {
    int ieta = int(hkpu->GetBinCenter(i));
    int j = ppu3->FindBin(ieta);
    double ediff = ppu3->GetBinContent(j)*2.*pmom;
    double ehcal = etot * puFactor(ieta, pmom, etot, ediff);
    // ehcal = esum - k * ediff
    double k = (ediff>0 ? (etot - ehcal) / ediff : 1);
    hkpu->SetBinContent(i, k);
  }

  hkpu->Draw("HIST");
  hpu3->Draw("SAME HIST");
  
} // unittest_puFactor
