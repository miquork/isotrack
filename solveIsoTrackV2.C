// Purpose: Construct and solve linear system of HcalRespCorr
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"

#include <iostream>

// Map between channel index and (band,width,depth)
void decodeChannel(const int channel, int &iband, int &iwidth, int &idepth) {
  const int nband = 2; // (core=0,side=1)
  const int nwidth = 3; // (delta_ieta=-1,0,1) plus 1
  const int ndepth = 7; // (depths 1-7) minus 1
  const int nch = nband*nwidth*ndepth;
  assert(channel>=0 && channel<nch);
  
  iband = channel/(nwidth*ndepth);
  int rest = channel%(nwidth*ndepth);
  iwidth = rest/ndepth;
  idepth = rest%ndepth;
}

void solveIsoTrackV2() {

  TFile *f = new TFile("IsoTrackV2.root","READ");
  assert(f && !f->IsZombie());

  TProfile2D *p2 = (TProfile2D*)f->Get("p2"); assert(p2);
  TProfile3D *p3 = (TProfile3D*)f->Get("p3"); assert(p3);

  // Store results back into arrays for cleaner mapping
  // Number of "channels" for each ieta to keep track of:
  const int nband = 2; // (core=0,side=1)
  const int nwidth = 3; // (delta_ieta=-1,0,1) plus 1
  const int ndepth = 7; // (depths 1-7) minus 1
  const int nch = nband*nwidth*ndepth;
  const int neta = 58; // ieta bins
  const int nc = neta*ndepth; // correction factors per RecHit (ieta,depth)

  //////////////////////////////////////////////////////////////////////
  // Read profile data into arrays for bit easier(?) index gymnastics //
  //////////////////////////////////////////////////////////////////////
  
  double esum[neta][nband][nwidth][ndepth];
  double esume[neta][nband][nwidth][ndepth];
  double eprod[neta][nch][nch];
  double eprode[neta][nch][nch];
  assert(p2->GetNbinsX()==neta);
  assert(p2->GetNbinsY()==nch);
  assert(p3->GetNbinsX()==neta);
  assert(p3->GetNbinsY()==nch);
  assert(p3->GetNbinsZ()==nch);
  for (int ieta = 0; ieta != neta; ++ieta) {

    for (int ch = 0; ch != nch; ++ch) {
      int band, width, depth;
      decodeChannel(ch, band, width, depth);
      esum[ieta][band][width][depth] = p2->GetBinContent(ieta+1, ch+1);
      esume[ieta][band][width][depth] = p2->GetBinError(ieta+1, ch+1);

      for (int ch2 = 0; ch2 != nch; ++ch2) {
	eprod[ieta][ch][ch2] = p3->GetBinContent(ieta+1, ch+1, ch2+1);
	eprode[ieta][ch][ch2] = p3->GetBinError(ieta+1, ch+1, ch2+1);
      } // for ch2
    } // for ch
  } // for ieta
  
  // Accessing esum[neta][nch] view
  double (*esum_flat)[nch] = reinterpret_cast<double (*)[nch]>(esum);
  
  // Covariances
  double as[neta][nch][nch];
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int ch1 = 0; ch1 != nch; ++ch1) {
      for (int ch2 = 0; ch2 != nch; ++ch2) {
	
	double Ex = esum_flat[ieta][ch1];
	double Ey = esum_flat[ieta][ch1];
	double Exy = eprod[ieta][ch1][ch2];
	as[ieta][ch1][ch2] = Exy - Ex*Ey;
      } // for ch2
    } // for ch1
  } // for ieta

  
  /////////////////////////////////////////////////////////////////////////
  // 2024-12-19: index gymnastics overwhelmed me, but o1  managed this:  //
  // https://chatgpt.com/share/67641149-10b4-800f-9461-cd80c7197f5f      //
  ////////////////////////////////////////////////////////////////////////

  // Assume these are defined/filled:
  // double esum[neta][nband][nwidth][ndepth];
  // double eprod[neta][nch][nch];
  // constants: nch = nband*nwidth*ndepth

  const int M = nc + neta; // total system size
  TMatrixD A(M, M);
  TVectorD b(M);
  A.Zero();
  b.Zero();

  // Helper lambda to get c-index:
  auto cIndex = [&](int ieta, int depth){ return ieta*ndepth + depth; };
  auto lambdaIndex = [&](int ieta){ return nc + ieta; };

  // We must map from (band,width,depth) to global depth index and ieta offset
  // ieta offsets: width in {0,1,2} corresponds to ieta offsets { -1,0,+1 }
  // channel index: channel = iband*nwidth*ndepth + iwidth*ndepth + idepth
  
  for (int i = 0; i < neta; ++i) {
    // Constraint equation: sum over alpha: E[i,alpha]*c(...) = 1
    // Loop over all channels alpha:
    for (int iband=0; iband<nband; iband++) {
      for (int iwidth=0; iwidth<nwidth; iwidth++) {
	int ieta_off = i + (iwidth - 1); // iwidth=0->i-1,1->i,2->i+1
	if (ieta_off<0 || ieta_off>=neta) continue; // outside range
	for (int idepth=0; idepth<ndepth; idepth++) {
	  int alpha = iband*nwidth*ndepth + iwidth*ndepth + idepth;
	  double E = esum[i][iband][iwidth][idepth]; // <E_i/p>
	  
	  int sign = (iband==0) ? +1 : -1;
	  int ci = cIndex(ieta_off, idepth);
	  A(lambdaIndex(i), ci) += sign*E; 
	}
      }
    }
    b[lambdaIndex(i)] = 1.0;
  }

  // Now the equations for partial derivatives w.r.t. c_{j,d}:
  // For each c_{j,d}, we form:
  for (int j=0; j<neta; j++) {
    for (int d=0; d<ndepth; d++) {
      int row = cIndex(j,d);

      // We want: sum_i sum_alpha sum_beta eprod[i][alpha][beta]*c(...) + sum_i lambda_i * E[i,alpha] = 2*sum_i E[i,alpha]
      // We'll accumulate these terms:

      // Accumulators for RHS:
      double rhs = 0.0;
      
      // Loop over all i (track ieta bins):
      for (int i=0; i<neta; i++) {
	// For each i, we must consider all channels that depend on c_{j,d}.
	// c_{j,d} appears in channels where the ieta of channel matches j and depth matches d.
	// channel indexing: we must check each alpha and see if it corresponds to (j,d).
	for (int iband=0; iband<nband; iband++) {
	  for (int iwidth=0; iwidth<nwidth; iwidth++) {
	    int ieta_off = i + (iwidth - 1);
	    if (ieta_off != j) continue; // only channels that point to ieta=j
	    if (ieta_off<0 || ieta_off>=neta) continue;
	    for (int idepth=0; idepth<ndepth; idepth++) {
	      int alpha = iband*nwidth*ndepth + iwidth*ndepth + idepth;
	      if (idepth != d) continue; // must match the depth d we differentiate w.r.t.

	      int signAlpha = (iband == 0) ? +1 : -1;
	      double Ealpha = esum[i][iband][iwidth][idepth];
	      rhs += 2.0 * (signAlpha*Ealpha); // from the -2 Ealpha*c + ... rearranged form
	      
	      // Now fill matrix elements involving other c_{...} and lambdas:
	      // sum_beta eprod[i][alpha][beta]*c(...) term:
	      for (int beta=0; beta<nch; beta++) {

		double val = eprod[i][alpha][beta]; 

		int iband2, iwidth2, idepth2;
		decodeChannel(beta, iband2, iwidth2, idepth2);
		
		int ieta_off2 = i + (iwidth2 - 1);
		if (ieta_off2<0 || ieta_off2>=neta) continue;
		int ccol = cIndex(ieta_off2, idepth2);
		int signBeta = (iband2 == 0) ? +1 : -1;
		A(row, ccol) += signAlpha*signBeta*val;
	      }
	      
	      // sum_i lambda_i * E[i,alpha] term:
	      // For each lambda_i corresponding to this i:
	      A(row, lambdaIndex(i)) += signAlpha*Ealpha;
	    }
	  }
	}
      }
      
      b[row] = rhs;
    }
  }

  // Now we have a (M x M) linear system: A * X = b
  // The solution is in b[nc] now
  // The lambda_i are in b[nc + ieta], if needed.

  
  ////////////////////////////////////////////////////////////
  // Solution failed because not all depths are at all ieta //
  // We need to compress them out first                     //
  ////////////////////////////////////////////////////////////
  
  // A and b are defined: TMatrixD A(M,M); TVectorD b(M);
  // M = nc + neta, for example.
  std::vector<int> oldToNew(M, -1);
  std::vector<int> newToOld; 
  newToOld.reserve(M);

  // Identify non-empty rows:
  for (int i = 0; i < M; ++i) {
    bool nonEmpty = false;
    for (int j = 0; j < M; ++j) {
      if (A(i,j) != 0.0) {
	nonEmpty = true;
	break;
      }
    }
    if (nonEmpty) {
      oldToNew[i] = (int)newToOld.size();
      newToOld.push_back(i);
    }
  }

  // Build reduced system
  int Mred = (int)newToOld.size();
  TMatrixD Ared(Mred, Mred);
  TVectorD bred(Mred);

  for (int irow = 0; irow < Mred; ++irow) {
    int old_i = newToOld[irow];
    bred[irow] = b[old_i];
    for (int icol = 0; icol < Mred; ++icol) {
      int old_j = newToOld[icol];
      Ared(irow, icol) = A(old_i, old_j);
    }
  }

  // Now solve Ared * Xred = bred
  TDecompLU solverRed(Ared);
  bool okRed = solverRed.Solve(bred);
  if (!okRed) {
    cout << "Reduced solve failed!" << endl;
  }
  else {
    cout << "Reduced solve succeeded" << endl;
  }

  // Map solution back to full size:
  TVectorD fullX(M);
  fullX.Zero();
  for (int irow = 0; irow < Mred; ++irow) {
    fullX[newToOld[irow]] = bred[irow];
  }

  // fullX now contains the solution for the previously non-empty rows.
  // The empty-row variables remain zero.

  
  /////////////////////////////////////
  // Save results for later plotting //
  /////////////////////////////////////
  
  TFile *fout = new TFile("solveIsoTrackV2.root","RECREATE");
  assert(fout && !fout->IsZombie());
  fout->cd();
  int color[ndepth] = {kBlue, kOrange+1, kGreen+1, kRed,
		       kYellow+1, kOrange-7, kGray+2};
  for (int depth = 0; depth != ndepth; ++depth) {
    TH1D *h = new TH1D(Form("hf_dd_%d",depth+1),";ieta;CorrFactV2",
		       61,-30.5,30.5);
    h->SetMarkerColor(color[depth]);
    h->SetLineColor(color[depth]);
    //for (int i = 0; i != M; ++i) {
    for (int i = 0; i != nc; ++i) {
      // int i = ieta*ndepth+depth
      if (i%ndepth!=depth) continue;
      int ieta = i/ndepth; // 0-57
      int jeta = ieta - 29; // -29-28
      if (jeta>=0) jeta += 1; // -29 to +29 with 0 excluded
      double val = fullX[i];
      int j = h->FindBin(jeta);
      h->SetBinContent(j, val);

      // Let's make this simple: corr uncertainty is dominated by
      // the relative uncertainty in core ieta at same depth
      // scaled up by fraction of energy in that ieta,depth over total energy
      // => not quite working for depth 6, 4, maybe 1
      //double ecore = esum[ieta][0][1][depth];
      //double err = esume[ieta][0][1][depth];
      // => try using full 3x5 for estimate of the depth
      // => better, although still large for depth 7 and maybe small for 2,3
      double ecore = (esum[ieta][0][0][depth] + esum[ieta][0][1][depth] + esum[ieta][0][2][depth]);
      double err = sqrt(pow(esume[ieta][0][0][depth],2) + pow(esume[ieta][0][1][depth],2) + pow(esume[ieta][0][2][depth],2));
      // Make error bit larger for depth2,3
      if (depth==1 || depth==2) {
	ecore = esum[ieta][0][1][depth];
	err = esume[ieta][0][1][depth];
      }
      if (ecore!=0 && err>0 && val>0)
	h->SetBinError(j, val * (err / ecore) * ((1./val) / ecore));
    }
  }
  fout->Write();
  fout->Close();

  
  /////////////////////////////////////
  // Figure out the uncertainties    //
  /////////////////////////////////////
  /*
  // Re-access the reduced matrix Ared and vector bred
  const int Mred = Ared.GetNrows();
  
  // Compute the inverse of the reduced matrix
  TMatrixD AredInv = Ared;
  AredInv.Invert();
  
  // Initialize covariance matrix for b (based on input errors)
  TMatrixD Cov_b(Mred, Mred);
  Cov_b.Zero();
  for (int i = 0; i < Mred; ++i) {
    for (int j = 0; j < Mred; ++j) {
      double sumEprod = 0.0;
      for (int k = 0; k < nch; ++k) {
	for (int l = 0; l < nch; ++l) {
	  sumEprod += eprode[i][k][l] * Ared(i, k) * Ared(j, l);
	}
      }
      double sumEsume = 0.0;
      for (int k = 0; k < nch; ++k) {
	sumEsume += esume[i][k] * Ared(i, k);
      }
      Cov_b(i, j) = sumEprod + sumEsume;
    }
  }
  
  // Propagate uncertainties: Cov_X = AredInv * Cov_b * (AredInv)^T
  TMatrixD Cov_X = AredInv * Cov_b * AredInv.T();
  
  // Extract uncertainties: Delta_X = sqrt(diag(Cov_X))
  TVectorD Delta_X(Mred);
  for (int i = 0; i < Mred; ++i) {
    Delta_X[i] = std::sqrt(Cov_X(i, i));
  }
  
  // Map uncertainties back to full space
  TVectorD fullDeltaX(M);
  fullDeltaX.Zero();
  for (int irow = 0; irow < Mred; ++irow) {
    fullDeltaX[newToOld[irow]] = Delta_X[irow];
  }
  
    // Save uncertainties in histograms
  TFile *fout = new TFile("solveIsoTrackV2_withUncertainties.root", "RECREATE");
  for (int depth = 0; depth != ndepth; ++depth) {
    TH1D *h = new TH1D(Form("hf_dd_unc_%d", depth + 1), ";ieta;Uncertainty",
		       61, -30.5, 30.5);
    for (int i = 0; i != M; ++i) {
      if (i % ndepth != depth) continue;
      int ieta = i / ndepth;
      ieta -= 29; // Adjust ieta range
      if (ieta >= 0) ieta += 1;
      double unc = fullDeltaX[i];
      int j = h->FindBin(ieta);
      h->SetBinContent(j, unc);
    }
    h->Write();
  }
  */  
  //fout->Write();
  //fout->Close();
  
} // solveIsoTrackV2

