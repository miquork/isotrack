// Purpose: Construct and solve linear system of HcalRespCorr
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"

#include <iostream>

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
  //double esume[neta][nband][nwidth][ndepth];
  double eprod[neta][nch][nch];
  //double eprode[neta][nch][nch];
  assert(p2->GetNbinsX()==neta);
  assert(p2->GetNbinsY()==nch);
  assert(p3->GetNbinsX()==neta);
  assert(p3->GetNbinsY()==nch);
  assert(p3->GetNbinsZ()==nch);
  for (int ieta = 0; ieta != neta; ++ieta) {

    // Fast memory access by using fact that arrays are consecutive blocks
    double *pe = &(esum[ieta][0][0][0]); // pointer to first element per ieta
  //double *pee = &(esume[ieta][0][0][0]); // pointer to first element per ieta
    // int ch = (band,i)*nwidth*ndepth + (width,j)*ndepth + (depth,k);
    // double etot = *(p+ch); // = esum[i][j][k]
    for (int ch = 0; ch != nch; ++ch) {
      pe[ch] = p2->GetBinContent(ieta+1, ch+1);
      //pee[ch] = p2->GetBinError(ieta+1, ch+1);

      for (int ch2 = 0; ch2 != nch; ++ch2) {
	eprod[ieta][ch][ch2] = p3->GetBinContent(ieta+1, ch+1, ch2+1);
	//eprod[ieta][ch][ch2] = p3->GetBinError(ieta+1, ch+1, ch2+1);
      } // for ch2
    } // for ch
  } // for ieta

  /*
  // Remapping to ieta x (jeta x depth) matrix for response
  // ([p_i] =) [1] = [R_ij] * [c_j]
  //double R[neta][nc];
  TMatrixD R(neta,nc);
  TVectorD E(nc);

  // Reset matrix elements
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int ic = 0; ic != nc; ++ic) {
      R[ieta][ic] = 0;
    } // for depth
  } // for ieta

  // Reset vector elements
  for (int ic = 0; ic != nc; ++ic) {
    E[ic] = 0;
  } // for ieta
  
  // Distribute and resum means to corresponding c_i (ieta,depth) element
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int depth = 0; depth != ndepth; ++depth) {

      // Loop over delta_ieta and distribute to correct ieta
      for (int iw = 0; iw != nwidth; ++iw) {

	// Check that we don't outside bounds for ic
	if ((ieta+(iw-1))<0 || (ieta+(iw-1))>=neta) continue;
	
	int ic = ndepth*(ieta+(iw-1))+depth;
	// Reminder:   esum[neta][nband=2(core==0,side==1)][nwidth=3][ndepth];
	double E_tot = esum[ieta][0][iw][depth] - esum[ieta][1][iw][depth];
	R[ieta][ic] += E_tot; 
	E[ic] += E_tot;
	
      }
      
    } // for depth
  } // for ieta

  // Remapping to (ieta x depth) x (ieta x depth) matrix for covariance
  // sigma_corr^2 = c_T * V * c
  //double V[nc][nc];
  TMatrixD V(nc,nc);

  // Reset matrix elements
  for (int ic = 0; ic != nc; ++ic) {
    for (int ic2 = 0; ic2 != nc; ++ic2) {
      V[ic][ic2] = 0;
    } // for ic2
  } // for ic

  // Distribute and resum E_i*E_j products to c_i (ieta,depth) element
  // Just one loop over ieta, because products non-zero only for ieta==jeta
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int idep = 0; idep != ndepth; ++idep) {
      for (int jdep = 0; jdep != ndepth; ++jdep) {

	// We now have so many width pairs, that make loops explicit
	for (int iw = 0; iw != nwidth; ++iw) {
	  for (int jw = 0; jw != nwidth; ++jw) {

	    // Check that we don't go outside bounds for ic, jc
	    if ((ieta+(iw-1))<0 || (ieta+(iw-1))>=neta ||
		(ieta+(jw-1))<0 || (ieta+(jw-1))>=neta) continue;

	    // Map ieta, depth, width to correct c_i
	    int ic = ndepth*(ieta+(iw-1))+idep;
	    int jc = ndepth*(ieta+(jw-1))+jdep;
	    // Reminder: ch = (band)*nwidth*ndepth + (width)*ndepth + (depth);
	    int ich0 = 0*nband*ndepth + iw*ndepth + idep;
	    int jch0 = 0*nband*ndepth + jw*ndepth + jdep;
	    int ich1 = 1*nband*ndepth + iw*ndepth + idep;
	    int jch1 = 1*nband*ndepth + jw*ndepth + idep;
	    // E_eff,i*E_eff,j = (E_i0-E_i1)*(E_j0-E_j1)
	    //                 = E_i0*E_j0-E_i0*E_j1-E_i1*E_j0+E_i1*E_j1
	    V[ic][jc] +=
	      eprod[ieta][ich0][jch0] - eprod[ieta][ich0][jch1] -
	      eprod[ieta][ich1][jch0] + eprod[ieta][ich1][jch1];
	
	  } // for jw
	} // for iw
	  
      } // for jdep
    } // for idep
  } // for ieta

  // Calculate covariance as E_ij-E_i*Ej
  for (int ic = 0; ic != nc; ++ic) {
    for (int jc = 0; jc != nc; ++jc) {
      V[ic][jc] = V[ic][jc] - E[ic]*E[jc];
    } // for jc
  } // for ic

  // Also store normalized matrix for easier checking
  TMatrixD V_norm(nc,nc);
  for (int ic = 0; ic != nc; ++ic) {
    for (int jc = 0; jc != nc; ++jc) {
      if (V[ic][ic]==0 || V[jc][jc]==0)
	V_norm[ic][jc] = 0;
      else
	V_norm[ic][jc] = V[ic][jc]/sqrt(V[ic][ic]*V[jc][jc]);
    } // for jc
  } // for ic
  
  // Draw matrix as a sanity check
  //R.Draw("COLZ");
  //V.Draw("COLZ");
  V_norm.Draw("COLZ");
  */

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

		int ibandAlpha = alpha/(nwidth*ndepth);
		int ibandBeta = beta/(nwidth*ndepth);
		
		double val = eprod[i][alpha][beta]; 
		// beta channel corresponds to some (iband2, iwidth2, idepth2)
		int iband2 = beta/(nwidth*ndepth);
		int rest = beta%(nwidth*ndepth);
		int iwidth2 = rest/ndepth;
		int idepth2 = rest%ndepth;
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
  // Solve it:
  TDecompLU solver(A);
  bool ok = solver.Solve(b);
  if (!ok) {
    cout << "Solve failed!" << endl;
  }
  
  // The solution is in b now:
  double corr[nc];
  for (int ic=0; ic<nc; ic++) {
    corr[ic] = b[ic]; 
  }

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
  int color[ndepth] = {kBlue, kOrange+1, kGreen+1, kRed,
		       kYellow+1, kOrange-7, kGray+2};
  for (int depth = 0; depth != ndepth; ++depth) {
    TH1D *h = new TH1D(Form("hf_dd_%d",depth+1),";ieta;CorrFactV2",
		       61,-30.5,30.5);
    h->SetMarkerColor(color[depth]);
    h->SetLineColor(color[depth]);
    for (int i = 0; i != M; ++i) {
      // int i = ieta*ndepth+depth
      if (i%ndepth!=depth) continue;
      int ieta = i/ndepth; // 0-57
      ieta -= 29; // -29-28
      if (ieta>=0) ieta += 1; // -29 to +29 with 0 excluded
      double val = fullX[i];
      int j = h->FindBin(ieta);
      h->SetBinContent(j, val);
    }
  }
  fout->Write();
  fout->Close();

  
  /*
  // Regions for (core,side)x(delta_ieta=-1,0,+1)x(depth)
  for (int i = 0; i != nband; ++i) { // (core,side)
    for (int j = 0; j != nwidth; ++j) {
      for (int k = 0; k != ndepth; ++k) {
	esum[i][j][k] = 0;
      }
    } // for j
  } // for i
  
  // p(i) = Sum_d,di c_ieta,depth (Ecore_ieta,depth - k*Eside_ieta,depth)
  
  // Take out empty bins and figure out mappings
  int bin(0);
  map<int, int> mp2toM, mMtop2;
  for (int i = 1; i != p2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != p2->GetNbinsY()+1; ++i) {
      if (p2->GetBinContent(i,j)!=0 && p2->GetBinError(i,j)!=0) {
	int k = p2->GetBin(i,j);
	mp2toM[k] = bin;
	mMtop2[bin] = k;
	++bin;
      }
    } // for j
  } // for i
  const int N_tot = bin;

  cout << "Found "<<N_tot<<" non-empty elements\n";

  TVectorD E_eff_total(N_tot);
  TMatrixD V_total(N_tot, N_tot);

  for (int i = 0; i != N_tot; ++i) {
    E_eff_total
  }
  */  
} // solveIsoTrackV2
