#include "SolveGrandMatrix.h"

void SolveGrandMatrix(char *fname, int nIndep){
    // open new tfile named grandFile and read in the matrices
    int nDets = readFile(fname, nIndep);
    if (!grandFile || grandFile->IsZombie()) {
        std::cout << "Could not open input file: " << fname << std::endl;
        return;
    }

    initMatices(nIndep, nDets - nIndep);
    solve(nIndep, nDets - nIndep);

    ofile = new TFile("~/Working/local_repos/japan-MOLLER/rootfiles/regression_output.root", "RECREATE", "regression output");
    WriteOutput(nIndep);
}


// returns number of independent + dependent variables
int readFile(char *fname, int nIndep){
    grandFile = new TFile(fname);
    if (!grandFile || grandFile->IsZombie()) {
        std::cout << "Could not open input file: " << fname << std::endl;
        return 0;
    }
    
    // i dont understand pointers
    TMatrixD _mNij = *grandFile->Get<TMatrixD>("N_ij");
    TMatrixD _mSij = *grandFile->Get<TMatrixD>("S_ij");
    TMatrixD _mMij = *grandFile->Get<TMatrixD>("M_ij");
    TMatrixD _mCij = *grandFile->Get<TMatrixD>("C_ij");
    TMatrixD _mVij = *grandFile->Get<TMatrixD>("V_ij");
    TMatrixD _mRij = *grandFile->Get<TMatrixD>("R_ij");
    TMatrixD _sigma_ij = *grandFile->Get<TMatrixD>("sigma_ij");
    TMatrixD _sigma_ji = *grandFile->Get<TMatrixD>("sigma_ji");

    mNij.ResizeTo(_mNij); mNij = _mNij;
    mSij.ResizeTo(_mSij); mSij = _mSij;
    mMij.ResizeTo(_mMij); mMij = _mMij;
    mCij.ResizeTo(_mCij); mCij = _mCij;
    mVij.ResizeTo(_mVij); mVij = _mVij;
    mRij.ResizeTo(_mRij); mRij = _mRij;
    sigma_ij.ResizeTo(_sigma_ij); sigma_ij = _sigma_ij;
    sigma_ji.ResizeTo(_sigma_ji); sigma_ji = _sigma_ji;

    // can probably assume that number of good events is maximum value in mNij
    fGoodEventNumber = mNij.Max();

    grandFile->GetObject("name_vector", allnames);

    int nDets = allnames->size();

    // we assume the independent variables are first for now

    grandFile->Close();

    return nDets;
}

void WriteOutput(int nIndep){

    // Write objects
    Axy.Write("slopes");
    dAxy.Write("sigSlopes");

    mRPP.Write("IV_IV_correlation");
    mRPY.Write("IV_DV_correlation");
    mRYY.Write("DV_DV_correlation");
    mRYYp.Write("DV_DV_correlation_prime");

    mMP.Write("IV_mean");
    mMY.Write("DV_mean");
    mMYp.Write("DV_mean_prime");


    // // number of events
    // TMatrixD Mstat(1,1);
    // Mstat(0,0)=getUsedEve();
    // Mstat.Write("MyStat");

    auto indepnames = new std::vector<std::string>(allnames->begin(), allnames->begin() + nIndep);
    auto depnames = new std::vector<std::string>(allnames->begin() + nIndep, allnames->end());

    ofile->WriteObjectAny(indepnames, "vector<string>", "IVnames");
    ofile->WriteObjectAny(depnames, "vector<string>", "DVnames");


    // sigmas
    mSP.Write("IV_sigma");
    mSY.Write("DV_sigma");
    mSYp.Write("DV_sigma_prime");



    // raw covariances
    mVPP.Write("IV_IV_rawVariance");
    mVPY.Write("IV_DV_rawVariance");
    mVYY.Write("DV_DV_rawVariance");
    mVYYp.Write("DV_DV_rawVariance_prime");
    TVectorD mVY2 = (TVectorD) TMatrixDDiag(mVYY);
    mVY2.Write("DV_rawVariance");
    TVectorD mVP2 = (TVectorD) TMatrixDDiag(mVPP);
    mVP2.Write("IV_rawVariance");
    TVectorD mVY2prime = (TVectorD) TMatrixDDiag(mVYYp);
    mVY2prime.Write("DV_rawVariance_prime");

    // normalized covariances
    mSPP.Write("IV_IV_normVariance");
    mSPY.Write("IV_DV_normVariance");
    mSYY.Write("DV_DV_normVariance");
    mSYYp.Write("DV_DV_normVariance_prime");
    TVectorD sigY2 = (TVectorD) TMatrixDDiag(mSYY);
    sigY2.Write("DV_normVariance");
    TVectorD sigX2 = (TVectorD) TMatrixDDiag(mSPP);
    sigX2.Write("IV_normVariance");
    TVectorD sigY2prime = (TVectorD) TMatrixDDiag(mSYYp);
    sigY2prime.Write("DV_normVariance_prime");


    Axy.Write("A_xy");
    Ayx.Write("A_yx");


    ofile->Write();
    ofile->Close();
}

void initMatices(int nP, int nY){
  mMP.ResizeTo(nP);
  mMY.ResizeTo(nY);
  mMYp.ResizeTo(nY);

  mVPP.ResizeTo(nP,nP);
  mVPY.ResizeTo(nP,nY);
  mVYP.ResizeTo(nY,nP);
  mVYY.ResizeTo(nY,nY);
  mVYYp.ResizeTo(nY,nY);
  mVP.ResizeTo(nP);
  mVY.ResizeTo(nY);
  mVYp.ResizeTo(nY);

  mSPP.ResizeTo(mVPP);
  mSPY.ResizeTo(mVPY);
  mSYP.ResizeTo(mVYP);
  mSYY.ResizeTo(mVYY);
  mSYYp.ResizeTo(mVYYp);

  Axy.ResizeTo(nP,nY);
  Ayx.ResizeTo(nY,nP);
  dAxy.ResizeTo(Axy);
  dAyx.ResizeTo(Ayx);

  mSP.ResizeTo(nP);
  mSY.ResizeTo(nY);
  mSYp.ResizeTo(nY);

  mRPP.ResizeTo(mVPP);
  mRPY.ResizeTo(mVPY);
  mRYP.ResizeTo(mVYP);
  mRYY.ResizeTo(mVYY);
  mRYYp.ResizeTo(mVYYp);

  mVFULL_clean.ResizeTo(nP+nY,nP+nY);
  mSFULL_clean.ResizeTo(nP+nY,nP+nY);
}

// np is independent, ny is dependent 
void solve(int nP, int nY){

    TMatrixD mVFULL_clean; mVFULL_clean.ResizeTo(mNij);
    TMatrixD mSFULL_clean; mSFULL_clean.ResizeTo(mNij);


    mVPY.ResizeTo(nP, nY); mVPY = mCij.GetSub(0,nP-1,nP,nP+nY-1);
    mVPP.ResizeTo(nP, nP); mVPP = mCij.GetSub(0,nP-1,0,nP-1);
    mVYY.ResizeTo(nY, nY); mVYY = mCij.GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mRPY.ResizeTo(nP, nY); mRPY = mRij.GetSub(0,nP-1,nP,nP+nY-1);
    mRPP.ResizeTo(nP, nP); mRPP = mRij.GetSub(0,nP-1,0,nP-1);
    mRYY.ResizeTo(nY, nY); mRYY = mRij.GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mSPY.ResizeTo(nP, nY); mSPY = mVij.GetSub(0,nP-1,nP,nP+nY-1);
    mSPP.ResizeTo(nP, nP); mSPP = mVij.GetSub(0,nP-1,0,nP-1);
    mSYY.ResizeTo(nY, nY); mSYY = mVij.GetSub(nP,nP+nY-1,nP,nP+nY-1);

    // off-diagonal raw covariance
    
    mVYP.Transpose(mVPY);
    // diagonal variances
    TMatrixD sigmaP = sigma_ij.GetSub(0,nP-1,0,nP-1);
    TMatrixD sigmaY = sigma_ij.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mVP = TMatrixDDiag(sigmaP);
    mVY = TMatrixDDiag(sigmaY);

    // "Clean" matrices
    for(int i = 0; i < allnames->size(); ++i){
        for(int j = i; j < allnames->size(); ++j){
            mVFULL_clean(i,j) = mRij(i,j) * sigma_ij(i,j) * sigma_ji(j,i) * (fGoodEventNumber - 1);
            mVFULL_clean(j,i) = mVFULL_clean(i,j);
            mSFULL_clean(i,j) = mRij(i,j) * sigma_ij(i,j) * sigma_ji(j,i);
            mSFULL_clean(j,i) = mSFULL_clean(i,j);
        }
    }

    TMatrixD mVPY_clean = mVFULL_clean.GetSub(0,nP-1,nP,nP+nY-1);
    TMatrixD mVPP_clean = mVFULL_clean.GetSub(0,nP-1,0,nP-1);
    TMatrixD mVYY_clean = mVFULL_clean.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    TMatrixD mVYP_clean(TMatrixD::kTransposed, mVPY_clean);



    TMatrixD mSPY_clean = mSFULL_clean.GetSub(0,nP-1,nP,nP+nY-1);
    TMatrixD mSPP_clean = mSFULL_clean.GetSub(0,nP-1,0,nP-1);
    TMatrixD mSYY_clean = mSFULL_clean.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    TMatrixD mSYP_clean(TMatrixD::kTransposed, mSPY_clean);

    TVectorD mSP_clean = TMatrixDDiag(mSPP_clean);
    mSP_clean.Sqrt();
    TVectorD mSY_clean = TMatrixDDiag(mSYY_clean);
    mSY_clean.Sqrt();



    // Check for goodness and then get rid of bad columns
    // Warn if correlation matrix determinant close to zero (heuristic)
    if (mRPP.Determinant() < std::pow(10,-(2*nP))) {
        std::cout << "LRB: correlation matrix nearly singular, "
                  << "determinant = " << mRPP.Determinant()
                  << " (set includes highly correlated variable pairs)"
                  << std::endl;
        if (fGoodEventNumber > 10*nP) {
            std::cout << fGoodEventNumber << " events" << std::endl;
            std::cout << "Covariance matrix: " << std::endl; mVPP_clean.Print();
            std::cout << "Correlation matrix: " << std::endl; mRPP.Print();
        }
        std::cout << "LRB: solving failed (this happens when only few events)."
                  << std::endl;
        return;
    }

    //==========================================================
    //Solve Step 3
    // slopes

    TMatrixD invRPP(TMatrixD::kInverted, mRPP);

    Axy = TMatrixD(invRPP, TMatrixD::kMult, mRPY);
    Axy.NormByColumn(mSP_clean); // divide
    Axy.NormByRow(mSY_clean, ""); // mult
    Ayx.Transpose(Axy);

    // new means
    mMYp = mMY - Ayx * mMP;

    // new raw covariance
    mVYYp = mVYY_clean + Ayx * mVPP_clean * Axy - (Ayx * mVPY_clean + mVYP_clean * Axy);

    // new variances
    mVYp = TMatrixDDiag(mVYYp); 

    for (int i = 0; i < mVYp.GetNrows(); i++) {
        if (mVYp(i) < 0 && fabs(mVYp(i)) < 1e-12) {
            mVYp(i) = 0;
        }
    }
    mVYp.Sqrt();

    // new normalized covariance
    mSYYp = mSYY_clean + Ayx * mSPP_clean * Axy - (Ayx * mSPY_clean + mSYP_clean * Axy);

    // uncertainties on the new means
    mSYp = TMatrixDDiag(mSYYp); 
    mSYp.Sqrt();

    // new correlation matrix
    mRYYp = mVYYp; 
    mRYYp.NormByColumn(mVYp); 
    mRYYp.NormByRow(mVYp);

    // slope uncertainties
    double norm = 1. / (fGoodEventNumber - nP - 1);
    dAxy.Zero();
    dAxy.Rank1Update(TMatrixDDiag(invRPP), TMatrixDDiag(mRYYp), norm); // diag mRYYp = row of ones
    dAxy.Sqrt();

    dAxy.NormByColumn(mSP_clean); // divide

    dAxy.NormByRow(mSYp, ""); // mult

    dAyx.Transpose(dAxy);
}