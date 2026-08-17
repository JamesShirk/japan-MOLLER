#include <fstream>
#include <sstream>
#include "Correlator.h"

// Makes temp matricies that we will change so that we don't interfere with the original values of the ROOT file
void MakeTemp(std::string name) {
    std::string tempName = name + "_temp";
    if (!(gDirectory->Get(tempName.c_str()) != nullptr)) {
        TMatrixT<double>* matrix = (TMatrixT<double>*)gDirectory->Get(name.c_str());
        matrix->Write(tempName.c_str(), TObject::kOverwrite);
        std::cout << "Made the temp and backup " << name << std::endl;
    }
}

// Makes the temp matricies needed for the solve function
void initialize() {
    MakeTemp("C_ij");
    MakeTemp("R_ij");
    MakeTemp("V_ij");
    MakeTemp("N_ij");
    MakeTemp("S_ij");
    MakeTemp("M_ij");
    MakeTemp("sigma_ij");
}

// This switches one row and column with another for a defined matrix
void SwitchM(TMatrixT<double> &matrix, int &rc1, int &rc2) {
    for (int col = 0; col < matrix.GetNcols(); col++) {
        double temp = matrix(rc1, col);
        matrix(rc1, col) = matrix(rc2, col);
        matrix(rc2, col) = temp;
    }
    for (int row = 0; row < matrix.GetNrows(); row++) {
        double temp = matrix(row, rc1);
        matrix(row, rc1) = matrix(row, rc2);
        matrix(row, rc2) = temp;
    }
}

// This switches to elements inside a vector
template <typename type>
void SwitchV(type &vector, int &elem1, int &elem2) {
    auto temp = vector[elem2];
    vector[elem2] = vector[elem1];
    vector[elem1] = temp;
}

// This grabs the temp matricies and then does the defined columns/row switching then saves the temp matricies back to the root file
// This can def be more efficient it currently opens, changes, then saves the file instead it should open do all the operations then save.
void SwitchRowCol (int rc1, int rc2) {

    TMatrixT<double>* mCij = (TMatrixT<double>*)gDirectory->Get("C_ij_temp");
    TMatrixT<double>* mRij = (TMatrixT<double>*)gDirectory->Get("R_ij_temp");
    TMatrixT<double>* mVij = (TMatrixT<double>*)gDirectory->Get("V_ij_temp");
    TMatrixT<double>* mNij = (TMatrixT<double>*)gDirectory->Get("N_ij_temp");
    TMatrixT<double>* mSij = (TMatrixT<double>*)gDirectory->Get("S_ij_temp");
    TMatrixT<double>* mMij = (TMatrixT<double>*)gDirectory->Get("M_ij_temp");
    TMatrixT<double>* sigma_ij = (TMatrixT<double>*)gDirectory->Get("sigma_ij_temp");

    SwitchM(*mCij, rc1, rc2);
    SwitchM(*mRij, rc1, rc2);
    SwitchM(*mVij, rc1, rc2);
    SwitchM(*mNij, rc1, rc2);
    SwitchM(*mSij, rc1, rc2);
    SwitchM(*mMij, rc1, rc2);
    SwitchM(*sigma_ij, rc1, rc2);

    mCij->Write("C_ij_temp", TObject::kOverwrite);
    mRij->Write("R_ij_temp", TObject::kOverwrite);
    mVij->Write("V_ij_temp", TObject::kOverwrite);
    mNij->Write("N_ij_temp", TObject::kOverwrite);
    mSij->Write("S_ij_temp", TObject::kOverwrite);
    mMij->Write("M_ij_temp", TObject::kOverwrite);
    sigma_ij->Write("sigma_ij_temp", TObject::kOverwrite);
}

// This switches 2 elements within the vectors that we need to swap with the matricies
void SwitchElem(int elem1, int elem2) {
    TVectorT<double>* mean = (TVectorT<double>*)gDirectory->Get("mean_temp");

    SwitchV(*mean, elem1, elem2);
    SwitchV(CombinedName, elem1, elem2);
    SwitchV(CombinedType, elem1, elem2);

    mean->Write("mean_temp", TObject::kOverwrite);
}

// This makes sure that we are not just switching two variables push pushing up or down a variable to maintain the original order
void insertRowCol(int rc1, int rc2, int type) {
    if (rc1 > rc2) {
        int temp = rc1;
        rc1 = rc2;
        rc2 = temp;
    }
// For moving something up (changing to iv)
if (type == 1) {
    for (int i = rc2; i > rc1; i--) {
        SwitchRowCol(i-1, i);
        SwitchElem(i-1, i);
        std::cout << i-1 << " " << i << std::endl;
    } }
// For moving something down (like removing)
if (type == 2) {
    for (int i = rc1; i < rc2; i++) {
        SwitchRowCol(i, i+1);
        SwitchElem(i, i+1);
        std::cout << i << " " << i+1 << std::endl;
    } }
}

// This clears the temp matricies
void Delete() {
    gDirectory->Delete("C_ij_temp;*");
    gDirectory->Delete("R_ij_temp;*");
    gDirectory->Delete("V_ij_temp;*");
    gDirectory->Delete("N_ij_temp;*");
    gDirectory->Delete("S_ij_temp;*");
    gDirectory->Delete("M_ij_temp;*");
    gDirectory->Delete("sigma_ij_temp;*");
    gDirectory->Delete("DV_mean_temp;*");
    gDirectory->Delete("IV_mean_temp;*");
    gDirectory->Delete("mean_temp;*");
    std::cout << "Matricies cleared " << std::endl;
}

// Initializes the matrix sizes
void resize(int nP, int nY)
{
    mMPP.ResizeTo(nP,nP);         mMYY.ResizeTo(nY,nY);
    mVPY.ResizeTo(nP,nY);         mVPP.ResizeTo(nP,nP);         mVYY.ResizeTo(nY,nY);
    mRPY.ResizeTo(nP,nY);         mRPP.ResizeTo(nP,nP);         mRYY.ResizeTo(nY,nY);
    mSPY.ResizeTo(nP,nY);         mSPP.ResizeTo(nP,nP);         mSYY.ResizeTo(nY,nY);
    sigmaP.ResizeTo(nP,nP);       sigmaY.ResizeTo(nY,nY);
    mVPY_clean.ResizeTo(nP,nY);   mVPP_clean.ResizeTo(nP,nP);   mVYY_clean.ResizeTo(nY,nY); mVYP_clean.ResizeTo(nY,nP);
    mSPY_clean.ResizeTo(nP,nY);   mSPP_clean.ResizeTo(nP,nP);   mSYY_clean.ResizeTo(nY,nY); mSYP_clean.ResizeTo(nY,nP);
    Axy.ResizeTo(nP,nY);          Ayx.ResizeTo(nY,nP);
    dAxy.ResizeTo(nP,nY);         dAyx.ResizeTo(nY,nP);
    mVYP.ResizeTo(nY,nP);         mSYP.ResizeTo(nY,nP);

    mMP.ResizeTo(nP);        mMY.ResizeTo(nY);
    mVP.ResizeTo(nP);        mVY.ResizeTo(nY);
    mSP_clean.ResizeTo(nP);  mSY_clean.ResizeTo(nY);
    mVYp.ResizeTo(nY);       mSYp.ResizeTo(nY);
    mMYp.ResizeTo(nY);
}

// The original solve function with some minor tweaks to make sure the original matricies are used with the row/column matricies
// This also saves the solved matricies into a new ROOT file which has an identicaly format the the original file
void solve(int nP, int nY){

    resize(nP, nY);

    TMatrixT<double>* mCij = (TMatrixT<double>*)gDirectory->Get("C_ij_temp");
    TMatrixT<double>* mRij = (TMatrixT<double>*)gDirectory->Get("R_ij_temp");
    TMatrixT<double>* mVij = (TMatrixT<double>*)gDirectory->Get("V_ij_temp");
    TMatrixT<double>* mNij = (TMatrixT<double>*)gDirectory->Get("N_ij_temp");
    TMatrixT<double>* mSij = (TMatrixT<double>*)gDirectory->Get("S_ij_temp");
    TMatrixT<double>* mMij = (TMatrixT<double>*)gDirectory->Get("M_ij_temp");
    TMatrixT<double>* sigma_ij = (TMatrixT<double>*)gDirectory->Get("sigma_ij_temp");

    if (!mCij||!mRij||!mVij||!mNij||!mSij||!mMij||!sigma_ij) {
        std::cout << "ERROR: one or more _temp matrices missing from file!" << std::endl;
        return;
    }

    TMatrixD Cij = mCij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Rij = mRij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Vij = mVij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Nij = mNij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Sij = mSij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Mij = mMij->GetSub(0,nP+nY-1,0,nP+nY-1);
    TMatrixD Sigmaij = sigma_ij->GetSub(0,nP+nY-1,0,nP+nY-1);

    mMPP = mMij->GetSub(0,nP-1,0,nP-1);
    mMYY = mMij->GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mMP = TMatrixDDiag(mMPP);
    mMY = TMatrixDDiag(mMYY);

    mVFULL_clean; mVFULL_clean.ResizeTo(*mNij);
    mSFULL_clean; mSFULL_clean.ResizeTo(*mNij);

    // still just assume the independent variables are first
    mVPY = mCij->GetSub(0,nP-1,nP,nP+nY-1);
    mVPP = mCij->GetSub(0,nP-1,0,nP-1);
    mVYY = mCij->GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mRPY = mRij->GetSub(0,nP-1,nP,nP+nY-1);
    mRPP = mRij->GetSub(0,nP-1,0,nP-1);
    mRYY = mRij->GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mSPY = mVij->GetSub(0,nP-1,nP,nP+nY-1);
    mSPP = mVij->GetSub(0,nP-1,0,nP-1);
    mSYY = mVij->GetSub(nP,nP+nY-1,nP,nP+nY-1);

    // off-diagonal raw covariance


    fGoodEventNumber = mNij->Max();
    mVYP.Transpose(mVPY);
    // diagonal variances
    sigmaP = sigma_ij->GetSub(0,nP-1,0,nP-1);
    sigmaY = sigma_ij->GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mVP = TMatrixDDiag(sigmaP);
    mVY = TMatrixDDiag(sigmaY);
    // "Clean" matrices
    for(int i = 0; i <CombinedName.size(); ++i){
        for(int j = i; j <CombinedName.size(); ++j){
            mVFULL_clean(i,j) = (*mRij)(i,j) * (*sigma_ij)(i,j) * (*sigma_ij)(j,i) * (fGoodEventNumber - 1);
            mVFULL_clean(j,i) = mVFULL_clean(i,j);
            mSFULL_clean(i,j) = (*mRij)(i,j) * (*sigma_ij)(i,j) * (*sigma_ij)(j,i);
            mSFULL_clean(j,i) = mSFULL_clean(i,j);
        }
    }

    mVPY_clean = mVFULL_clean.GetSub(0,nP-1,nP,nP+nY-1);
    mVPP_clean = mVFULL_clean.GetSub(0,nP-1,0,nP-1);
    mVYY_clean = mVFULL_clean.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mVYP_clean.Transpose(mVPY_clean);

    mSPY_clean = mSFULL_clean.GetSub(0,nP-1,nP,nP+nY-1);
    mSPP_clean = mSFULL_clean.GetSub(0,nP-1,0,nP-1);
    mSYY_clean = mSFULL_clean.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mSYP_clean.Transpose(mSPY_clean);

    mSP_clean = TMatrixDDiag(mSPP_clean);
    mSP_clean.Sqrt();
    mSY_clean = TMatrixDDiag(mSYY_clean);
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
    TMatrixD mVYYp = mVYY_clean + Ayx * mVPP_clean * Axy - (Ayx * mVPY_clean + mVYP_clean * Axy);

    // new variances
    mVYp = TMatrixDDiag(mVYYp); 

    for (int i = 0; i < mVYp.GetNrows(); i++) {
        if (mVYp(i) < 0 && fabs(mVYp(i)) < 1e-12) {
            mVYp(i) = 0;
        }
        if (mVYp(i) < 0) {
        std::cout << "mVYp(" << i << ") = " << mVYp(i) << std::endl;
        }
    }

    mVYp.Sqrt();

    // new normalized covariance
    TMatrixD mSYYp = mSYY_clean + Ayx * mSPP_clean * Axy - (Ayx * mSPY_clean + mSYP_clean * Axy);

    // uncertainties on the new means
    mSYp = TMatrixDDiag(mSYYp); 
    mSYp.Sqrt();

    // new correlation matrix
    TMatrixD mRYYp = mVYYp; 
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

//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------

  TFile* f = new TFile("newblueR_grand.slope.root", "RECREATE");
  Axy.Write("slopes");
  dAxy.Write("sigSlopes");

  mRPP.Write("IV_IV_correlation");
  mRPY.Write("IV_DV_correlation");
  mRYY.Write("DV_DV_correlation");
  mRYYp.Write("DV_DV_correlation_prime");

  mMP.Write("IV_mean");
  mMY.Write("DV_mean");
  mMYp.Write("DV_mean_prime");

  // number of events
//   TMatrixD Mstat(1,1);
//   Mstat(0,0)=getUsedEve();
//   Mstat.Write("MyStat");

  // sigmas
  mSP.Write("IV_sigma");
  mSY.Write("DV_sigma");
  mSYp.Write("DV_sigma_prime");

  // raw covariances
  mVPP.Write("IV_IV_rawVariance");
  mVPY.Write("IV_DV_rawVariance");
  mVYY.Write("DV_DV_rawVariance");
  mVYYp.Write("DV_DV_rawVariance_prime");
  TVectorD mVY2((TMatrixDDiag(mVYY)));
  mVY2.Write("DV_rawVariance");
  TVectorD mVP2((TMatrixDDiag(mVPP)));
  mVP2.Write("IV_rawVariance");
  TVectorD mVY2prime((TMatrixDDiag(mVYYp)));
  mVY2prime.Write("DV_rawVariance_prime");

  // normalized covariances
  mSPP.Write("IV_IV_normVariance");
  mSPY.Write("IV_DV_normVariance");
  mSYY.Write("DV_DV_normVariance");
  mSYYp.Write("DV_DV_normVariance_prime");
  TVectorD sigY2((TMatrixDDiag(mSYY)));
  sigY2.Write("DV_normVariance");
  TVectorD sigX2((TMatrixDDiag(mSPP)));
  sigX2.Write("IV_normVariance");
  TVectorD sigY2prime((TMatrixDDiag(mSYYp)));
  sigY2prime.Write("DV_normVariance_prime");
 
  Axy.Write("A_xy");
  Ayx.Write("A_yx");

  Nij.Write("N_ij");
  Sij.Write("S_ij");
  Mij.Write("M_ij");
  Cij.Write("C_ij");
  Vij.Write("V_ij");
  Rij.Write("R_ij");
  Sigmaij.Write("sigma_ij");

  //... IVs
  TH1D hiv("IVname","names of IVs",nP,-0.5,nP-0.5);
  for (int i=0;i<nP;i++) hiv.Fill(CombinedName[i].c_str(),i);
  hiv.Write();

  //... DVs
  TH1D hdv("DVname","names of IVs",nY,-0.5,nY-0.5);
  for (int i=nP;i<nP+nY;i++) hdv.Fill(CombinedName[i].c_str(),i);
  hdv.Write();

  f->Close();
}

// This is run to confirm the temp matricies and run the solve function and then create a new ROOT file
void Commit() {
    for (int i =0; i < Operations.size(); i++) {
        std::cout << Operations[i] << std::endl;
    }
    std::string response;
    std::cout << "Does this sound right? (y/n):" << std::endl;
    std::cin >> response;
    if (response == "y") { 
    solve(nP, nY);
    std::cout << "The new file has been made" << std::endl;
    Delete();
    TFile *_file0 = TFile::Open("blueR_grand4new.slope.root", "UPDATE");
} else {
    std::cout << "Operation Canceled" << std::endl;
} }

// what you actually run to change the rows and columns of the temp matricies
void ChangeIVDV () {

    CombinedName.clear();
    CombinedType.clear();
    InputName.clear();
    InputType.clear();
    IndepVar.clear();
    RemoveVar.clear();
    Operations.clear();
    Delete();
    initialize();

// grabs the raw names and categorizes them
    std::vector<std::string> *depNames = (std::vector<std::string>*)gDirectory->Get("Dependent_Names");
    std::vector<std::string> *indepNames = (std::vector<std::string>*)gDirectory->Get("Independent_Names");
    for (int i = 0; i < indepNames->size(); i++) {
        CombinedName.push_back((*indepNames)[i]); 
        CombinedType.push_back("iv"); }
    for (int i = 0; i < depNames->size(); i++) {
        CombinedName.push_back((*depNames)[i]);
        CombinedType.push_back("dv"); }

// grabs the input file and makes two lists for type of variables ("iv" or "dv") and the name for each variable
    std::ifstream file("rootScripts/Correlator.txt");
    if (!file.is_open()) {
        std::cout << "FAILED TO OPEN FILE!" << std::endl;
        return;
    }
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> type >> name;
        name.erase(0, 5);
        InputName.push_back(name);
        InputType.push_back(type); 
    }

    //Making sure the the variable actually existed in the original matrix and that they have the correct type
    // for (int i = 0; i < InputName.size(); i++) {
    //     bool locate = false;
    //     for (int j = 0; j < CombinedName.size(); j++) {
    //         if (InputName[i] == CombinedName[j]) {
    //             locate = true;
    //         }
    //     }
    //     if (locate == false) {
    //         std::cout << InputType[i] << "/" << InputName[i] << " is not listed in the original matrix. Aborting program" << std::endl;
    //         return;
    //     }
    // }

    // This checks the input file and tells us what we are removing and changing in the original matrix 
    //also it edits the CombinedType so we can work independent of the Input file
    for (int i =0; i < CombinedName.size(); i++) {
        bool found = false;
        for (int j = 0; j < InputName.size(); j++) {
            if (CombinedName[i] == InputName[j]) {
                found = true;
                if (CombinedType[i] == InputType[j]) {
                    continue;
                }
                else {
                    std::ostringstream oss;
                    oss << "Changing " << CombinedName[i] << " from " << CombinedType[i] << " to " << InputType[j];
                    std::cout << oss.str() << std::endl;
                    Operations.push_back(oss.str());
                    CombinedType[i] = InputType[j];
                }
                if (InputType[j] == "rm") {
                    found = false;
                }
            }
        }
        if (!found) {
            std::cout << "Removing " << CombinedName[i] << " from original matrix" << std::endl;
            RemoveVar.push_back(i);
        }
    }

    //This removes variables that do not appear in the input file but do appear in the original matrix or have rm as the type
    //It doesn't actually remove anything but moves the unwanted variables on the bottom so we can filter them out in the actual solver code
    //Note that we can remove from the original matrix but we can't add new variables
    int t = CombinedName.size() - 1;
    for (int i = RemoveVar.size() - 1; i >= 0; i--) {
        insertRowCol(t, RemoveVar[i], 2);
        t -= 1;
    }

    //IndepVar lists the # of independent variables so that the sorting loop knows how many row swaps to do
    for (int i = 0; i < InputType.size(); i++) {
        if (CombinedType[i] == "iv") {
            IndepVar.push_back(CombinedName[i]);
        } }
    
    //defines # of dependent and independent variables I hope that this can change nY and nP in the actual solve funtion to define the submatrix cuts
    nP = IndepVar.size();
    nY = InputType.size() - IndepVar.size() - RemoveVar.size();
    nR = RemoveVar.size();

    //Simple output to make sure that the # of dependent and independent variables are what we expect
    std::cout << "Original indep and dep sizes are " << indepNames->size() << " " << depNames->size() << std::endl;
    std::cout << "New indep and dep sizes are " << nP << " " << nY << std::endl;
    std::cout << "We are removing " << nR << " variables" << std::endl;

    // A simple sorting function that puts any variable marked with iv in the first nP rows and columns
    // It looks for the iv with the smallest position in the vector (example 0 would be the lowest) and sorts it into the 0th row/columns
    // It does this for every iv but won't do row operations if it finds that an iv is already in an early row
    // it then records the row/column operations done so it can undo them when someone plugs in a new set of iv and dv
    j = 0;
    for (int i = 0; i < IndepVar.size(); i++) {
        int earlyI = -1;
        for (int h = j; h < InputType.size(); h++) {
            if (CombinedType[h] == "iv") {
                earlyI = h;
                break;
            }
        }
        rc1 = j;
        rc2 = earlyI;
        j += 1;
        if (rc1 == rc2) {
            continue;
        }
        insertRowCol(rc1, rc2, 1);
    }
}