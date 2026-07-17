#include "SolveGrandMatrix.h"

void SolveGrandMatrix(char *fname, char *configFile = "/Users/jamesshirk/Working/local_repos/japan-MOLLER/rootScripts/correlator/grandcorrelator.conf") {

    // open new tfile named grandFile and read in the matrices
    int nDets = readFile(fname);

    // check if the file actually opened
    if (!grandFile || grandFile->IsZombie()) {
        std::cout << "Could not open input file: " << fname << std::endl;
        return;
    }

    // get number of independent variables from the config file and move them to the front of the matrix (if they arent there already)
    int nIndep = getDVIV(configFile);

    for(auto name : allnames){
        std::cout << name << std::endl;
    }

    // cant initialize matrices until we know how many dets, which we get from size of the grand matrices
    initMatices(nIndep, nDets - nIndep);
    solve(nIndep, nDets - nIndep);

    ofile = new TFile("~/Working/local_repos/japan-MOLLER/rootfiles/regression_output.root", "RECREATE", "regression output");
    WriteOutput(nIndep);

}

// returns number of independent + dependent variables
int readFile(char *fname){
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

    // we dont know how big the matrices are until we read them
    // moot if we declare the variable here too but then we couldnt use them in sovle()
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

    std::vector<std::string> *_allnames;
    grandFile->GetObject("names", _allnames);

    allnames = (*_allnames);
    // allnames = new std::vector<std::string>(_allnames->begin(), _allnames->end());

    int nDets = allnames.size();

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

    auto indepnames = new std::vector<std::string>(allnames.begin(), allnames.begin() + nIndep);
    auto depnames = new std::vector<std::string>(allnames.begin() + nIndep, allnames.end());

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

// tmatrices and tvectors default initialize to 0x0 (or length 0 given tvectors)
// if we want them to simultanously be global and be able to assign them later we have to initialize them to the right size
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

// returns number of IVs
// need to add checks to see if variable exists
int getDVIV(char *fname){
    
    std::ifstream infile(fname);
    std::string line;
    std::vector<std::string> indepvars;

    while (std::getline(infile, line)) {
        // ignore comments
        if (line.find("#") != std::string::npos) {
            continue;
        }
        else if (line.find("iv") != std::string::npos) {
            indepvars.push_back(line.erase(0, line.find("iv ") + 3));
        }
    }

    std::vector<int> indep_indices;

    for (const auto& name : indepvars) {
        auto pos = find(allnames.begin(), allnames.end(), name) - allnames.begin();
        if ((int) pos <allnames.size()) indep_indices.push_back(pos);
    }

    // return number of ivs and also move them to the front of the matrix if they arent there already
    return insertIV(indep_indices);
}

// move any IVs that arent at the front and then return number of IVs
int insertIV(std::vector<int> indep_indices) {
    int nIV = indep_indices.size();
    int infront = 0;
    // we dont need to move variables to the front if they're already there
    for (int idx : indep_indices) {
        if (idx > infront && idx < nIV) {
            infront = idx;
        }
    }

    infront++;

    for (int idx : indep_indices) {
        if (idx > infront) {
            std::cout << "Moving variable " << allnames[idx] << " to " << infront << std::endl;
            // much slower but this will preserve the order of the variables
            insertRowCol(infront, idx);
            //SwitchMatrices(idx, infront);
            std::cout << "New variable at " << idx << " is " << allnames[idx] << std::endl;
            infront++;
        }
    }
    return nIV;
}

// in grand matrix, this swaps the row and column for one variable with another
// we also swap positions in the allnames vector so that we can keep track of which variable is which
void SwitchRowCol (int rc1, int rc2, TMatrixD *A) {

    const TVectorD col_1 = TMatrixDColumn_const((*A),rc1);
    const TVectorD row_1 = TMatrixDRow_const((*A),rc1);

    TMatrixDRow((*A),rc1) = TMatrixDRow_const((*A),rc2);
    TMatrixDRow((*A),rc2) = row_1;

    TMatrixDColumn((*A),rc1) = TMatrixDColumn_const((*A),rc2);
    TMatrixDColumn((*A),rc2) = col_1;
    
}

// void SwitchRowCol (int rc1, int rc2, TMatrixD *A) {
//     for (int col = 0; col < A->GetNcols(); col++) {
//         double temp = (*A)(rc1, col);
//         (*A)(rc1, col) = (*A)(rc2, col);
//         (*A)(rc2, col) = temp;
//     }
//     for (int row = 0; row < A->GetNrows(); row++) {
//         double temp = (*A)(row, rc1);
//         (*A)(row, rc1) = (*A)(row, rc2);
//         (*A)(row, rc2) = temp;
//     }
// }

// much slower but this will preserve the order of the variables
void insertRowCol(int rc1, int rc2) {
    if (rc1 > rc2) {
        int temp = rc1;
        rc1 = rc2;
        rc2 = rc1;
    }
    for (int i = rc2; i > rc1; i--) {
        SwitchMatrices(i - 1, i);
    }
}

void SwitchMatrices(int rc1, int rc2) {
    SwitchRowCol(rc1, rc2, &mCij);
    SwitchRowCol(rc1, rc2, &mRij);
    SwitchRowCol(rc1, rc2, &mVij);
    SwitchRowCol(rc1, rc2, &mNij);
    SwitchRowCol(rc1, rc2, &mMij);
    SwitchRowCol(rc1, rc2, &mSij);
    SwitchRowCol(rc1, rc2, &sigma_ij);
    SwitchRowCol(rc1, rc2, &sigma_ji);
    std::iter_swap(allnames.begin() + rc1, allnames.begin() + rc2);
}

// np is independent, ny is dependent 
void solve(int nP, int nY){

    TMatrixD mVFULL_clean; mVFULL_clean.ResizeTo(mNij);
    TMatrixD mSFULL_clean; mSFULL_clean.ResizeTo(mNij);


    // still just assume the independent variables are first
    mVPY = mCij.GetSub(0,nP-1,nP,nP+nY-1);
    mVPP = mCij.GetSub(0,nP-1,0,nP-1);
    mVYY = mCij.GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mRPY = mRij.GetSub(0,nP-1,nP,nP+nY-1);
    mRPP = mRij.GetSub(0,nP-1,0,nP-1);
    mRYY = mRij.GetSub(nP,nP+nY-1,nP,nP+nY-1);


    mSPY = mVij.GetSub(0,nP-1,nP,nP+nY-1);
    mSPP = mVij.GetSub(0,nP-1,0,nP-1);
    mSYY = mVij.GetSub(nP,nP+nY-1,nP,nP+nY-1);

    // off-diagonal raw covariance
    
    mVYP.Transpose(mVPY);
    // diagonal variances
    TMatrixD sigmaP = sigma_ij.GetSub(0,nP-1,0,nP-1);
    TMatrixD sigmaY = sigma_ij.GetSub(nP,nP+nY-1,nP,nP+nY-1);
    mVP = TMatrixDDiag(sigmaP);
    mVY = TMatrixDDiag(sigmaY);

    // "Clean" matrices
    for(int i = 0; i <allnames.size(); ++i){
        for(int j = i; j <allnames.size(); ++j){
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

// // vector I define to help "remember" the actions taken in one ChangeIVDV run so that it can automatically reset when you run the program again
// // if a better method of this exists please let me know
// std::vector<int> inst;

// // what you actually run to change the rows and columns
// // the input file needs to have the same format as the matricies with iv and dv in from and you simply change some from iv to dv or vice versa to make changes to the variables
// void ChangeIVDV () {

// // grabs the raw names and categorizes them
//     std::vector<std::string> *depNames = (std::vector<std::string>*)gDirectory->Get("Dependent_Names");
//     std::vector<std::string> *indepNames = (std::vector<std::string>*)gDirectory->Get("Independent_Names");
//     std::vector<std::string> CombinedName;
//     std::vector<std::string> CombinedType;
//     for (int i = 0; i < indepNames->size(); i++) {
//         CombinedName.push_back((*indepNames)[i]); 
//         CombinedType.push_back("iv"); }
//     for (int i = 0; i < depNames->size(); i++) {
//         CombinedName.push_back((*depNames)[i]);
//         CombinedType.push_back("dv"); }

// // grabs the input file and makes two lists for type of variables ("iv" or "dv") and the name for each variable
//     std::ifstream file("rootScripts/Correlator.txt");
//     std::string line;
//     std::vector<std::string> InputType, InputName, NP;
//     if (!file.is_open()) {
//         std::cout << "FAILED TO OPEN FILE!" << std::endl;
//         return;
//     }
//     int row1, row2;
//     while (std::getline(file, line)) {
//         std::istringstream iss(line);
//         std::string type, name;
//         iss >> type >> name;
//         name.erase(0, 5);
//         InputName.push_back(name);
//         InputType.push_back(type); 
//     }

//     // This undoes the last row and column swaps done so you don't have to run the program twice to get the original matrix back
//     if (inst.size() >= 1) {
//         for (int i = 0; i < inst.size(); i += 2) {
//             int cr1 = inst[inst.size() - i - 1]; int cr2 = inst[inst.size() - i - 2];
//             SwitchRowCol(cr1, cr2);
//         }
//     }
//     inst.clear();

//     //Making sure the the variable actually existed in the original matrix and that they have the correct type
//     for (int i = 0; i < InputName.size(); i++) {
//         bool locate = false;
//         for (int j = 0; j < CombinedName.size(); j++) {
//             if (InputName[i] == CombinedName[j]) {
//                 locate = true;
//             }
//         }
//         if (InputType[i] != "dv" && InputType[i] != "iv" && InputType[i] != "rm") {
//             std::cout << InputType[i] << " is an unknown type. Please use either dv, iv, or rm. Aborting Program" << std::endl;
//             return;
//         }
//         if (locate == false) {
//             std::cout << InputType[i] << "/" << InputName[i] << " is not listed in the original matrix. Aborting program" << std::endl;
//             return;
//         }
//     }
//     std::vector<int> RemoveVar;
//     std::vector<std::string> IndepVar;

//     // This checks the input file and tells us what we are removing and changing in the original matrix 
//     //also it edits the CombinedType so we can work independent of the Input file
//     for (int i =0; i < CombinedName.size(); i++) {
//         bool found = false;
//         for (int j = 0; j < InputName.size(); j++) {
//             if (CombinedName[i] == InputName[j]) {
//                 found = true;
//                 if (CombinedType[i] == InputType[j]) {
//                     continue;
//                 }
//                 else {
//                     std::cout << "Changing " << CombinedName[i] << " from " << CombinedType[i] << " to " << InputType[j] << std::endl;
//                     CombinedType[i] = InputType[j];
//                 }
//                 if (InputType[j] == "rm") {
//                     found = false;
//                 }
//             }
//         }
//         if (!found) {
//             std::cout << "Removing " << CombinedName[i] << " from original matrix" << std::endl;
//             RemoveVar.push_back(i);
//         }
//     }

//     //This removes variables that do not appear in the input file but do appear in the original matrix or have rm as the type
//     //It doesn't actually remove anything but moves the unwanted variables on the bottom so we can filter them out in the actual solver code
//     //Note that we can remove from the original matrix but we can't add new variables
//     int t = CombinedName.size() - 1;
//     for (int i = 0; i < RemoveVar.size(); i++) {
//         int rc1 = t;
//         int rc2 = RemoveVar[i];
//         std::string tempType = CombinedType[rc2];
//         CombinedType[rc2] = CombinedType[rc1];
//         CombinedType[rc1] = tempType;
//         std::string tempName = CombinedName[rc2];
//         CombinedName[rc2] = CombinedName[rc1];
//         CombinedName[rc1] = tempName;
//         SwitchRowCol(rc1, rc2);
//         inst.push_back(rc1); inst.push_back(rc2);
//         t -= 1;
//     }

//     //IndepVar lists the # of independent variables so that the sorting loop knows how many row swaps to do
//     for (int i = 0; i < InputType.size(); i++) {
//         if (CombinedType[i] == "iv") {
//             IndepVar.push_back(CombinedName[i]);
//         } }
    
//     //defines # of dependent and independent variables I hope that this can change nY and nP in the actual solve funtion to define the submatrix cuts
//     int nP = IndepVar.size();
//     int nY = InputType.size() - IndepVar.size();

//     //Simple output to make sure that the # of dependent and independent variables are what we expect
//     std::cout << "Original dep and indep sizes are " << indepNames->size() << " " << depNames->size() << std::endl;
//     std::cout << "New dep and indep sizes are " << nP << " " << nY << std::endl;
//     std::cout << "-------------------------------------------------------------" << std::endl;

//     // A simple sorting function that puts any variable marked with iv in the first nP rows and columns
//     // It looks for the iv with the smallest position in the vector (example 0 would be the lowest) and sorts it into the 0th row/columns
//     // It does this for every iv but won't do row operations if it finds that an iv is already in an early row
//     // it then records the row/column operations done so it can undo them when someone plugs in a new set of iv and dv
//     int j = 0;
//     for (int i = 0; i < IndepVar.size(); i++) {
//         int lastI, earlyD, earlyI = -1;
//         for (int k = 0; k < InputType.size(); k++) {
//             if (CombinedType[k] == "iv") {
//                 lastI = k;
//             }
//         }
//         for (int h = j; h < InputType.size(); h++) {
//             if (CombinedType[h] == "iv") {
//                 earlyI = h;
//                 break;
//             }
//         }
//         for (int m = 0; m < InputType.size(); m++) {
//             if (CombinedType[m] == "dv") {
//                 earlyD = m;
//                 break;
//             }
//         }
//         int rc1 = j;
//         int rc2 = earlyI;
//         j += 1;
//         if (rc1 == rc2) {
//             continue;
//         }
//         std::string tempType = CombinedType[rc2];
//         CombinedType[rc2] = CombinedType[rc1];
//         CombinedType[rc1] = tempType;
//         std::string tempName = CombinedName[rc2];
//         CombinedName[rc2] = CombinedName[rc1];
//         CombinedName[rc1] = tempName;
//         SwitchRowCol(rc1, rc2);
//         inst.push_back(rc1); inst.push_back(rc2);
//     }
// }
