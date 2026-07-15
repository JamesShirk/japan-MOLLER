#include <fstream>
#include <sstream>
// use this to output the Matrix if you want to view them, do not use "C_ij->Print()"
// Make sure that you set the txt file to the default order before quitting out
// ((TMatrixT<double>*)gDirectory->Get("C_ij"))->Print();


// This switches one row/column with another and writes it to the gDirectory file
void SwitchRowCol (int rc1, int rc2) {
    TMatrixT<double>* mCij = (TMatrixT<double>*)gDirectory->Get("C_ij");
    TMatrixT<double>* mRij = (TMatrixT<double>*)gDirectory->Get("R_ij");
    TMatrixT<double>* mVij = (TMatrixT<double>*)gDirectory->Get("V_ij");
    std::cout << rc1 << " " << rc2 << std::endl;
    for (int col = 0; col < mCij->GetNcols(); col++) {
        double temp = (*mCij)(rc1, col);
        (*mCij)(rc1, col) = (*mCij)(rc2, col);
        (*mCij)(rc2, col) = temp;
    }
    for (int row = 0; row < mCij->GetNrows(); row++) {
        double temp = (*mCij)(row, rc1);
        (*mCij)(row, rc1) = (*mCij)(row, rc2);
        (*mCij)(row, rc2) = temp;
    }
    for (int col = 0; col < mRij->GetNcols(); col++) {
        double temp = (*mRij)(rc1, col);
        (*mRij)(rc1, col) = (*mRij)(rc2, col);
        (*mRij)(rc2, col) = temp;
    }
    for (int row = 0; row < mRij->GetNrows(); row++) {
        double temp = (*mRij)(row, rc1);
        (*mRij)(row, rc1) = (*mRij)(row, rc2);
        (*mRij)(row, rc2) = temp;
    }
    for (int col = 0; col < mVij->GetNcols(); col++) {
        double temp = (*mVij)(rc1, col);
        (*mVij)(rc1, col) = (*mVij)(rc2, col);
        (*mVij)(rc2, col) = temp;
    }
    for (int row = 0; row < mVij->GetNrows(); row++) {
        double temp = (*mVij)(row, rc1);
        (*mVij)(row, rc1) = (*mVij)(row, rc2);
        (*mVij)(row, rc2) = temp;
    }
    mCij->Write("C_ij", TObject::kOverwrite);
    mRij->Write("R_ij", TObject::kOverwrite);
    mVij->Write("V_ij", TObject::kOverwrite);
}

// vector I define to help "remember" the actions taken in one ChangeIVDV run so that it can automatically reset when you run the program again
// if a better method of this exists please let me know
std::vector<int> inst;

// what you actually run to change the rows and columns
// the input file needs to have the same format as the matricies with iv and dv in from and you simply change some from iv to dv or vice versa to make changes to the variables
void ChangeIVDV () {

// grabs the raw names and categorizes them
    std::vector<std::string> *depNames = (std::vector<std::string>*)gDirectory->Get("Dependent_Names");
    std::vector<std::string> *indepNames = (std::vector<std::string>*)gDirectory->Get("Independent_Names");
    std::vector<std::string> CombinedName;
    std::vector<std::string> CombinedType;
    for (int i = 0; i < indepNames->size(); i++) {
        CombinedName.push_back((*indepNames)[i]); 
        CombinedType.push_back("iv"); }
    for (int i = 0; i < depNames->size(); i++) {
        CombinedName.push_back((*depNames)[i]);
        CombinedType.push_back("dv"); }

// grabs the input file and makes two lists for type of variables ("iv" or "dv") and the name for each variable
    std::ifstream file("rootScripts/Correlator.txt");
    std::string line;
    std::vector<std::string> InputType, InputName, NP;
    if (!file.is_open()) {
        std::cout << "FAILED TO OPEN FILE!" << std::endl;
        return;
    }
    int row1, row2;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type, name;
        iss >> type >> name;
        name.erase(0, 5);
        InputName.push_back(name);
        InputType.push_back(type); 
    }

    // This undoes the last row and column swaps done so you don't have to run the program twice to get the original matrix back
    if (inst.size() >= 1) {
        for (int i = 0; i < inst.size(); i += 2) {
            int cr1 = inst[inst.size() - i - 1]; int cr2 = inst[inst.size() - i - 2];
            SwitchRowCol(cr1, cr2);
        }
    }
    inst.clear();

    //Making sure the the variable actually existed in the original matrix and that they have the correct type
    for (int i = 0; i < InputName.size(); i++) {
        bool locate = false;
        for (int j = 0; j < CombinedName.size(); j++) {
            if (InputName[i] == CombinedName[j]) {
                locate = true;
            }
        }
        if (InputType[i] != "dv" && InputType[i] != "iv" && InputType[i] != "rm") {
            std::cout << InputType[i] << " is an unknown type. Please use either dv, iv, or rm. Aborting Program" << std::endl;
            return;
        }
        if (locate == false) {
            std::cout << InputType[i] << "/" << InputName[i] << " is not listed in the original matrix. Aborting program" << std::endl;
            return;
        }
    }
    std::vector<int> RemoveVar;
    std::vector<std::string> IndepVar;

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
                    std::cout << "Changing " << CombinedName[i] << " from " << CombinedType[i] << " to " << InputType[j] << std::endl;
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
    for (int i = 0; i < RemoveVar.size(); i++) {
        int rc1 = t;
        int rc2 = RemoveVar[i];
        std::string tempType = CombinedType[rc2];
        CombinedType[rc2] = CombinedType[rc1];
        CombinedType[rc1] = tempType;
        std::string tempName = CombinedName[rc2];
        CombinedName[rc2] = CombinedName[rc1];
        CombinedName[rc1] = tempName;
        SwitchRowCol(rc1, rc2);
        inst.push_back(rc1); inst.push_back(rc2);
        t -= 1;
    }

    //IndepVar lists the # of independent variables so that the sorting loop knows how many row swaps to do
    for (int i = 0; i < InputType.size(); i++) {
        if (CombinedType[i] == "iv") {
            IndepVar.push_back(CombinedName[i]);
        } }
    
    //defines # of dependent and independent variables I hope that this can change nY and nP in the actual solve funtion to define the submatrix cuts
    int nP = IndepVar.size();
    int nY = InputType.size() - IndepVar.size();

    //Simple output to make sure that the # of dependent and independent variables are what we expect
    std::cout << "Original dep and indep sizes are " << indepNames->size() << " " << depNames->size() << std::endl;
    std::cout << "New dep and indep sizes are " << nP << " " << nY << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;

    // A simple sorting function that puts any variable marked with iv in the first nP rows and columns
    // It looks for the iv with the smallest position in the vector (example 0 would be the lowest) and sorts it into the 0th row/columns
    // It does this for every iv but won't do row operations if it finds that an iv is already in an early row
    // it then records the row/column operations done so it can undo them when someone plugs in a new set of iv and dv
    int j = 0;
    for (int i = 0; i < IndepVar.size(); i++) {
        int lastI, earlyD, earlyI = -1;
        for (int k = 0; k < InputType.size(); k++) {
            if (CombinedType[k] == "iv") {
                lastI = k;
            }
        }
        for (int h = j; h < InputType.size(); h++) {
            if (CombinedType[h] == "iv") {
                earlyI = h;
                break;
            }
        }
        for (int m = 0; m < InputType.size(); m++) {
            if (CombinedType[m] == "dv") {
                earlyD = m;
                break;
            }
        }
        int rc1 = j;
        int rc2 = earlyI;
        j += 1;
        if (rc1 == rc2) {
            continue;
        }
        std::string tempType = CombinedType[rc2];
        CombinedType[rc2] = CombinedType[rc1];
        CombinedType[rc1] = tempType;
        std::string tempName = CombinedName[rc2];
        CombinedName[rc2] = CombinedName[rc1];
        CombinedName[rc1] = tempName;
        SwitchRowCol(rc1, rc2);
        inst.push_back(rc1); inst.push_back(rc2);
    }
}
