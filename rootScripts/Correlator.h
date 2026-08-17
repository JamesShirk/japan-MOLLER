#pragma once

// cpp includes
#include <algorithm>
#include <fstream>
#include <sstream>
#include <utility>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string>

// root includes
#include "TFile.h"
#include "TVectorD.h"
#include "TMatrixD.h"

TFile *_file0 = TFile::Open("blueR_grand4new.slope.root", "UPDATE");

void MakeTemp(std::string name);
void initialize();
void SwitchM(TMatrixT<double> &matrix, int &rc1, int &rc2);
template <typename type>
void SwitchV(type &vector, int &elem1, int &elem2);
void SwitchRowCol (int rc1, int rc2);
void SwitchElem(int elem1, int elem2);
void insertRowCol(int rc1, int rc2);
void Delete();
void Reset();
void resize();
void solve(int nP, int nY);

TMatrixD mVFULL, mRFULL,mSFULL;
TMatrixD mVFULL_clean, mSFULL_clean;

Long64_t fGoodEventNumber;    ///< accumulated so far

/// correlations
TMatrixD mRPY, mRYP;
TMatrixD mRPP, mRYY;
TMatrixD mRYYp;

TMatrixD mMPP, mMYY;

/// unnormalized covariances
TMatrixD mVPY, mVYP;
TMatrixD mVPP, mVYY;
TMatrixD mVYYp;

TMatrixD sigmaP, sigmaY;
TMatrixD mVPY_clean, mVPP_clean, mVYY_clean, mVYP_clean;
TMatrixD mSPY_clean, mSPP_clean, mSYY_clean, mSYP_clean;

TVectorD mSP_clean, mSY_clean;

/// variances
TVectorD mVP, mVY;
TVectorD mVYp;

/// normalized covariances
TMatrixD mSPY, mSYP;
TMatrixD mSPP, mSYY;
TMatrixD mSYYp;

/// sigmas
TVectorD mSP, mSY;
TVectorD mSYp;

/// mean values
TVectorD mMP, mMY, mMYp;
TMatrixD Axy, Ayx, dAxy, dAyx;

std::vector<std::string> allnames;

TFile *ofile;

void Commit();
void ChangeIVDV ();

std::vector<std::string> Operations;

std::vector<std::string> CombinedName, CombinedType;
std::string line;
std::vector<std::string> InputType, InputName, NP;

int row1, row2;

std::string type, name;

std::vector<int> RemoveVar;
std::vector<std::string> IndepVar;

int nP, nY, nR, rc1, rc2, j;

void List();