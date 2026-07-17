#pragma once

// cpp includes
#include <algorithm>
#include <fstream>
#include <utility>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string>

// root includes
#include "TFile.h"
#include "TVectorD.h"
#include "TMatrixD.h"

void SolveGrandMatrix(char *fname, char *configFile);

int readFile(char *fname);
void openOutputFile(std::string path, std::string fname);
void initMatices(int nP, int nY);
void WriteOutput(int nIndep);

int getDVIV(char *fname);
int insertIV(std::vector<int> indep_indices);

void insertRowCol(int rc1, int rc2);
void SwitchRowCol (int rc1, int rc2, TMatrixD A);
void SwitchMatrices(int rc1, int rc2);

void solve(int nP, int nY);

TFile *grandFile;

TMatrixD mNij, mSij, mMij;
TMatrixD mCij, mVij, mRij;
TMatrixD mVFULL, mRFULL,mSFULL;
TMatrixD sigma_ij, sigma_ji;
TMatrixD mVFULL_clean, mSFULL_clean;

Long64_t fGoodEventNumber;    ///< accumulated so far

/// correlations
TMatrixD mRPY, mRYP;
TMatrixD mRPP, mRYY;
TMatrixD mRYYp;

/// unnormalized covariances
TMatrixD mVPY, mVYP;
TMatrixD mVPP, mVYY;
TMatrixD mVYYp;

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