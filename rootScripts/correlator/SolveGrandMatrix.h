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

void SolveGrandMatrix(char *fname, int nIndep);
int readFile(char *fname, int nIndep);
void WriteOutput(int nIndep);
void solve(int nP, int nY);
void openOutputFile(std::string path, std::string fname = "regression_output.root");
void initMatices(int nP, int nY);

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

std::vector<std::string> *allnames;

TFile *ofile;