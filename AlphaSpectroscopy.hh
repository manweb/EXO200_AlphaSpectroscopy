#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TStyle.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TMath.h"

#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#include "EXOUtilities/EXOWaveformData.hh"
#include "EXOUtilities/EXOWaveform.hh"

using namespace std;

TH1F *hAPDSumWaveform;
TH1F *hAPDSumWaveformBLSub;
TH1F *hWireSumWaveform;
TH1F *hWireSumWaveformBLSub;
TH1F *hWireSumWaveformInv;
TH1F *hWireSumWaveformBLSubInv;

double eLifeTime;

bool IsNoiseEvent(EXOEventData *WFED);
void SubtractBaseline(int *DataSamples, int n_sample, double *bl, double *stdev);
double CorrectZ(double csc, double z);
