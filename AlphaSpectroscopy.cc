// ************************************************************************************************
// Name:	AlphaSpectroscopy.C
// Author:	Manuel Weber
// Input:	EXO200 Root data files
// Ouput:	Root file containing information on found alpha events. Each entry corresponds
//		to one found alpha
// Description:	This script searches for alpha events in background runs. The energy of the
//		charge clusters are corrected for purity and a w value of 20 eV. Optionally
//		it also corrects for the individual u-wire gains (relative). To exclude the
//		events on the cathode a fiducial cut (-180 < z < -10, 10 < z < 180) is applied.
//
// ************************************************************************************************

#include "AlphaSpectroscopy.hh"

int main(int argc, char* argv[])
{
  if (argc != 2) {cout << "Usage: ./AlphaSpectroscopy output.root" << endl; return 0;}

  char *oname = argv[1];
  //double eLifeTime = atof(argv[3]);

  // create TChain for reconstructed data
  TChain *AllFiles = new TChain("tree","data tree");

  //char *fname = "/home/exodaq/EXO200Analysis/DataFiles/PhysicsRunsGood.dat";
  char *fname = "/nfs/slac/g/exo/maweber/EXO200Analysis/AlphaSpectroscopy/scripts/AlphaSpectroscopy.dat";

  // Get data files
  string fnameTMP;
  ifstream ifs(fname);
  while (getline(ifs, fnameTMP)) {
     std::cout << "Adding " << fnameTMP << " to the chain" << std::endl;
     AllFiles->Add(fnameTMP.c_str());
  }

  int nentries = AllFiles->GetEntries();

  cout << "Total number of events in run: " << nentries << endl;

  // output tree
  TTree *oTree = new TTree("t","Alpha Spectroscopy");

  int fnr;
  int fne;
  int fncl;
  int ftrigsec;
  int ftrigsub;
  int ftrigoff;
  double ftsc;
  double ftcl;
  double fX;
  double fY;
  double fZ;
  double fXsc;
  double fYsc;
  double fZsc;
  double fc1sc;
  double fc2sc;
  double fcsc;
  double fcscCorr;
  double feccl;
  double fercl;
  double fCL;
  double fAlphaEnergy;

  oTree->Branch("fnr",&fnr,"fnr/I");
  oTree->Branch("fne",&fne,"fne/I");
  oTree->Branch("fncl",&fncl,"fncl/I");
  oTree->Branch("ftrigsec",&ftrigsec,"ftrigsec/I");
  oTree->Branch("ftrigsub",&ftrigsub,"ftrigsub/I");
  oTree->Branch("ftrigoff",&ftrigoff,"ftrigoff/I");
  oTree->Branch("ftsc",&ftsc,"ftsc/D");
  oTree->Branch("ftcl",&ftcl,"ftcl/D");
  oTree->Branch("fX",&fX,"fX/D");
  oTree->Branch("fY",&fY,"fY/D");
  oTree->Branch("fZ",&fZ,"fZ/D");
  oTree->Branch("fXsc",&fXsc,"fXsc/D");
  oTree->Branch("fYsc",&fYsc,"fYsc/D");
  oTree->Branch("fZsc",&fZsc,"fZsc/D");
  oTree->Branch("fc1sc",&fc1sc,"fc1sc/D");
  oTree->Branch("fc2sc",&fc2sc,"fc2sc/D");
  oTree->Branch("fcsc",&fcsc,"fcsc/D");
  oTree->Branch("fcscCorr",&fcscCorr,"fcscCorr/D");
  oTree->Branch("feccl",&feccl,"feccl/D");
  oTree->Branch("fercl",&fercl,"fercl/D");
  oTree->Branch("fCL",&fCL,"fCL/D");
  oTree->Branch("fAlphaEnergy",&fAlphaEnergy,"fAlphaEnergy/D");

  // set branch address for event data
  EXOEventData* ED = 0;
  AllFiles->SetBranchAddress("EventBranch", &ED);

  // loop over all entries
  for (int i = 0; i < nentries; i++) {
     //cout << "begin of event " << i << endl;
     if (i%1000 == 0) {std::cout << i << " events processed   eLifeTime = " << eLifeTime << std::endl;}

     AllFiles->GetEntry(i);

     // check if this is a noise event
     if (ED->fEventHeader.fTaggedAsNoise == 1) {continue;}
     if (ED->fEventHeader.fIndividualTriggerRequest == 0 && ED->fEventHeader.fSumTriggerRequest == 0) {continue;}

     // get electron lifetime according to fit
     //double purFitP0 = 293.8;
     //double purFitP1 = -1.68;
     double purFitP0 = 0.0;
     double purFitP1 = 0.0;
     double purFitP2 = 0.0;
     double purFitP3 = 0.0;
     double purFitP4 = 0.0;

     double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.0;

    if (purTime < 23) {
      purFitP0 = 250.0;
      purFitP1 = 0.0;
      purFitP2 = 0.0;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }
    if (purTime >= 23 && purTime < 58) {
      purFitP0 = -284.596;
      purFitP1 = 53.6978;
      purFitP2 = -1.88664;
      purFitP3 = 0.0269101;
      purFitP4 = -0.000133772;
    }
    if (purTime >= 58 && purTime < 81.6) {
      purFitP0 = 14068.5;
      purFitP1 = -908.011;
      purFitP2 = 21.8864;
      purFitP3 = -0.230994;
      purFitP4 = 0.00090631;
    }
    if (purTime >= 81.6 && purTime < 94.0) {
      purFitP0 = -9011.55;
      purFitP1 = 115.417;
      purFitP2 = 0.0;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }
    if (purTime >= 94.0 && purTime < 102.5) {
      purFitP0 = 2000.0;
      purFitP1 = 0.0;
      purFitP2 = 0.0;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }
    if (purTime >= 102.5 && purTime < 113.0) {
      purFitP0 = -1208000.0;
      purFitP1 = 34380.0;
      purFitP2 = -325.9;
      purFitP3 = 1.03;
      purFitP4 = 0.0;
    }
    if (purTime >= 113.0 && purTime < 129.6) {
      purFitP0 = -48740.0;
      purFitP1 = 805.0;
      purFitP2 = -3.259;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }
    if (purTime >= 129.6 && purTime < 142.0) {
      purFitP0 = -29510.0;
      purFitP1 = 230.1;
      purFitP2 = 0.0;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }
    if (purTime >= 142.0) {
      purFitP0 = 4000.0;
      purFitP1 = 0.0;
      purFitP2 = 0.0;
      purFitP3 = 0.0;
      purFitP4 = 0.0;
    }

     eLifeTime = purFitP4*purTime*purTime*purTime*purTime + purFitP3*purTime*purTime*purTime + purFitP2*purTime*purTime + purFitP1*purTime + purFitP0;

     int nsc = ED->GetNumScintillationClusters();

     //if (nsc != 1) {continue;}

     // loop over scintillation clusters
     for (int sID = 0; sID < nsc; sID++) {
        EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(sID);

        double tsc1 = scint_cluster->fTime;

        bool GoodScintCluster = true;
        // make sure no other scintillation cluster is within 110us
        for (int i = 0; i < nsc; i++) {
           if (i == sID) {continue;}

           double tsc2 = ED->GetScintillationCluster(i)->fTime;
           if (TMath::Abs(tsc2 - tsc1) < 110000) {GoodScintCluster = false; break;}
        }

        if (!GoodScintCluster) {continue;}
        if (scint_cluster->fTime > 1928000) {continue;}

        double Scint1 = scint_cluster->fCountsSumOnAPDPlaneOne;
        double Scint2 = scint_cluster->fCountsSumOnAPDPlaneTwo;
        double SumScint = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;
        double SumScintCorr = SumScint;

        // cut on number photon counts
        //if (SumScint < 6000) {continue;}

        // get number of associated charge cluster
        int ncl = scint_cluster->GetNumChargeClusters();

        // apply single site cuts
        if (ncl > 1) {continue;}
        //if (scint_cluster->fTime > 1900000) {continue;}

        fncl = 0;
        double clX = -999.9;
        double clY = -999.9;
        double clZ = -999.9;
        double ChargeCorr = 0.0;
        double ChargeRaw = 0.0;
        ftcl = -999.9;
        EXOChargeCluster *charge_cluster = 0;
        if (ncl == 1) {
           fncl = 1;

           // get associated charge cluster
           charge_cluster = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0);

           if (!(charge_cluster->fIs3DCluster)) {continue;} // all three coordinates must be reconstructed

           // fiducial cut
           clX = charge_cluster->fX;
           clY = charge_cluster->fY;
           clZ = charge_cluster->fZ;
           if (clX < -200 || clX > 200) {continue;}
           if (clY < -200 || clY > 200) {continue;}
           if (clZ < -200 || clZ > 200) {continue;}
           //if (clZ > -12 && clZ < 12) {continue;}
           //if (clZ < -150) {continue;}
           //if (clZ > 150) {continue;}

           if (TMath::Sqrt(clX*clX + clY*clY) < 163 && ((clZ > -172 && clZ < -20) || (clZ > 20 && clZ < 172))) {SumScintCorr = CorrectZ(SumScint, clZ);}

           double dt = charge_cluster->fDriftTime / 1000; // drift time of the charge cluster
           //ChargeCorr = charge_cluster->fCorrectedEnergy * TMath::Exp(dt/eLifeTime); // purity corrected energy
           ChargeCorr = charge_cluster->fPurityCorrectedEnergy;
           ChargeRaw = charge_cluster->fRawEnergy;
        }

        // Fill histograms and write tree
        double CLThresh = 0.05;
        if (ED->fRunNumber >= 2379) {CLThresh = 0.01;}
        if (SumScint > 0.0) {
           double CL = ChargeRaw/SumScint;
           if (CL > CLThresh) {continue;}

           double factor = 0.0;
           if (clZ > -150 && clZ < -12 || clZ > 12 && clZ < 150) {factor = 0.2606;}
           else {if (clZ < -150 || clZ > 150) {factor = 0.40205;}
                 else {factor = 0.40986;}
           }

           fnr = ED->fRunNumber;
           fne = ED->fEventNumber;
           ftrigsec = ED->fEventHeader.fTriggerSeconds;
           ftrigsub = ED->fEventHeader.fTriggerMicroSeconds;
           ftrigoff = ED->fEventHeader.fTriggerOffset;
           ftsc = scint_cluster->fTime;
           fX = clX;
           fY = clY;
           fZ = clZ;
           fXsc = scint_cluster->fX;
           fYsc = scint_cluster->fY;
           fZsc = scint_cluster->fZ;
           fc1sc = Scint1;
           fc2sc = Scint2;
           fcsc = SumScint;
           fcscCorr = SumScintCorr;
           feccl = ChargeCorr;
           fercl = ChargeRaw;
           fCL = CL;
           fAlphaEnergy = ChargeCorr + factor*SumScint;
           if (ncl == 1) {ftcl = charge_cluster->fCollectionTime;}
           else {ftcl = -999.9;}

           oTree->Fill();
        }
     }
  }

  // write tree to file
  TFile *oFile = new TFile(oname,"RECREATE");
  oTree->Write();

  oFile->Close();

  return 1;
}

double CorrectZ(double csc, double z)
{
/*  double a1 = -11.16;
  double a2 = 29.46;
  double b1 = 21112.1;
  double b2 = 21200.6;
*/

  double a1 = -11.4924;
  double a2 = 30.624;
  double b1 = 21099.6;
  double b2 = 21223.5;

  double z0 = 96.0; 
  double csc_middle = ((a2 - a1)*z0 + b2 + b1) / 2.0;
  double corr1 = csc_middle / (-1.0*a1*z0 + b1);
  double corr2  = csc_middle / (a2*z0 + b2);

  double b = 0.0;
  double c = 0.0;
  double csc_corr = 0.0;

  if (z < 0) {
     b = csc - a1 * z;
     c = csc / (-1.0*a1*z0 + b);
     csc_corr = csc / c * corr1;
  }
  else {
     b = csc - a2 * z;
     c = csc / (a2*z0 + b);
     csc_corr = csc / c * corr2;
  }

  return csc_corr;
}

