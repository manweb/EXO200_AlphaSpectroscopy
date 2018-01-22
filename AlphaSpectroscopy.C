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

TH1F *hAPDSumWaveform = new TH1F("hAPDSumWaveform", "APD Sum Waveform",2048,0,2048);
TH1F *hAPDSumWaveformBLSub = new TH1F("hAPDSumWaveformBLSub", "APD Sum Waveform - baseline > threshold",2048,0,2048);
TH1F *hWireSumWaveform = new TH1F("hWireSumWaveform", "Wire Sum Waveform",2048,0,2048);
TH1F *hWireSumWaveformBLSub = new TH1F("hWireSumWaveformBLSub", "Wire Sum Waveform - baseline > threshold",2048,0,2048);
TH1F *hWireSumWaveformInv = new TH1F("hWireSumWaveformInv", "Inverted Wire Sum Waveform",2048,0,2048);
TH1F *hWireSumWaveformBLSubInv = new TH1F("hWireSumWaveformBLSubInv", "Inverted Wire Sum Waveform - baseline > threshold",2048,0,2048);

void AlphaSpectroscopy()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  gSystem->Load("libEXOUtilities");

  double w_correction = 1.0;
  //double eLifeTime = 290;
  double eLifeTime = 240;

  // initialize wavform histogram
  //hAPDSumWaveform = new TH1F("hAPDSumWaveform", "APD Sum Waveform",2048,0,2048);
  //hAPDSumWaveformBLSub = new TH1F("hAPDSumWaveformBLSub", "APD Sum Waveform - baseline > threshold",2048,0,2048);
  //hWireSumWaveform = new TH1F("hWireSumWaveform", "Wire Sum Waveform",2048,0,2048);
  //hWireSumWaveformBLSub = new TH1F("hWireSumWaveformBLSub", "Wire Sum Waveform - baseline > threshold",2048,0,2048);
  //hWireSumWaveformInv = new TH1F("hWireSumWaveformInv", "Inverted Wire Sum Waveform",2048,0,2048);
  //hWireSumWaveformBLSubInv = new TH1F("hWireSumWaveformBLSubInv", "Inverted Wire Sum Waveform - baseline > threshold",2048,0,2048);

  // create TChain for reconstructed data
  TChain *AllFiles = new TChain("tree","t");

  // create TChain for corresponding waveform data
  TChain *WFFiles = new TChain("tree","t");

  //AllFiles->Add("/EXO200Data/Disk4/processed/838/run00000838-000.root");
  //AllFiles->Add("/EXO200Data/Disk4/processed/858/run00000858-000.root");
  //AllFiles->Add("/EXO200Data/Disk4/processed/858/run00000858-001.root");

  AllFiles->Add("/EXO200Data/Disk5/processed/1569/recon00001569-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1570/recon00001570-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1581/recon00001581-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1591/recon00001591-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1596/recon00001596-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1606/recon00001606-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1620/recon00001620-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1621/recon00001621-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1622/recon00001622-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1628/recon00001628-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1635/recon00001635-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1639/recon00001639-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1644/recon00001644-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1655/recon00001655-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1660/recon00001660-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1661/recon00001661-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1666/recon00001666-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1670/recon00001670-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1679/recon00001679-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1683/recon00001683-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1687/recon00001687-000.root");
  AllFiles->Add("/EXO200Data/Disk5/processed/1691/recon00001691-000.root");

  WFFiles->Add("/EXO200Data/Disk5/processed/1569/run00001569-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1570/run00001570-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1581/run00001581-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1591/run00001591-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1596/run00001596-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1606/run00001606-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1620/run00001620-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1621/run00001621-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1622/run00001622-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1628/run00001628-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1635/run00001635-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1639/run00001639-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1644/run00001644-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1655/run00001655-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1660/run00001660-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1661/run00001661-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1666/run00001666-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1670/run00001670-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1679/run00001679-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1683/run00001683-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1687/run00001687-000.root");
  WFFiles->Add("/EXO200Data/Disk5/processed/1691/run00001691-000.root");

  int nentries = AllFiles->GetEntries();

  cout << "Total number of events in run: " << nentries << endl;

  // output tree
  TTree *oTree = new TTree("t","Alpha Spectroscopy");

  int fnr;
  int fne;
  int ftrigsec;
  int ftrigsub;
  int ftrigoff;
  double ftsc;
  double fX;
  double fY;
  double fZ;
  double fcsc;
  double feccl;
  double fCL;
  double fAlphaEnergy;

  oTree->Branch("fnr",&fnr,"fnr/I");
  oTree->Branch("fne",&fne,"fne/I");
  oTree->Branch("ftrigsec",&ftrigsec,"ftrigsec/I");
  oTree->Branch("ftrigsub",&ftrigsub,"ftrigsub/I");
  oTree->Branch("ftrigoff",&ftrigoff,"ftrigoff/I");
  oTree->Branch("ftsc",&ftsc,"ftsc/D");
  oTree->Branch("fX",&fX,"fX/D");
  oTree->Branch("fY",&fY,"fY/D");
  oTree->Branch("fZ",&fZ,"fZ/D");
  oTree->Branch("fcsc",&fcsc,"fcsc/D");
  oTree->Branch("feccl",&feccl,"feccl/D");
  oTree->Branch("fCL",&fCL,"fCL/D");
  oTree->Branch("fAlphaEnergy",&fAlphaEnergy,"fAlphaEnergy/D");

  // set branch address for event data
  EXOEventData* ED = 0;
  AllFiles->SetBranchAddress("EventBranch", &ED);

  EXOEventData *WFED = 0;
  WFFiles->SetBranchAddress("EventBranch",&WFED);

  WFFiles->BuildIndex("EventBranch.fEventNumber","EventBranch.fRunNumber");

  // create histograms
  TH1F *hECorrB = new TH1F("hECorrB","Corrected Cluster Energy (bulk)",200,0,1500);
  TH1F *hECorrC = new TH1F("hECorrC","Corrected Cluster Energy (cathode)",200,0,1500);
  TH1F *hECorrS = new TH1F("hECorrS","Corrected Cluster Energy",200,0,1500);
  TH1F *hECorrA = new TH1F("hECorrA","Corrected Cluster Energy (anode)",200,0,1500);

  hECorrB->SetLineColor(kBlack);
  hECorrB->SetLineStyle(2);
  hECorrC->SetLineColor(kRed);
  hECorrC->SetLineStyle(2);
  hECorrA->SetLineColor(kBlue);
  hECorrA->SetLineStyle(2);
  hECorrS->SetFillColor(17);
  hECorrS->GetXaxis()->SetTitle("rec. energy [keV]");

  TH1F *hScintB = new TH1F("hScintB","Scintillation (bulk)",200,0,40000);
  TH1F *hScintC = new TH1F("hScintC","Scintillation (cathode)",200,0,40000);
  TH1F *hScintS = new TH1F("hScintS","Scintillation",200,0,40000);
  TH1F *hScintA = new TH1F("hScintA","Scintillation (anode",200,0,40000);

  hScintB->SetLineColor(kBlack);
  hScintB->SetLineStyle(2);
  hScintC->SetLineColor(kRed);
  hScintC->SetLineStyle(2);
  hScintA->SetLineColor(kBlue);
  hScintA->SetLineStyle(2);
  hScintS->SetFillColor(17);
  hScintS->GetXaxis()->SetTitle("scintillation [photon counts]");

  TH1F *hIonScintCorrB = new TH1F("hIonSintCorrB","Ion/scint (bulk)",100,0,1);
  TH1F *hIonScintCorrC = new TH1F("hIonSintCorrC","Ion/scint (cathode)",100,0,1);
  TH1F *hIonScintCorrS = new TH1F("hIonSintCorrS","Ion/scint",100,0,1);
  TH1F *hIonScintCorrA = new TH1F("hIonSintCorrA","Ion/scint (anode)",100,0,1);

  hIonScintCorrB->SetLineColor(kBlack);
  hIonScintCorrB->SetLineStyle(2);
  hIonScintCorrC->SetLineColor(kRed);
  hIonScintCorrC->SetLineStyle(2);
  hIonScintCorrA->SetLineColor(kBlue);
  hIonScintCorrA->SetLineStyle(2);
  hIonScintCorrS->SetFillColor(17);

  TH1F *hAlphaSpectrumB = new TH1F("hAlphaSpectrumB","Alpha energy spectrum (bulk)",200,0,10000);
  TH1F *hAlphaSpectrumC = new TH1F("hAlphaSpectrumC","Alpha energy spectrum (cathode)",200,0,10000);
  TH1F *hAlphaSpectrumS = new TH1F("hAlphaSpectrumS","Alpha energy spectrum",200,0,10000);
  TH1F *hAlphaSpectrumA = new TH1F("hAlphaSpectrumA","Alpha energy spectrum (anode)",200,0,10000);

  hAlphaSpectrumB->SetLineColor(kBlack);
  hAlphaSpectrumB->SetLineStyle(2);
  hAlphaSpectrumC->SetLineColor(kRed);
  hAlphaSpectrumC->SetLineStyle(2);
  hAlphaSpectrumA->SetLineColor(kBlue);
  hAlphaSpectrumA->SetLineStyle(2);
  hAlphaSpectrumS->SetFillColor(17);
  hAlphaSpectrumS->GetXaxis()->SetTitle("energy [keV]");

  TH2F *hIonScintECorr = new TH2F("hIonScintECorr","Ion/Scint vs energy",200,0,1500,100,0,1);

  hScintC->SetLineColor(kRed);

  TH1F *hX = new TH1F("hX","X position of alpha events",40,-200,200);
  TH1F *hY = new TH1F("hY","Y position of alpha events",40,-200,200);
  TH1F *hZ = new TH1F("hZ","Z position of alpha events",40,-200,200);
  TH2F *hXY = new TH2F("hXY","X-Y position of alhpa events",40,-200,200,40,-200,200);

  hX->GetXaxis()->SetTitle("x [mm]");
  hY->GetXaxis()->SetTitle("y [mm]");
  hZ->GetXaxis()->SetTitle("z [mm]");
  hXY->GetXaxis()->SetTitle("x [mm]");
  hXY->GetYaxis()->SetTitle("y [mm]");

  TH1F *hUHITPZ = new TH1F("hUHITPZ","U wire hits (PZ)",38,0,37);
  TH1F *hUHITNZ = new TH1F("hUHITNZ","U wire hits (NZ)",38,76,113);

  hUHITPZ->GetXaxis()->SetTitle("channel");
  hUHITNZ->GetXaxis()->SetTitle("channel");

  TGraph *grB = new TGraph(1000);
  grB->GetXaxis()->SetTitle("energy [keV]");
  grB->GetXaxis()->SetRangeUser(0,1000);
  grB->GetYaxis()->SetTitle("ion/scint [a.u.]");
  grB->GetYaxis()->SetRangeUser(0,1);
  grB->SetMarkerStyle(1);
  //grB->SetMarkerSize(0.5);

  TGraph *grC = new TGraph(1000);
  grC->GetXaxis()->SetTitle("energy [keV]");
  grC->GetXaxis()->SetRangeUser(0,1000);
  grC->GetYaxis()->SetTitle("ion/scint [a.u.]");
  grC->GetYaxis()->SetRangeUser(0,1);
  grC->SetMarkerStyle(1);
  //grC->SetMarkerSize(0.5);
  grC->SetMarkerColor(kRed);

  TGraph *grA = new TGraph(1000);
  grA->GetXaxis()->SetTitle("energy [keV]");
  grA->GetXaxis()->SetRangeUser(0,1000);
  grA->GetYaxis()->SetTitle("ion/scint [a.u.]");
  grA->GetYaxis()->SetRangeUser(0,1);
  grA->SetMarkerStyle(1);
  //grA->SetMarkerSize(0.5);
  grA->SetMarkerColor(kBlue);

  EXO3DView *v = new EXO3DView();

  int nAlphas = 0;
  int nPointsB = 0;
  int nPointsC = 0;
  int nPointsA = 0;
  // loop over all entries
  for (int i = 0; i < nentries; i++) {
     if (i%1000 == 0) {cout << i << " events processed" << endl;}

     AllFiles->GetEntry(i);

     int runID = ED->fRunNumber;
     int evtID = ED->fEventNumber;

     WFFiles->GetEntryWithIndex(evtID,runID);

     if (IsNoiseEvent(WFED)) {continue;}

     int nsc = ED->GetNumScintillationClusters();
     int nclAll = ED->GetNumChargeClusters();

     // cut on total number of charge clusters to reject noise events
     if (nclAll > 20) {continue;}

     // loop over scintillation clusters
     for (int sID = 0; sID < nsc; sID++) {
        EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(sID);

        double SumScint = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;

        // cut on number photon counts
        if (SumScint < 6000) {continue;}

        // loop over associated charge clusters and sum their energy
        int ncl = scint_cluster->GetNumChargeClusters();

        // apply single site cuts
        if (ncl != 1) {continue;}
        if (scint_cluster->fTime > 1900000) {continue;}

        // get associated charge cluster
        EXOChargeCluster *charge_cluster = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0);

        if (!(charge_cluster->fIs3DCluster)) {continue;} // all three coordinates must be reconstructed

        // fiducial cut
        double clX = charge_cluster->fX;
        double clY = charge_cluster->fY;
        double clZ = charge_cluster->fZ;
        if (clX < -200 || clX > 200) {continue;}
        if (clY < -200 || clY > 200) {continue;}
        if (clZ < -200 || clZ > 200) {continue;}
        //if (clZ > -12 && clZ < 12) {continue;}
        //if (clZ < -150) {continue;}
        //if (clZ > 150) {continue;}

        double dt = charge_cluster->fDriftTime / 1000; // drift time of the charge cluster

        /*int nUWires = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0)->GetNumUWireSignals();
        if (nUWires != 1) {continue;} // number u wires must be 1

        int chID = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0)->GetUWireSignalChannelAt(0);
        int dHalf = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0)->fDetectorHalf;
        if (!(dHalf == 0 || dHalf == 1)) {continue;} // detector half must either be 0 or 1

        if (dHalf == 1) {chID = chID - 76;}
        double gain_corr = GetWireGain(dHalf,chID);
        if (gain_corr == -999) {continue;}*/

        double ChargeCorr = charge_cluster->fCorrectedEnergy * TMath::Exp(dt/eLifeTime) * w_correction; // corrected energy for 290 us electron lifetime

        //double ChargeCorr = ED->GetScintillationCluster(sID)->GetChargeClusterAt(0)->fCorrectedEnergy / gain_corr * TMath::Exp(dt/eLifeTime) * w_correction; // corrected energy for 290 us electron lifetime and wire gains

        // Fill histograms and write tree
        if (SumScint > 0.0) {
           double CL = ChargeCorr/SumScint;
           if (CL > 0.05) {continue;}

           double factor = 0.0;
           if (clZ > -150 && clZ < -12 || clZ > 12 && clZ < 150) {
              hECorrB->Fill(ChargeCorr);
              hScintB->Fill(SumScint);
              hIonScintCorrB->Fill(CL);
              hScintB->Fill(SumScint);
              hAlphaSpectrumB->Fill(ChargeCorr + 0.2606*SumScint);
              hAlphaSpectrumS->Fill(ChargeCorr + 0.2606*SumScint);
              factor = 0.2606;

              grB->SetPoint(nPointsB,ChargeCorr,CL);
              nPointsB++;
           }
           else {if (clZ < -150 || clZ > 150) {
              hECorrA->Fill(ChargeCorr);
              hScintA->Fill(SumScint);
              hIonScintCorrA->Fill(CL);
              hScintA->Fill(SumScint);
              hAlphaSpectrumA->Fill(ChargeCorr + 0.40986*SumScint);
              hAlphaSpectrumS->Fill(ChargeCorr + 0.40986*SumScint);
              factor = 0.40205;

              grA->SetPoint(nPointsA,ChargeCorr,CL);
              nPointsA++;
              }
              else {
              hECorrC->Fill(ChargeCorr);
              hScintC->Fill(SumScint);
              hIonScintCorrC->Fill(CL);
              hScintC->Fill(SumScint);
              hAlphaSpectrumC->Fill(ChargeCorr + 0.40986*SumScint);
              hAlphaSpectrumS->Fill(ChargeCorr + 0.40986*SumScint);
              factor = 0.40986;

              grC->SetPoint(nPointsC,ChargeCorr,CL);
              nPointsC++;
              }
           }

           hECorrS->Fill(ChargeCorr);
           hScintS->Fill(SumScint);
           hIonScintCorrS->Fill(CL);
           hScintS->Fill(SumScint);
           hIonScintECorr->Fill(ChargeCorr,CL);

           hX->Fill(clX);
           hY->Fill(clY);
           hZ->Fill(clZ);
           hXY->Fill(clX,clY);

           // get u wire hits
           int nUWire = charge_cluster->GetNumUWireSignals();
           for (int UID = 0; UID < nUWire; UID++) {
              int chID = charge_cluster->GetUWireSignalChannelAt(UID);

              if (chID < 38) {hUHITPZ->Fill(chID);}
              else {hUHITNZ->Fill(chID);}
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
           fcsc = SumScint;
           feccl = ChargeCorr;
           fCL = CL;
           fAlphaEnergy = ChargeCorr + factor*SumScint;

           v->StackData(charge_cluster);

           nAlphas++;

           oTree->Fill();
        }
     }
  }

  /*TCanvas *c1 = new TCanvas("c1","Scintillation");
  hScintS->Draw();
  hScintB->Draw("same");
  hScintC->Draw("same");
  hScintA->Draw("same");

  TCanvas *c2 = new TCanvas("c2","ionization energy");
  hECorrS->Draw();
  hECorrB->Draw("same");
  hECorrC->Draw("same");
  hECorrA->Draw("same");

  TCanvas *c3 = new TCanvas("c3","Ion/Scint");
  hIonScintCorrS->Draw();
  hIonScintCorrB->Draw("same");
  hIonScintCorrC->Draw("same");
  hIonScintCorrA->Draw("same");

  TCanvas *c4 = new TCanvas("c4","Alpha spectrum");
  hAlphaSpectrumS->Draw();
  hAlphaSpectrumB->Draw("same");
  hAlphaSpectrumC->Draw("same");
  hAlphaSpectrumA->Draw("same");

  TCanvas *c5 = new TCanvas("c5","Ion/Scint vs energy");
  grB->Draw("A*");
  grC->Draw("*same");
  grA->Draw("*same");
  c5->SetLogy();

  TCanvas *c6 = new TCanvas("c6","X position");
  hX->Draw();

  TCanvas *c7 = new TCanvas("c7","y position");
  hY->Draw();

  TCanvas *c8 = new TCanvas("c8","Z position");
  hZ->Draw();

  TCanvas *c9 = new TCanvas("c9","X-Y position");
  hXY->Draw("colz");

  TCanvas *c10 = new TCanvas("c10","U wire hits (PZ)");
  hUHITPZ->Draw();

  TCanvas *c11 = new TCanvas("c11","U wire hits (NZ)");
  hUHITNZ->Draw();

  v->DrawStackedData();
  v->Save("Alpha3D.root");*/

  // write tree to file
  TFile *oFile = new TFile("../analysis/AlphaSpectroscopy_Result.root","RECREATE");
  oTree->Write();

  oFile->Close();
}

bool IsNoiseEvent(EXOEventData *WFED)
{
  // reset histograms
  hAPDSumWaveform->Reset();
  hAPDSumWaveformBLSub->Reset();
  hWireSumWaveform->Reset();
  hWireSumWaveformBLSub->Reset();
  hWireSumWaveformInv->Reset();
  hWireSumWaveformBLSubInv->Reset();

  // define noise threshold
  double ChargeSigmaThreshold = 4.0;
  double ChargeThreshold = 300;
  double ChargeThresholdInv = 300;
  double LightSigmaThreshold = 4.0;
  double LightThreshold = 250;

  // get waveform data and decompress them
  EXOWaveformData* wf_data = WFED->GetWaveformData();
  wf_data->Decompress();

  // get number of wavefrom sample
  int n_sample = wf_data->fNumSamples;

  // loop over all channels
  int *sumAPD = new int[n_sample];
  int *sumWire = new int[n_sample];
  for (int k = 0; k < n_sample; k++) {sumAPD[k] = 0; sumWire[k] = 0;}

  int nrSig = wf_data->GetNumWaveforms();
  for (int i = 0; i < nrSig; i++) {
     EXOWaveform* wf = wf_data->GetWaveform(i);
     int chID = wf->fChannel;

     if (chID >= 152) {
        for (int k = 0; k < n_sample; k++) {sumAPD[k] += wf->At(k);}
     }
     else {
        for (int k = 0; k < n_sample; k++) {sumWire[k] += wf->At(k);}
     }
  }

  // fill histograms
  for (int k = 0; k < n_sample; k++) {
     hAPDSumWaveform->SetBinContent(k+1,sumAPD[k]);
     hWireSumWaveform->SetBinContent(k+1,sumWire[k]);
     hWireSumWaveformInv->SetBinContent(k+1,-1 * sumWire[k]);
  }

  // define variables for baseline and standard deviation
  double bl_APD;
  double bl_Wire;
  double stdev_APD;
  double stdev_Wire;

  // subtract the baseline form the waveforms
  SubtractBaseline(sumAPD,n_sample,&bl_APD,&stdev_APD);
  SubtractBaseline(sumWire,n_sample,&bl_Wire,&stdev_Wire);

  // define threshold based on the number of sigmas
  double noiseThreshold_APD_TMP = bl_APD + LightSigmaThreshold * TMath::Sqrt(stdev_APD);
  int noiseThreshold_ADP = int(noiseThreshold_APD_TMP);

  double noiseThreshold_Wire_TMP = bl_Wire + ChargeSigmaThreshold * TMath::Sqrt(stdev_Wire);
  int noiseThreshold_Wire = int(noiseThreshold_Wire_TMP);

  // get signals above threshold
  for (int k = 0; k < n_sample; k++) {
     if (sumAPD[k] > LightThreshold) {hAPDSumWaveformBLSub->SetBinContent(k+1,sumAPD[k]);}
     else {hAPDSumWaveformBLSub->SetBinContent(k+1,0);}

     if (sumWire[k] > ChargeThreshold) {hWireSumWaveformBLSub->SetBinContent(k+1,sumWire[k]);}
     else {hWireSumWaveformBLSub->SetBinContent(k+1,0);}

     if (-1 * sumWire[k] > ChargeThresholdInv) {hWireSumWaveformBLSubInv->SetBinContent(k+1,-1*sumWire[k]);}
     else {hWireSumWaveformBLSubInv->SetBinContent(k+1,0);}
  }

  // find the peaks in the spectrum
  TSpectrum *specAPD = new TSpectrum();
  specAPD->Search(hAPDSumWaveformBLSub, 2, "goff");

  int nPeaksAPD = specAPD->GetNPeaks();
  float *peakPositionAPDX = specAPD->GetPositionX();
  float *peakPositionAPDY = specAPD->GetPositionY();

  TSpectrum *specWire = new TSpectrum();
  specWire->Search(hWireSumWaveformBLSub, 2, "goff");

  int nPeaksWire = specWire->GetNPeaks();
  float *peakPositionWireX = specWire->GetPositionX();
  float *peakPositionWireY = specWire->GetPositionY();

  TSpectrum *specWireInv = new TSpectrum();
  specWireInv->Search(hWireSumWaveformBLSubInv, 2, "goff");

  int nPeaksWireInv = specWireInv->GetNPeaks();
  float *peakPositionWireInvX = specWireInv->GetPositionX();
  float *peakPositionWireInvY = specWireInv->GetPositionY();

  // find matches of peak positions in APD and wire spectrum
  int nMatchesP = 0;
  for (int i = 0; i < nPeaksWire; i++) {
     for (int k = 0; k < nPeaksAPD; k++) {if (TMath::Abs(peakPositionAPDX[k]-peakPositionWireX[i]) > 10) {nMatchesP++; break;}}
  }

  int nMatchesN = 0;
  for (int i = 0; i < nPeaksWireInv; i++) {
     for (int k = 0; k < nPeaksAPD; k++) {if (TMath::Abs(peakPositionAPDX[k]-peakPositionWireInvX[i]) > 10) {nMatchesN++; break;}}
  }

  // clean up memory
  delete sumAPD;
  delete sumWire;
  delete specAPD;
  delete specWire;
  delete specWireInv;

  // apply cuts
  //if (nPeaksWire > 0 || nPeaksWireInv > 0) {return true;}
  if (nMatchesP > 0 || nMatchesN > 0) {return true;}

  return false;
}

void SubtractBaseline(int *DataSamples, int n_sample, double *bl, double *stdev)
{
  double blTMP = 0.0;
  double nSumSampleTMP = double(n_sample) / 6.0;
  int nSumSample = int(nSumSampleTMP);
  double *DataMean = new double[nSumSample];

  for (int i = 0; i < nSumSample; i++) {
     blTMP += double(DataSamples[i]);
     DataMean[i] = double(DataSamples[i]);
  }

  // mean of first nsample/6 samples
  if (nSumSample > 0) {blTMP /= nSumSample;}

  // standard deviation of first nsample/6 samples
  double stdevTMP = 0.0;
  for (int i = 0; i < nSumSample; i++) {
     stdevTMP += (DataMean[i] - blTMP)*(DataMean[i] - blTMP);
  }

  if (nSumSample > 0) {stdevTMP = TMath::Sqrt(stdevTMP / nSumSample);}

  for (int i = 0; i < n_sample; i++) {DataSamples[i] = DataSamples[i] - int(blTMP);}

  *bl = blTMP;
  *stdev = stdevTMP;

  // clean up memory
  delete DataMean;

  return;
}
