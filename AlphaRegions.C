void AlphaRegions()
{
  gROOT->SetStyle("Plain");

  gSystem->Load("libEXOUtilities");

  TFile *f = new TFile("../analysis/FirstDataPeriod/AlphaSpectroscopy_838_858_fidCut_10_180.root","READ");
  TTree *t = (TTree*)f->Get("t");

  double fX;
  double fY;
  double fZ;
  double fCL;
  double feccl;
  double fAlphaEnergy;

  t->SetBranchAddress("fX",&fX);
  t->SetBranchAddress("fY",&fY);
  t->SetBranchAddress("fZ",&fZ);
  t->SetBranchAddress("fCL",&fCL);
  t->SetBranchAddress("feccl",&feccl);
  t->SetBranchAddress("fAlphaEnergy",&fAlphaEnergy);

  TH1F *hAlphaSpectrum = new TH1F("hAlphaSpectrum","Alpha energy spectrum",200,0,10000);
  TH1F *hAlphaRegion1 = new TH1F("hAlpha1Coinc","Energy of first alpha in coincidence",200,0,10000);
  TH1F *hAlphaRegion2 = new TH1F("hAlpha2Coinc","Energy of second alpha in coincidence",200,0,10000);

  TH1F *hX1 = new TH1F("hX1","X positions",50,-200,200);
  TH1F *hY1 = new TH1F("hY1","Y positions",50,-200,200);
  TH1F *hZ1 = new TH1F("hZ1","Z positions",50,-200,200);
  TH2F *hXY1 = new TH2F("hXY1","X-Y positions",50,-200,200,50,-200,200);

  TH1F *hX2 = new TH1F("hX2","X positions",50,-200,200);
  TH1F *hY2 = new TH1F("hY2","Y positions",50,-200,200);
  TH1F *hZ2 = new TH1F("hZ2","Z positions",50,-200,200);
  TH2F *hXY2 = new TH2F("hXY2","X-Y positions",50,-200,200,50,-200,200);

  hAlphaRegion1->SetLineColor(kRed);
  hAlphaRegion1->SetFillColor(kRed);
  hAlphaRegion1->SetFillStyle(3004);

  hAlphaRegion2->SetLineColor(kBlue);
  hAlphaRegion2->SetFillColor(kBlue);
  hAlphaRegion2->SetFillStyle(3005);

  hX1->SetLineColor(kRed);
  hX1->SetFillColor(kRed);
  hX1->SetFillStyle(3004);

  hY1->SetLineColor(kRed);
  hY1->SetFillColor(kRed);
  hY1->SetFillStyle(3004);

  hZ1->SetLineColor(kRed);
  hZ1->SetFillColor(kRed);
  hZ1->SetFillStyle(3004);

  hXY1->SetMarkerColor(kRed);
  hXY1->SetMarkerStyle(20);
  hXY1->SetMarkerSize(0.5);

  hX2->SetLineColor(kBlue);
  hX2->SetFillColor(kBlue);
  hX2->SetFillStyle(3005);

  hY2->SetLineColor(kBlue);
  hY2->SetFillColor(kBlue);
  hY2->SetFillStyle(3005);

  hZ2->SetLineColor(kBlue);
  hZ2->SetFillColor(kBlue);
  hZ2->SetFillStyle(3005);

  hXY2->SetMarkerColor(kBlue);
  hXY2->SetMarkerStyle(20);
  hXY2->SetMarkerSize(0.5);

  TGraph *gr = new TGraph(1000);
  gr->GetXaxis()->SetTitle("energy [keV]");
  gr->GetXaxis()->SetRangeUser(0,3000);
  gr->GetYaxis()->SetTitle("ion/scint [a.u.]");
  gr->GetYaxis()->SetRangeUser(0,1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);

  TGraph *gr2 = new TGraph(1000);
  gr2->GetXaxis()->SetTitle("energy [keV]");
  gr2->GetXaxis()->SetRangeUser(0,3000);
  gr2->GetYaxis()->SetTitle("ion/scint [a.u.]");
  gr2->GetYaxis()->SetRangeUser(0,1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.5);
  gr2->SetMarkerColor(kRed);

  TGraph *gr3 = new TGraph(1000);
  gr3->GetXaxis()->SetTitle("energy [keV]");
  gr3->GetXaxis()->SetRangeUser(0,3000);
  gr3->GetYaxis()->SetTitle("ion/scint [a.u.]");
  gr3->GetYaxis()->SetRangeUser(0,1);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.5);
  gr3->SetMarkerColor(kBlue);

  int nentries = t->GetEntries();
  cout << "Number of entries: " << nentries << endl;

  // loop over all alpha events
  int nPoints1 = 0;
  int nPoints2 = 0;
  for (int i = 0; i < nentries; i++) {
     t->GetEntry(i);

     hAlphaSpectrum->Fill(fAlphaEnergy);
     gr->SetPoint(i,feccl,fCL);

     if (fCL < 0.000075*feccl) {
        hAlphaRegion1->Fill(fAlphaEnergy);
        hX1->Fill(fX);
        hY1->Fill(fY);
        hZ1->Fill(fZ);
        hXY1->Fill(fX,fY);
        gr2->SetPoint(nPoints1,feccl,fCL);
        nPoints1++;
     }
     else {
        hAlphaRegion2->Fill(fAlphaEnergy);
        hX2->Fill(fX);
        hY2->Fill(fY);
        hZ2->Fill(fZ);
        hXY2->Fill(fX,fY);
        gr3->SetPoint(nPoints2,feccl,fCL);
        nPoints2++;
     }
  }

  TCanvas *c1 = new TCanvas("c1","Ion/Scint vs ionization");
  gr->Draw("AP");

  TCanvas *c2 = new TCanvas("c2","Ion/Scint vs ionization");
  gr2->Draw("AP");
  gr3->Draw("Psame");

  TCanvas *c3 = new TCanvas("c3","Energy Spectrum");
  hAlphaSpectrum->Draw();
  hAlpha1Coinc->Draw("same");
  hAlpha2Coinc->Draw("same");

  TCanvas *c4 = new TCanvas("c4","X positions");
  hX1->Draw();
  hX2->Draw("same");

  TCanvas *c5 = new TCanvas("c5","Y positions");
  hY1->Draw();
  hY2->Draw("same");

  TCanvas *c6 = new TCanvas("c6","Z positions");
  hZ1->Draw();
  hZ2->Draw("same");

  TCanvas *c7 = new TCanvas("c7","X-Y positions");
  hXY1->Draw();
  hXY2->Draw("same");

  return;
}

