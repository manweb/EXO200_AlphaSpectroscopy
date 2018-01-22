void ChargeLight_vs_eccl()
{
  gROOT->SetStyle("Plain");

  //TFile *f = new TFile("../analysis/FirstDataPeriod/AlphaSpectroscopy_838_858.root","READ");
  TFile *f = new TFile("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1782.root","READ");

  TTree *t = (TTree*)f->Get("t");

  TH2F *h1 = new TH2F("h1","Alpha events",100,0,3000,100,0,2);

  double feccl;
  double fcsc;
  double fX;
  double fY;
  double fZ;
  double fCL;

  t->SetBranchAddress("feccl",&feccl);
  t->SetBranchAddress("fcscCorr",&fcsc);
  t->SetBranchAddress("fX",&fX);
  t->SetBranchAddress("fY",&fY);
  t->SetBranchAddress("fZ",&fZ);
  t->SetBranchAddress("fCL",&fCL);

  int nentries = t->GetEntries();

  TGraph *gr = new TGraph(nentries);

  gr->GetXaxis()->SetTitle("energy [keV]");
  gr->GetXaxis()->SetRangeUser(0,3000);
  gr->GetYaxis()->SetTitle("ion/scint [a.u.]");
  gr->GetYaxis()->SetRangeUser(0,2);
  gr->SetMarkerStyle(6);
  gr->SetMarkerSize(1.0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);

  h1->GetXaxis()->SetTitle("ionization energy [keV]");
  h1->GetYaxis()->SetTitle("ion/scint [keV/photon counts]");

  for (int i = 0; i < nentries; i++) {
     t->GetEntry(i);

     double R = TMath::Sqrt(fX*fX + fY*fY);
     if (R > 163) {continue;}
     if (!(fZ > -180 && fZ < -20) && !(fZ > 20 && fZ < 180)) {continue;}

     gr->SetPoint(i,feccl,fCL);
  }

  // lines and area for exclusion region
  TLine *lH = new TLine(0,0.05,1600,0.05);
  TLine *lV = new TLine(600,0,600,1.7);
  lH->SetLineStyle(2);
  lV->SetLineStyle(2);

  int n = 5;
  double x[5] = {600,1600,1600,600,600};
  double y[5] = {0.05,0.05,1.7,1.7,0.05};
  TPolyLine *p = new TPolyLine(n,x,y);
  p->SetFillStyle(3005);
  p->SetFillColor(kGray+1);

  TCanvas *c1 = new TCanvas();
  gr->Draw("AP");
  lH->Draw("same");
  lV->Draw("same");
  p->Draw("fsame");
  c1->SetLogy();
}

