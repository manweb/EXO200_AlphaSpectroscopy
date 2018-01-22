void Charge_vs_Light()
{
  gROOT->SetStyle("Plain");

  //TFile *f = new TFile("../analysis/FirstDataPeriod/AlphaSpectroscopy_838_858.root","READ");
  TFile *f = new TFile("../analysis/SecondDataPeriod/AlphaSpectroscopy_1333.root","READ");

  TTree *t = (TTree*)f->Get("t");

  TH1F *hZ = new TH1F("hZ","Z position of alpha events",80,-200,200);
  TH2F *h1 = new TH2F("h1","Alpha events",100,0,500,100,0,40000);
  TH2F *h2 = new TH2F("h2","Alpha events",100,0,500,100,0,40000);
  TH2F *h3 = new TH2F("h3","Scintillation vs Z",80,-200,200,100,0,30000);

  t->Draw("fZ>>hZ");
  t->Draw("fcsc:feccl>>h1","fZ>-10 && fZ<10");
  t->Draw("fcsc:feccl>>h2","fZ>-150 && fZ<=-10 || fZ>=10 && fZ<150");
  t->Draw("fcsc:fZ>>h3");

  hZ->GetXaxis()->SetTitle("z [mm]");

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);

  h1->GetXaxis()->SetTitle("ionization energy [keV]");
  h1->GetYaxis()->SetTitle("scintillaion [photon counts]");

  h2->SetLineColor(kRed);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.5);

  h3->GetXaxis()->SetTitle("z [mm]");
  h3->GetYaxis()->SetTitle("scintillation [photon counts]");

  TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
  l->AddEntry(h1,"Events on the cathode");
  l->AddEntry(h2,"Events in the bulk");

  TCanvas *c1 = new TCanvas();
  h1->Draw();
  h2->Draw("same");
  l->Draw("same");

  TCanvas *c2 = new TCanvas();
  hZ->Draw();

  TCanvas *c2 = new TCanvas();
  h3->Draw();
}

