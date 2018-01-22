double fitFunction(double *x, double *par)
{
   if (x[0] > 16000 && x[0] < 23000) {
      TF1::RejectPoint();
      return 0;
   }

   return TMath::Exp(par[0] + par[1]*x[0]);
}

void ChargeCollectionEfficiency()
{
  TChain *t = new TChain("t");
  t->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1782.root");

  // 214Po data
  TChain *t2 = new TChain("t");
  t2->Add("../../BiPo/analysis/SecondDataPeriod/New/bga_all.root");

  // define histograms
  TH1F *hNCL_EQ1 = new TH1F("hNCL_EQ1","ncl = 1",100,0,40000);
  TH1F *hNCL_LT1 = new TH1F("hNCL_LT1","ncl <= 1",100,0,40000);
  TH1F *hNCL_EQ0 = new TH1F("hNCL_EQ0","ncl == 0",100,0,40000);
  TH1F *hNCL_EQ1_fid = new TH1F("hNCL_EQ1_fid","ncl = 1 (fiducial volume)",100,0,40000);
  TH1F *h214Po_LT1 = new TH1F("h214Po_LT1","214Po spectrum ncl <= 1",100,0,40000);
  TH1F *h214Po_EQ1 = new TH1F("h214Po_EQ1","214Po spectrum ncl == 1",100,0,40000);
  TH1F *h214Po_fid = new TH1F("h214Po_fid","214Po spectrum (fiducial)",100,0,40000);

  hNCL_EQ1->SetLineColor(kRed);
  //hNCL_EQ1->SetLineStyle(2);

  hNCL_EQ0->SetLineColor(kBlue);
  //hNCL_EQ0->SetLineStyle(2);

  hNCL_EQ1_fid->SetLineColor(kBlue);
  //hNCL_EQ1_fid->SetLineStyle(2);

  t->Draw("fcsc>>hNCL_EQ1","fncl == 1");
  t->Draw("fcsc>>hNCL_LT1");
  t->Draw("fcsc>>hNCL_EQ0","fncl == 0");
  t->Draw("fcsc>>hNCL_EQ1_fid","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -172 && fZ < -20) || (fZ > 20 && fZ < 172))");

  t2->Draw("fcsc2>>h214Po_LT1");
  t2->Draw("fcsc2>>h214Po_EQ1","fepcl2 > 0");
  t2->Draw("fcsc2>>h214Po_fid","TMath::Sqrt(fx2*fx2 + fy2*fy2) < 163 && ((fz2 > -172 && fz2 < -20) || (fz2 > 20 && fz2 < 172))");

  // subtract 214Po spectrum
  //hNCL_EQ1->Add(h214Po_EQ1,-1);
  //hNCL_LT1->Add(h214Po_LT1,-1);
  //hNCL_EQ1_fid->Add(h214Po_fid,-1);

  TF1 *fit = new TF1("fit",fitFunction,3500,30000,2);

  fit->SetLineColor(kGreen);
  fit->SetLineWidth(2);

  hNCL_EQ0->Fit("fit","rn");

  double par[2];
  fit->GetParameters(par);

  TF1 *background = new TF1("background","expo",3500,30000);
  background->SetParameters(par);

  background->SetLineColor(kGreen);
  background->SetLineWidth(2);

  // create clone of sum spectrum and subtract background
  TH1F *hNCL_LT1_copy = (TH1F*)hNCL_LT1->Clone("hNCL_LT1_copy");
  //hNCL_LT1_copy->Add(fit,-1);
  hNCL_LT1_copy->Add(background,-1);

  // create clone of diff spectrum and divide by ncl == 1 spectrum
  TH1F *hNCL_LT1_div = (TH1F*)hNCL_LT1_copy->Clone("hNCL_LT1_div");
  //hNCL_LT1_div->Divide(hNCL_EQ1);
  hNCL_LT1_div->Add(hNCL_EQ1,-1);

  // create copy of fiducial spectrum and multiply by div spectrum
  TH1F *hNCL_EQ1_fid_copy = (TH1F*)hNCL_EQ1_fid->Clone("hNCL_EQ1_fid_copy");
  //hNCL_EQ1_fid_copy->Multiply(hNCL_LT1_div);
  hNCL_EQ1_fid_copy->Add(hNCL_LT1_div);

  hNCL_EQ1_fid_copy->SetLineColor(kGray);
  hNCL_EQ1_fid_copy->SetFillColor(kGray);
  //hNCL_EQ1_fid_copy->SetLineStyle(2);

  // create clone of ncl == 0 spectrum and subtract background
  TH1F *hNCL_EQ0_copy = (TH1F*)hNCL_EQ0->Clone("hNCL_EQ0_copy");
  //hNCL_EQ0_copy->Add(fit,-1);
  hNCL_EQ0_copy->Add(background,-1);

  // zero out negative bins
  for (int i = 0; i < 100; i++) {
     if (i < 38) {hNCL_EQ1_fid_copy->SetBinContent(i+1,0);}
     double BinContent = hNCL_EQ1_fid_copy->GetBinContent(i+1);
     if (BinContent < 0) {hNCL_EQ1_fid_copy->SetBinContent(i+1,0);}
  }

  double IntegralError;
  double Integral = hNCL_EQ1_fid_copy->IntegralAndError(1,100,IntegralError);

  cout << hNCL_EQ1_fid_copy->Integral() << endl;
  cout << Integral << " +- " << IntegralError << endl;
  cout << (hNCL_EQ1_fid->Integral())/hNCL_EQ1_fid_copy->Integral() << endl;

  TLegend *l = new TLegend(0,0,0.2,0.2);
  l->AddEntry(hNCL_LT1_copy,"ncl <= 1");
  l->AddEntry(hNCL_EQ1,"ncl == 1");
  l->AddEntry(hNCL_EQ1_fid,"ncl == 1 (fiducial)");
  l->AddEntry(hNCL_EQ1_fid_copy,"ncl <= 1 (fiducial)");

  TLegend *l2 = new TLegend(0,0,0.2,0.2);
  l2->AddEntry(hNCL_LT1,"ncl <= 1");
  l2->AddEntry(hNCL_EQ1,"ncl == 1");
  l2->AddEntry(hNCL_EQ0,"ncl == 0");
  l2->AddEntry(fit,"fit to background");

  TCanvas *c1 = new TCanvas();
  hNCL_LT1->Draw();
  hNCL_EQ1->Draw("same");
  hNCL_EQ0->Draw("same");
  l2->Draw("same");
  background->Draw("same");

  TCanvas *c2 = new TCanvas();
  hNCL_LT1_copy->Draw();
  hNCL_EQ1->Draw("same");
  hNCL_EQ1_fid_copy->Draw("same");
  hNCL_EQ1_fid->Draw("same");
  l->Draw("same");

  TCanvas *c3 = new TCanvas();
  hNCL_EQ0_copy->Draw();

  TCanvas *c4 = new TCanvas();
  hNCL_EQ1_fid_copy->Draw();

  TCanvas *c5 = new TCanvas();
  hNCL_LT1->Draw();
  h214Po_LT1->Draw("same");

  return;
}
