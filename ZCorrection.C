void ZCorrection()
{
  TChain *t1 = new TChain("t");
  TChain *t2 = new TChain("t");

  t1->Add("../../BiPo/analysis/SecondDataPeriod/Reprocessed_110817/bga_all_1569_1882_physics.root");
  t2->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1882_physics.root");

  int nentries1 = t1->GetEntries();
  int nentries2 = t2->GetEntries();

  double fx1;
  double fy1;
  double fz1;
  double fx2;
  double fy2;
  double fz2;
  double fcsc1;
  double fcsc2;

  t1->SetBranchAddress("fcsc2",&fcsc1);
  t1->SetBranchAddress("fx2",&fx1);
  t1->SetBranchAddress("fy2",&fy1);
  t1->SetBranchAddress("fz2",&fz1);

  t2->SetBranchAddress("fcsc",&fcsc2);
  t2->SetBranchAddress("fX",&fx2);
  t2->SetBranchAddress("fY",&fy2);
  t2->SetBranchAddress("fZ",&fz2);

  TGraph *gr1 = new TGraph(nentries1);
  TGraph *gr2 = new TGraph(nentries2);

  TGraph *gr3 = new TGraph(nentries2);
  TGraph *gr4 = new TGraph(nentries2);

  gr1->SetMarkerStyle(4);
  gr2->SetMarkerStyle(6);
  gr1->SetMarkerColor(kRed);

  gr3->SetMarkerStyle(6);
  gr4->SetMarkerStyle(6);

  TH1F *h1 = new TH1F("h1","No correction",100,0,40000);
  TH1F *h2 = new TH1F("h2","Corrected",100,0,40000);

  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  for (int i = 0; i < nentries1; i++) {
     t1->GetEntry(i);

     if (TMath::Sqrt(fx1*fx1 + fy1*fy1) > 163) {continue;}
     if (fz1 < -172 || fz1 > 172) {continue;}
     if (fz1 > -20 && fz1 < 20) {continue;}

     gr1->SetPoint(n1,fz1,fcsc1);
     n1++;
  }

  for (int i = 0; i < nentries2; i++) {
     t2->GetEntry(i);

     if (TMath::Sqrt(fx2*fx2 + fy2*fy2) > 163) {continue;}
     if (fz2 < -172 || fz2 > 172) {continue;}
     if (fz2 > -20 && fz2 < 20) {continue;}

     gr2->SetPoint(n2,fz2,fcsc2);
     h1->Fill(fcsc2);

     if (fz2 < 0) {
        if (fcsc2 > -3000.0/158.0*fz2 + 14772.2 && fcsc2 < -3000.0/158.0*fz2 + 24772.2) {gr3->SetPoint(n3,fz2,fcsc2); n3++;}
     }
     else {
        if (fcsc2 > 8000.0/158.0*fz2 + 14696.2 && fcsc2 < 8000.0/158.0*fz2 + 23392.4) {gr3->SetPoint(n3,fz2,fcsc2); n3++;}
     }

     n2++;
  }

  TF1 *fit1 = new TF1("fit1","[0]*x+[1]",-180,-20);
  TF1 *fit2 = new TF1("fit2","[0]*x+[1]",20,180);

  fit1->SetLineWidth(1);
  fit1->SetLineColor(kRed);

  fit2->SetLineWidth(1);
  fit2->SetLineColor(kRed);

  gr3->Fit("fit1","rn");
  gr3->Fit("fit2","rn");

  double a1;
  double a2;
  double b1;
  double b2;

  a1 = fit1->GetParameter(0);
  a2 = fit2->GetParameter(0);
  b1 = fit1->GetParameter(1);
  b2 = fit2->GetParameter(1);

  cout << "a1=" << a1 << endl;
  cout << "a2=" << a2 << endl;
  cout << "b1=" << b1 << endl;
  cout << "b2=" << b2 << endl;

  double z0 = 96.0;
  double csc_middle = ((a2 - a1)*z0 + b2 + b1) / 2.0;
  double corr1 = csc_middle / (-1.0*a1*z0 + b1);
  double corr2  = csc_middle / (a2*z0 + b2);

  double b = 0.0;
  double c = 0.0;
  double csc_corr = 0.0;
  int n2_corr = 0.0;
  for (int i = 0; i < nentries2; i++) {
     t2->GetEntry(i);

     if (TMath::Sqrt(fx2*fx2 + fy2*fy2) > 163) {continue;}
     if (fz2 < -172 || fz2 > 172) {continue;}
     if (fz2 > -20 && fz2 < 20) {continue;}

     if (fz2 < 0) {
        b = fcsc2 - a1 * fz2;
        c = fcsc2 / (-1.0*a1*z0 + b);
        csc_corr = fcsc2 / c * corr1;
     }
     else {
        b = fcsc2 - a2 * fz2;
        c = fcsc2 / (a2*z0 + b);
        csc_corr = fcsc2 / c * corr2;
     }

     gr4->SetPoint(n2_corr,fz2,csc_corr);
     h2->Fill(csc_corr);

     n2_corr++;
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);

  TCanvas *c1 = new TCanvas();
  mg->Draw("AP");

  TCanvas *c2 = new TCanvas();
  gr3->Draw("AP");
  fit1->Draw("same");
  fit2->Draw("same");

  TCanvas *c3 = new TCanvas();
  gr4->Draw("AP");

  TCanvas *c4 = new TCanvas();
  c4->Divide(1,2);
  c4->cd(1);
  h1->Draw();
  c4->cd(2);
  h2->Draw();

  return;
}
