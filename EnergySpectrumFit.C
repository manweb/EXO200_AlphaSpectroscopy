double fitEnergy(double *x, double *par)
{
  //if (par[6] < 0.2*par[3]) {par[6] = 0.2*par[3];}
  //if (par[6] > 0.8*par[3]) {par[6] = 0.8*par[3];}
  //double val = par[0]*TMath::Gaus(x[0],par[1],par[2]) + par[3]*TMath::Gaus(x[0],par[4],par[5]) + par[6]*TMath::Gaus(x[0],par[7],par[8]) + par[9]*TMath::Gaus(x[0],par[10],par[11]);
  double val = par[0]*TMath::Gaus(x[0],par[1],par[2]) + par[3]*TMath::Gaus(x[0],par[4],par[5]) + par[6]*TMath::Gaus(x[0],par[7],par[5]) + par[9]*TMath::Gaus(x[0],par[10],par[11]);

  return val;
}

void EnergySpectrumFit()
{
  TChain *t = new TChain("t");
  t->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1882_physics.root");

  TH1F *hNCL_EQ1_fid_energy = new TH1F("hNCL_EQ1_fid_energy","ncl = 1 (fiducial volume, converted)",100,0,10000);
  hNCL_EQ1_fid_energy->GetXaxis()->SetTitle("energy [keV]");
  hNCL_EQ1_fid_energy->GetYaxis()->SetTitle("entries/100 keV");

  // 0.255
  //t->Draw("0.2496*fcsc+feccl>>hNCL_EQ1_fid_energy","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -172 && fZ < -20) || (fZ > 20 && fZ < 172))");

  //double RunTime = 23.08;
  //double RunTime = 40.79;
  double RunTime = 31.4025;

  int nentries = t->GetEntries();

  int fncl;
  double fx;
  double fy;
  double fz;
  double feccl;
  double fcsc;

  t->SetBranchAddress("fncl",&fncl);
  t->SetBranchAddress("feccl",&feccl);
  t->SetBranchAddress("fcsc",&fcsc);
  t->SetBranchAddress("fX",&fx);
  t->SetBranchAddress("fY",&fy);
  t->SetBranchAddress("fZ",&fz);

  double a1 = -11.16;
  double a2 = 29.46;
  double b1 = 21112.1;
  double b2 = 21200.6;

  double z0 = 96.0;
  double csc_middle = ((a2 - a1)*z0 + b2 + b1) / 2.0;
  double corr1 = csc_middle / (-1.0*a1*z0 + b1);
  double corr2  = csc_middle / (a2*z0 + b2);

  double b = 0.0;
  double c = 0.0;
  double csc_corr = 0.0;
  for (int i = 0; i < nentries; i++) {
     t->GetEntry(i);

     if (TMath::Sqrt(fx*fx + fy*fy) > 163) {continue;}
     if (fz < -172 || fz > 172) {continue;}
     if (fz > -20 && fz < 20) {continue;}

     if (fz < 0) {
        b = fcsc - a1 * fz;
        c = fcsc / (-1.0*a1*z0 + b);
        csc_corr = fcsc / c * corr1;
     }
     else {
        b = fcsc - a2 * fz;
        c = fcsc / (a2*z0 + b);
        csc_corr = fcsc / c * corr2;
     }

     hNCL_EQ1_fid_energy->Fill(0.2423*csc_corr + feccl);
  }

  // fit to energy spectrum
  TF1 *fitE = new TF1("fitE",fitEnergy,4000,10000,12);
  fitE->SetLineColor(kBlue);
  fitE->SetLineWidth(1);

  double params[12] = {2,4270,500,30,5590,500,30,6115,500,8,7883,500};
  fitE->SetParameters(params);

  /*fitE->FixParameter(1,4270);
  fitE->FixParameter(4,5590);
  fitE->FixParameter(6,6115);
  fitE->FixParameter(8,5407);
  fitE->FixParameter(11,7883)*/

  fitE->SetParLimits(0,0,10);
  fitE->SetParLimits(1,4300,4500);
  fitE->SetParLimits(2,0,250);
  fitE->SetParLimits(3,0,200);
  fitE->SetParLimits(4,5400,5700);
  fitE->SetParLimits(5,0,1000);
  fitE->SetParLimits(6,0,200);
  fitE->SetParLimits(7,5900,6300);
  fitE->SetParLimits(8,0,1000);
  fitE->SetParLimits(9,0,50);
  fitE->SetParLimits(10,7700,8000);
  fitE->SetParLimits(11,0,1000);

  hNCL_EQ1_fid_energy->Fit("fitE","rnL");

  double par[12];
  fitE->GetParameters(par);

  //double parErrors[13];
  double *parErrors = fitE->GetParErrors();

  TF1 *rn222 = new TF1("rn222","gaus",3000,10000);
  TF1 *po218 = new TF1("po218","gaus",3000,10000);
  //TF1 *po210 = new TF1("po210","gaus",3000,10000);
  TF1 *u238 = new TF1("u238","gaus",3000,10000);
  TF1 *po214 = new TF1("po214","gaus",3000,10000);

  rn222->SetLineColor(kRed);
  rn222->SetLineWidth(1);
  rn222->SetLineStyle(2);

  po218->SetLineColor(kGreen);
  po218->SetLineWidth(1);
  po218->SetLineStyle(2);

  /*po210->SetLineColor(kGray+2);
  po210->SetLineWidth(1);
  po210->SetLineStyle(2);*/

  rn222->SetParameter(0,par[3]);
  rn222->SetParameter(1,par[4]);
  rn222->SetParameter(2,par[5]);

  po218->SetParameter(0,par[6]);
  po218->SetParameter(1,par[7]);
  po218->SetParameter(2,par[5]);

  /*po210->SetParameter(0,par[7]);
  po210->SetParameter(1,par[8]);
  po210->SetParameter(2,par[9]);*/

  u238->SetParameter(0,par[0]);
  u238->SetParameter(1,par[1]);
  u238->SetParameter(2,par[2]);

  po214->SetParameter(0,par[9]);
  po214->SetParameter(1,par[10]);
  po214->SetParameter(2,par[11]);

  u238->SetLineColor(kYellow);
  u238->SetLineWidth(1);
  u238->SetLineStyle(2);

  po214->SetLineColor(kBlack);
  po214->SetLineWidth(1);
  po214->SetLineStyle(2);

  double IntegralError;
  double Inegral = hNCL_EQ1_fid_energy->IntegralAndError(38,100,IntegralError);
  cout << Inegral << " +- " << IntegralError << endl;

  cout << rn222->Integral(3000,10000)/100.0/RunTime << endl;
  cout << po218->Integral(3000,10000)/100.0/RunTime << endl;
  //cout << po210->Integral(3000,10000)/100.0/RunTime << endl;
  cout << u238->Integral(3000,10000)/100.0/RunTime << endl;
  cout << po214->Integral(3000,10000)/100.0/RunTime << endl;

  TLegend *l = new TLegend(0.0,0.0,0.2,0.2);

  char rn222Label[50];
  char po218Label[50];
  //char po210Label[50];
  char u238Label[50];
  char po214Label[50];
  /*sprintf(rn222Label,"Rn222 (%.1f +- %.1f)",rn222->Integral(3000,10000)/100.0/RunTime,rn222->IntegralError(3000,10000)/100.0/RunTime);
  sprintf(po218Label,"Po218 (%.1f +- %.1f)",po218->Integral(3000,10000)/100.0/RunTime,po218->IntegralError(3000,10000)/100.0/RunTime);
  sprintf(po210Label,"Po210 (%.1f +- %.1f)",po210->Integral(3000,10000)/100.0/RunTime,po210->IntegralError(3000,10000)/100.0/RunTime);
  sprintf(u238Label,"U238 (%.1f +- %.1f)",u238->Integral(3000,10000)/100.0/RunTime,u238->IntegralError(3000,10000)/100.0/RunTime);
  sprintf(po214Label,"Po214 (%.1f +- %.1f)",po214->Integral(3000,10000)/100.0/RunTime,po214->IntegralError(3000,10000)/100.0/RunTime);*/
  sprintf(rn222Label,"Rn222 (%.1f +- %.1f)",par[3]*par[5]*2.5066/100.0/RunTime,TMath::Sqrt(par[5]*par[5]*parErrors[3]*parErrors[3] + par[3]*par[3]*parErrors[5]*parErrors[5])*1.77245/100.0/RunTime);
  sprintf(po218Label,"Po218 (%.1f +- %.1f)",par[6]*par[5]*2.5066/100.0/RunTime,TMath::Sqrt(par[5]*par[5]*parErrors[6]*parErrors[6] + par[6]*par[6]*parErrors[5]*parErrors[5])*1.77245/100.0/RunTime);
  //sprintf(po210Label,"Po210 (%.1f +- %.1f)",par[7]*par[9]*2.5066/100.0/RunTime,TMath::Sqrt(par[9]*par[9]*parErrors[7]*parErrors[7] + par[7]*par[7]*parErrors[9]*parErrors[9])*1.77245/100.0/RunTime);
  sprintf(u238Label,"U238 (%.1f +- %.1f)",par[0]*par[2]*2.5066/100.0/RunTime,TMath::Sqrt(par[2]*par[2]*parErrors[0]*parErrors[0] + par[0]*par[0]*parErrors[2]*parErrors[2])*1.77245/100.0/RunTime);
  sprintf(po214Label,"Po214 (%.1f +- %.1f)",par[9]*par[11]*2.5066/100.0/RunTime,TMath::Sqrt(par[9]*par[11]*parErrors[9]*parErrors[11] + par[9]*par[9]*parErrors[11]*parErrors[11])*1.77245/100.0/RunTime);

  l->AddEntry(rn222,rn222Label);
  l->AddEntry(po218,po218Label);
  //l->AddEntry(po210,po210Label);
  l->AddEntry(u238,u238Label);
  l->AddEntry(po214,po214Label);

  TCanvas *c1 = new TCanvas();
  hNCL_EQ1_fid_energy->Draw("EP");
  fitE->Draw("same");
  rn222->Draw("same");
  po218->Draw("same");
  //po210->Draw("same");
  u238->Draw("same");
  po214->Draw("same");
  l->Draw("same");

  return;
}
