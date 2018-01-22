// define histograms
TH1F *hNCL_EQ1 = new TH1F("hNCL_EQ1","ncl = 1",100,0,40000);
TH1F *hNCL_LT1 = new TH1F("hNCL_LT1","ncl <= 1",100,0,40000);
TH1F *hNCL_EQ0 = new TH1F("hNCL_EQ0","ncl == 0",100,0,40000);
TH1F *hNCL_EQ1_fid = new TH1F("hNCL_EQ1_fid","ncl = 1 (fiducial volume)",100,0,40000);
TH1F *h214Po_LT1 = new TH1F("h214Po_LT1","214Po spectrum ncl <= 1",100,0,40000);
TH1F *h214Po_EQ1 = new TH1F("h214Po_EQ1","214Po spectrum ncl == 1",100,0,40000);
TH1F *h214Po_fid = new TH1F("h214Po_fid","Scintillation of 214Po",100,0,40000);
TH1F *h214Po_eccl = new TH1F("h214Po_eccl","Charge of 214Po",50,0,500);
TH1F *h214Po_epcl = new TH1F("h214Po_epcl","Purity corrected charge of 214Po",50,0,500);
TH1F *h214Po_CL = new TH1F("h214Po_CL","Ion/Scint of 214Po",50,0,1000);

TH1F *hCSC1 = new TH1F("hCSC1","Scintillation of 222Rn",100,0,40000);
TH1F *hCSC2 = new TH1F("hCSC2","Scintillation of 218Po",100,0,40000);
TH1F *hECCL1 = new TH1F("hECCL1","Charge first alpha",100,0,500);
TH1F *hEPCL1 = new TH1F("hEPCL1","Purity corrected charge first alpha",100,0,500);
TH1F *hECCL2 = new TH1F("hECCL2","Charge second alpha",100,0,500);
TH1F *hEPCL2 = new TH1F("hEPCL2","Purity corrected charge second alpha",100,0,500);
TH1F *hCL1 = new TH1F("hCL1","Ion/Scint first alpha",100,0,1000);
TH1F *hCL2 = new TH1F("hCL2","Ion/Scint second alpha",100,0,1000);

TH1F *hAlphaDT = new TH1F("hAlphaDT","Alpha coincidence dt",100,0,300000000);
TH2F *hXY = new TH2F("hXY","Position of coincidence alphas",40,-200,200,40,-200,200);

TChain *t = new TChain("t");

void GetConcidences();

double fitFunction(double *x, double *par)
{
   if (x[0] > 16000 && x[0] < 23000) {
      TF1::RejectPoint();
      return 0;
   }

   return TMath::Exp(par[0] + par[1]*x[0]);
}

double ErrFunction(double *x, double *par)
{
  //return par[0]*(1 + TMath::Erf((x[0]-13581)/5754));
  //return par[0]*(1 + TMath::Erf((x[0]-14518)/2631));
  return par[0]*(1 + TMath::Erf((x[0]-par[1])/par[2]));
}

double ErrFunctionCharge(double *x, double *par)
{
  return par[0]*(1 + TMath::Erf((x[0]-57.5)/15.81));
}

double LineFitError(double *x, double *par)
{
  double sdx1 = par[5] * par[1];
  double sdx2 = par[6] * par[1];
  double sdx3 = par[7] * par[1];

  double dy1 = par[8];
  double dy2 = par[9];
  double dy3 = par[10];

  double sy1 = TMath::Sqrt(sdx1*sdx1 + dy1*dy1);
  double sy2 = TMath::Sqrt(sdx2*sdx2 + dy2*dy2);
  double sy3 = TMath::Sqrt(sdx3*sdx3 + dy3*dy3);

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];

  double bracketX = 1.0/3.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3)); // [x]
  double bracketX2 = 1.0/3.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3)); // [x2]
  double bracket1 = 1.0/3.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/3.0 / bracket1;
  //double sa2 = 1.0/3.0 * bracketX2 / (bracketX2*bracket1 - bracketX*bracketX);
  double sb2 = 1.0/3.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  //cout << sa2 << "  " << sb2 << "  " << cov << endl;

  //cout <<x[0]<<"  "<<sb2<<"  "<<sa2<<"  "<<cov<<"  "<< (x[0]+cov)*(x[0]+cov)*sb2 + 2*(x[0]+cov)*cov + sa2 << endl;
  //return par[0] + par[1]*x[0] + par[2]*TMath::Sqrt((x[0]+par[4])*(x[0]+par[4])*par[3] + 2*(x[0]+par[4])*par[4] + par[5]);
  //return par[0] + par[1]*x[0] + par[11]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + 2*(x[0]+cov)*cov + sa2);
  return par[0] + par[1]*x[0] + par[11]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
  //return par[0] + par[1]*x[0] + par[11]*TMath::Sqrt(x[0]*x[0]*sb2 + 2*x[0]*cov + sa2);
}

double LineFit4Error(double *x, double *par)
{
  double sdx1 = par[6] * par[1];
  double sdx2 = par[7] * par[1];
  double sdx3 = par[8] * par[1];
  double sdx4 = par[9] * par[1];

  double dy1 = par[10];
  double dy2 = par[11];
  double dy3 = par[12];
  double dy4 = par[13];

  double sy1 = TMath::Sqrt(sdx1*sdx1 + dy1*dy1);
  double sy2 = TMath::Sqrt(sdx2*sdx2 + dy2*dy2);
  double sy3 = TMath::Sqrt(sdx3*sdx3 + dy3*dy3);
  double sy4 = TMath::Sqrt(sdx4*sdx4 + dy4*dy4);

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];

  double bracketX = 1.0/4.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4)); // [x]
  double bracketX2 = 1.0/4.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4)); // [x2]
  double bracket1 = 1.0/4.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/4.0 / bracket1;
  //double sa2 = 1.0/3.0 * bracketX2 / (bracketX2*bracket1 - bracketX*bracketX);
  double sb2 = 1.0/4.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  //cout << sa2 << "  " << sb2 << "  " << cov << endl;

  //cout <<x[0]<<"  "<<sb2<<"  "<<sa2<<"  "<<cov<<"  "<< (x[0]+cov)*(x[0]+cov)*sb2 + 2*(x[0]+cov)*cov + sa2 << endl;
  //return par[0] + par[1]*x[0] + par[2]*TMath::Sqrt((x[0]+par[4])*(x[0]+par[4])*par[3] + 2*(x[0]+par[4])*par[4] + par[5]);
  //return par[0] + par[1]*x[0] + par[11]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + 2*(x[0]+cov)*cov + sa2);
  return par[0] + par[1]*x[0] + par[14]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
  //return par[0] + par[1]*x[0] + par[11]*TMath::Sqrt(x[0]*x[0]*sb2 + 2*x[0]*cov + sa2);
}
double ErrFunctionConv(double *x, double *par)
{
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  // Control constants
  Double_t np = 1000.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t fErrf;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of erf and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fErrf = par[0]*(1 + TMath::Erf((xx-par[1])/par[2]));
     sum += fErrf * TMath::Gaus(xx,x[0],par[3]);

     xx = xupp - (i-.5) * step;
     fErrf = par[0]*(1 + TMath::Erf((xx-par[1])/par[2]));
     sum += fErrf * TMath::Gaus(xx,x[0],par[3]);
  }

  return (step * sum * invsq2pi / par[3]);
}

void ScintillationSpectrum()
{
  t->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1782.root");

  // 214Po data (23 days)
  TChain *t2 = new TChain("t");
  t2->Add("../../BiPo/analysis/SecondDataPeriod/bga_all.root");

  // 214Po data (41 days)
  TChain *t3 = new TChain("t");
  t3->Add("../../BiPo/analysis/SecondDataPeriod/Reprocessed_110817/bga_all_1569_2025.root");

  h214Po_epcl->SetLineColor(kRed);
  hEPCL1->SetLineColor(kRed);
  hEPCL2->SetLineColor(kRed);

// ************************************************************************************************
// This part generates the scintillation spectra for ncl = 0, ncl = 1, ncl <= 1 and the fit function
// to the background

  hNCL_EQ1->SetLineColor(kRed);
  //hNCL_EQ1->SetLineStyle(2);

  hNCL_EQ0->SetLineColor(kBlue);
  //hNCL_EQ0->SetLineStyle(2);

  hNCL_EQ1_fid->SetLineColor(kBlue);
  //hNCL_EQ1_fid->SetLineStyle(2);

  t->Draw("fcsc>>hNCL_EQ1","fncl == 1");
  t->Draw("fcsc>>hNCL_LT1");
  t->Draw("fcsc>>hNCL_EQ0","fncl == 0");
  t->Draw("fcscCorr>>hNCL_EQ1_fid","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -172 && fZ < -20) || (fZ > 20 && fZ < 172))");

  // fit function to the background excluding the region of the peak
  TF1 *fit = new TF1("fit",fitFunction,3500,40000,2);
  TF1 *background = new TF1("background","expo",3500,30000);

  hNCL_EQ0->Fit("fit","rn");

  double par[2];
  fit->GetParameters(par);
  background->SetParameters(par);

  background->SetLineColor(kGreen);
  background->SetLineWidth(2);

// ************************************************************************************************
// Subtract background from ncl <= 0 spectrum in order to determine the scaling between the ncl = 1
// and ncl <= 1 spectra. Alos subtract the background from the ncl = 0 spectrum to determine the events
// with no charge cluster

  // create clone of sum spectrum and subtract background
  TH1F *hNCL_LT1_copy = (TH1F*)hNCL_LT1->Clone("hNCL_LT1_copy");
  hNCL_LT1_copy->Add(background,-1);

  // create clone of diff spectrum and divide by ncl == 1 spectrum
  TH1F *hNCL_LT1_div = (TH1F*)hNCL_LT1_copy->Clone("hNCL_LT1_div");
  hNCL_LT1_div->Divide(hNCL_EQ1);

  // create copy of fiducial spectrum and multiply by div spectrum
  TH1F *hNCL_EQ1_fid_copy = (TH1F*)hNCL_EQ1_fid->Clone("hNCL_EQ1_fid_copy");
  hNCL_EQ1_fid_copy->Multiply(hNCL_LT1_div);

  hNCL_EQ1_fid_copy->SetLineColor(kGray);
  hNCL_EQ1_fid_copy->SetFillColor(kGray);
  //hNCL_EQ1_fid_copy->SetLineStyle(2);

  // create clone of ncl == 0 spectrum and subtract background
  TH1F *hNCL_EQ0_copy = (TH1F*)hNCL_EQ0->Clone("hNCL_EQ0_copy");
  hNCL_EQ0_copy->Add(background,-1);

// ************************************************************************************************
// Get alpha-alpha coincidences within 300s

  GetConcidences();

  // zero out negative bins
  for (int i = 0; i < 100; i++) {
     if (i < 38) {hNCL_EQ1_fid_copy->SetBinContent(i+1,0);}
     double BinContent1 = hNCL_EQ1_fid_copy->GetBinContent(i+1);
     if (BinContent1 < 0) {hNCL_EQ1_fid_copy->SetBinContent(i+1,0);}

     double BinContent2 = hCSC_sum->GetBinContent(i+1);
     if (BinContent2 < 0) {hCSC_sum->SetBinContent(i+1,0);}
  }

  // zero out bins below 10000 photon counts
  for (int i = 0; i < 30; i++) {hNCL_EQ1_fid->SetBinContent(i+1,0);}

// ************************************************************************************************
// This part creates the 214Po spectrum form the Bi-Po analysis and corrects the scintillation signal
// for z-position. It also generates the raw ionization spectrum.
  //h214Po_fid->SetLineColor(kGreen);
  //h214Po_fid->SetFillColor(kGreen);
  //h214Po_fid->SetFillStyle(3004);

  t2->Draw("fcsc2>>h214Po");
  t2->Draw("fcsc2>>h214Po","fepcl2 > 0");
  //t2->Draw("fcsc2>>h214Po_fid","TMath::Sqrt(fx2*fx2 + fy2*fy2) < 163 && ((fz2 > -172 && fz2 < -20) || (fz2 > 20 && fz2 < 172))");
  t2->Draw("fercl2>>h214Po_eccl","TMath::Sqrt(fx2*fx2 + fy2*fy2) < 163&& ((fz2 > -172 && fz2 < -20) || (fz2 > 20 && fz2 < 172))");
  t2->Draw("fepcl2>>h214Po_epcl","TMath::Sqrt(fx2*fx2 + fy2*fy2) < 163&& ((fz2 > -172 && fz2 < -20) || (fz2 > 20 && fz2 < 172))");

  double fx2;
  double fy2;
  double fz2;
  double fercl2;
  double fepcl2;
  double fcsc2;
  t2->SetBranchAddress("fx2",&fx2);
  t2->SetBranchAddress("fy2",&fy2);
  t2->SetBranchAddress("fz2",&fz2);
  t2->SetBranchAddress("fercl2",&fercl2);
  t2->SetBranchAddress("fepcl2",&fepcl2);
  t2->SetBranchAddress("fcsc2",&fcsc2);

  // apply z-correction
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
  for (int i = 0; i < t2->GetEntries(); i++) {
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

     h214Po_fid->Fill(csc_corr);
     h214Po_CL->Fill(csc_corr/fercl2);
  }

// ************************************************************************************************
// This part calculates the efficiency curve using the three population of alpha events (222Rn, 218Po, 214Po)

  // fit the charge and scintillation spectra of the three alpha population to a Gaussian
  TF1 *fitAlpha1_eccl = new TF1("fitAlpha1_eccl","gaus",0,200);
  TF1 *fitAlpha2_eccl = new TF1("fitAlpha2_eccl","gaus",0,200);
  TF1 *fit214Po_eccl = new TF1("fitPo214_eccl","gaus",0,300);
  TF1 *fitAlpha1_csc = new TF1("fitAlpha1_csc","gaus",15000,30000);
  TF1 *fitAlpha2_csc = new TF1("fitAlpha2_csc","gaus",15000,30000);
  TF1 *fit214Po_csc = new TF1("fitPo214_csc","gaus",25000,40000);

  TF1 *fitAlpha1_cl = new TF1("fitAlpha1_cl","gaus",100,400);
  TF1 *fitAlpha2_cl = new TF1("fitAlpha2_cl","gaus",100,400);
  TF1 *fit214Po_cl = new TF1("fit214Po_cl","gaus",100,400);

  TF1 *effCharge = new TF1("effCharge",ErrFunction,0,180,3);
  effCharge->SetParameters(0.5,57.5,15.81);
  effCharge->SetLineWidth(1);
  effCharge->SetLineColor(kBlack);

  hECCL1->Divide(effCharge);
  hECCL2->Divide(effCharge);
  h214Po_eccl->Divide(effCharge);

  hECCL1->Fit("fitAlpha1_eccl","rL");
  hECCL2->Fit("fitAlpha2_eccl","rL");
  h214Po_eccl->Fit("fitPo214_eccl","rL");
  hCSC1->Fit("fitAlpha1_csc","rL");
  hCSC2->Fit("fitAlpha2_csc","rL");
  h214Po_fid->Fit("fitPo214_csc","rL");

  hCL1->Fit("fitAlpha1_cl","rnL");
  hCL2->Fit("fitAlpha2_cl","rnL");
  h214Po_CL->Fit("fit214Po_cl","rnL");

  // create graph of the mean scintillation vs mean ionization of the three populations
  TGraphErrors *gr = new TGraphErrors(4);
  gr->SetMarkerStyle(7);
  gr->SetPoint(0,fitAlpha1_eccl->GetParameter(1),fitAlpha1_csc->GetParameter(1));
  gr->SetPoint(1,fitAlpha2_eccl->GetParameter(1),fitAlpha2_csc->GetParameter(1));
  gr->SetPoint(2,fit214Po_eccl->GetParameter(1),fit214Po_csc->GetParameter(1));
  gr->SetPoint(3,0,0);

  gr->SetPointError(0,fitAlpha1_eccl->GetParError(1),fitAlpha1_csc->GetParError(1));
  gr->SetPointError(1,fitAlpha2_eccl->GetParError(1),fitAlpha2_csc->GetParError(1));
  gr->SetPointError(2,fit214Po_eccl->GetParError(1),fit214Po_csc->GetParError(1));

  // create graph of mean ionization vs true energy of the alphas
  TGraphErrors *grCharge = new TGraphErrors(3);
  grCharge->SetMarkerStyle(7);
  grCharge->SetPoint(0,5590,fitAlpha1_eccl->GetParameter(1));
  grCharge->SetPoint(1,6115,fitAlpha2_eccl->GetParameter(1));
  grCharge->SetPoint(2,7833,fit214Po_eccl->GetParameter(1));

  grCharge->SetPointError(0,0,fitAlpha1_eccl->GetParError(1));
  grCharge->SetPointError(1,0,fitAlpha2_eccl->GetParError(1));
  grCharge->SetPointError(2,0,fit214Po_eccl->GetParError(1));

  // create graph of mean scintillation vs true energy of the alphas
  TGraphErrors *grScint = new TGraphErrors(3);
  grScint->SetMarkerStyle(7);
  grScint->SetPoint(0,5590,fitAlpha1_csc->GetParameter(1));
  grScint->SetPoint(1,6115,fitAlpha2_csc->GetParameter(1));
  grScint->SetPoint(2,7833,fit214Po_csc->GetParameter(1));

  grScint->SetPointError(0,0,fitAlpha1_csc->GetParError(1));
  grScint->SetPointError(1,0,fitAlpha2_csc->GetParError(1));
  grScint->SetPointError(2,0,fit214Po_csc->GetParError(1));

  // create graph of mean ionization vs true energy of the alphas
  TGraphErrors *grIonScint = new TGraphErrors(4);
  grIonScint->SetMarkerStyle(7);
  grIonScint->SetPoint(0,353.9,0); // 4.27 MeV
  grIonScint->SetPoint(1,307.1,fitAlpha1_cl->GetParameter(1)); // 5.59 MeV
  grIonScint->SetPoint(2,291.6,fitAlpha2_cl->GetParameter(1)); // 6.115 MeV
  grIonScint->SetPoint(3,252.2,fit214Po_cl->GetParameter(1)); // 7.833 MeV

  grIonScint->SetPointError(0,0,0);
  grIonScint->SetPointError(1,0,fitAlpha1_cl->GetParError(1));
  grIonScint->SetPointError(2,0,fitAlpha2_cl->GetParError(1));
  grIonScint->SetPointError(3,0,fit214Po_cl->GetParError(1));

  // fit functions to the three graphs
  TF1 *fitLinChargeScint = new TF1("fitLinChargeScint","[0]+[1]*x",2,500);
  TF1 *fitLinCharge = new TF1("fitLinCharge","[0]+[1]*x");
  TF1 *fitLinScint = new TF1("fitLinScint","[0]+[1]*x");
  TF1 *fitLinIonScint = new TF1("fitLinIonScint","[0]+[1]*x");

  fitLinChargeScint->SetLineWidth(1);
  fitLinChargeScint->SetLineStyle(2);
  fitLinCharge->SetLineWidth(1);
  fitLinCharge->SetLineStyle(2);
  fitLinScint->SetLineWidth(1);
  fitLinScint->SetLineStyle(2);
  fitLinIonScint->SetLineWidth(1);
  fitLinIonScint->SetLineStyle(2);

  gr->Fit("fitLinChargeScint","r");
  grCharge->Fit("fitLinCharge");
  grScint->Fit("fitLinScint");
  grIonScint->Fit("fitLinIonScint");

  TMarker *U238DataPoint = new TMarker(353.9,fitLinIonScint->Eval(353.9),20);
  U238DataPoint->SetMarkerSize(0.8);
  U238DataPoint->SetMarkerColor(kRed);

  // create same graph with zero x error but instead convert x error into y error
  TGraphErrors *gr2 = new TGraphErrors(4);
  gr2->SetMarkerStyle(7);
  gr2->SetMarkerColor(kRed);
  gr2->SetLineColor(kRed);
  gr2->SetPoint(0,fitAlpha1_eccl->GetParameter(1),fitAlpha1_csc->GetParameter(1));
  gr2->SetPoint(1,fitAlpha2_eccl->GetParameter(1),fitAlpha2_csc->GetParameter(1));
  gr2->SetPoint(2,fit214Po_eccl->GetParameter(1),fit214Po_csc->GetParameter(1));
  gr2->SetPoint(3,0,0);

  // create error curve on the efficiency fit
  TF1 *fitLinError1 = new TF1("fitLinError1",LineFitError,0,180,12);
  fitLinError1->SetParameter(0,fitLinChargeScint->GetParameter(0));
  fitLinError1->SetParameter(1,fitLinChargeScint->GetParameter(1));
  fitLinError1->SetParameter(2,fitAlpha1_eccl->GetParameter(1));
  fitLinError1->SetParameter(3,fitAlpha2_eccl->GetParameter(1));
  fitLinError1->SetParameter(4,fit214Po_eccl->GetParameter(1));
  fitLinError1->SetParameter(5,fitAlpha1_eccl->GetParError(1));
  fitLinError1->SetParameter(6,fitAlpha2_eccl->GetParError(1));
  fitLinError1->SetParameter(7,fit214Po_eccl->GetParError(1));
  fitLinError1->SetParameter(8,fitAlpha1_csc->GetParError(1));
  fitLinError1->SetParameter(9,fitAlpha2_csc->GetParError(1));
  fitLinError1->SetParameter(10,fit214Po_csc->GetParError(1));
  fitLinError1->SetParameter(11,1);

  fitLinError1->SetLineWidth(1);

  TF1 *fitLinError2 = new TF1("fitLinError2",LineFitError,0,180,12);
  fitLinError2->SetParameter(0,fitLinChargeScint->GetParameter(0));
  fitLinError2->SetParameter(1,fitLinChargeScint->GetParameter(1));
  fitLinError2->SetParameter(2,fitAlpha1_eccl->GetParameter(1));
  fitLinError2->SetParameter(3,fitAlpha2_eccl->GetParameter(1));
  fitLinError2->SetParameter(4,fit214Po_eccl->GetParameter(1));
  fitLinError2->SetParameter(5,fitAlpha1_eccl->GetParError(1));
  fitLinError2->SetParameter(6,fitAlpha2_eccl->GetParError(1));
  fitLinError2->SetParameter(7,fit214Po_eccl->GetParError(1));
  fitLinError2->SetParameter(8,fitAlpha1_csc->GetParError(1));
  fitLinError2->SetParameter(9,fitAlpha2_csc->GetParError(1));
  fitLinError2->SetParameter(10,fit214Po_csc->GetParError(1));
  fitLinError2->SetParameter(11,-1);

  fitLinError2->SetLineWidth(1);

  // create error curve on the efficiency fit
  TF1 *fitIonScintError1 = new TF1("fitIonScintError1",LineFitError,240,360,12);
  fitIonScintError1->SetParameter(0,fitLinIonScint->GetParameter(0));
  fitIonScintError1->SetParameter(1,fitLinIonScint->GetParameter(1));
  fitIonScintError1->SetParameter(2,252.2);
  fitIonScintError1->SetParameter(3,291.6);
  fitIonScintError1->SetParameter(4,307.1);
  fitIonScintError1->SetParameter(5,0);
  fitIonScintError1->SetParameter(6,0);
  fitIonScintError1->SetParameter(7,0);
  fitIonScintError1->SetParameter(8,fit214Po_cl->GetParError(1));
  fitIonScintError1->SetParameter(9,fitAlpha2_cl->GetParError(1));
  fitIonScintError1->SetParameter(10,fitAlpha1_cl->GetParError(1));
  fitIonScintError1->SetParameter(11,1);

  fitIonScintError1->SetLineWidth(1);

  TF1 *fitIonScintError2 = new TF1("fitIonScintError2",LineFitError,240,360,12);
  fitIonScintError2->SetParameter(0,fitLinIonScint->GetParameter(0));
  fitIonScintError2->SetParameter(1,fitLinIonScint->GetParameter(1));
  fitIonScintError2->SetParameter(2,252.2);
  fitIonScintError2->SetParameter(3,291.6);
  fitIonScintError2->SetParameter(4,307.1);
  fitIonScintError2->SetParameter(5,0);
  fitIonScintError2->SetParameter(6,0);
  fitIonScintError2->SetParameter(7,0);
  fitIonScintError2->SetParameter(8,fit214Po_cl->GetParError(1));
  fitIonScintError2->SetParameter(9,fitAlpha2_cl->GetParError(1));
  fitIonScintError2->SetParameter(10,fitAlpha1_cl->GetParError(1));
  fitIonScintError2->SetParameter(11,-1);

  fitIonScintError2->SetLineWidth(1);

  TGraph *grTMP1 = new TGraph(4);
  grTMP1->SetMarkerStyle(7);
  grTMP1->SetPoint(0,79.3-8.1,17719.4);
  grTMP1->SetPoint(1,fitAlpha1_eccl->GetParameter(1),fitAlpha1_csc->GetParameter(1));
  grTMP1->SetPoint(2,fitAlpha2_eccl->GetParameter(1),fitAlpha2_csc->GetParameter(1));
  grTMP1->SetPoint(3,fit214Po_eccl->GetParameter(1),fit214Po_csc->GetParameter(1));

  TF1 *fitTMP1 = new TF1("fitTMP1","[0]+[1]*x+[2]*x*x",0,180);
  fitTMP1->SetLineColor(kBlack);
  fitTMP1->SetLineWidth(1);
  fitTMP1->SetLineStyle(3);

  grTMP1->Fit("fitTMP1","rn");

  TGraph *grTMP2 = new TGraph(4);
  grTMP2->SetMarkerStyle(7);
  grTMP2->SetPoint(0,79.3+8.1,16271.5);
  grTMP2->SetPoint(1,fitAlpha1_eccl->GetParameter(1),fitAlpha1_csc->GetParameter(1));
  grTMP2->SetPoint(2,fitAlpha2_eccl->GetParameter(1),fitAlpha2_csc->GetParameter(1));
  grTMP2->SetPoint(3,fit214Po_eccl->GetParameter(1),fit214Po_csc->GetParameter(1));

  TF1 *fitTMP2 = new TF1("fitTMP2","[0]+[1]*x+[2]*x*x",0,180);
  fitTMP2->SetLineColor(kBlack);
  fitTMP2->SetLineWidth(1);
  fitTMP2->SetLineStyle(3);

  grTMP2->Fit("fitTMP2","rn"),

  fitTMP1->SetRange(0,105);
  fitTMP2->SetRange(0,105);

/*
  grTMP1->SetPointError(0,0.0,0.0);
  grTMP1->SetPointError(1,fitAlpha1_eccl->GetParError(1),fitAlpha1_csc->GetParError(1));
  grTMP1->SetPointError(2,fitAlpha2_eccl->GetParError(1),fitAlpha2_csc->GetParError(1));
  grTMP1->SetPointError(3,fit214Po_eccl->GetParError(1),fit214Po_csc->GetParError(1));
*/
  TF1 *fitLinErrorTMP1 = new TF1("fitLinErrorTMP1",LineFit4Error,0,180,15);
  fitLinErrorTMP1->SetParameter(0,fitLinChargeScint->GetParameter(0));
  fitLinErrorTMP1->SetParameter(1,fitLinChargeScint->GetParameter(1));
  fitLinErrorTMP1->SetParameter(2,79.3);
  fitLinErrorTMP1->SetParameter(3,fitAlpha1_eccl->GetParameter(1));
  fitLinErrorTMP1->SetParameter(4,fitAlpha2_eccl->GetParameter(1));
  fitLinErrorTMP1->SetParameter(5,fit214Po_eccl->GetParameter(1));
  fitLinErrorTMP1->SetParameter(6,8.1);
  fitLinErrorTMP1->SetParameter(7,fitAlpha1_eccl->GetParError(1));
  fitLinErrorTMP1->SetParameter(8,fitAlpha2_eccl->GetParError(1));
  fitLinErrorTMP1->SetParameter(9,fit214Po_eccl->GetParError(1));
  fitLinErrorTMP1->SetParameter(10,724.0);
  fitLinErrorTMP1->SetParameter(11,fitAlpha1_csc->GetParError(1));
  fitLinErrorTMP1->SetParameter(12,fitAlpha2_csc->GetParError(1));
  fitLinErrorTMP1->SetParameter(13,fit214Po_csc->GetParError(1));
  fitLinErrorTMP1->SetParameter(14,1);

  fitLinErrorTMP1->SetLineWidth(1);
  fitLinErrorTMP1->SetLineColor(kRed);

  TF1 *fitLinErrorTMP2 = new TF1("fitLinErrorTMP2",LineFit4Error,0,180,15);
  fitLinErrorTMP2->SetParameter(0,fitLinChargeScint->GetParameter(0));
  fitLinErrorTMP2->SetParameter(1,fitLinChargeScint->GetParameter(1));
  fitLinErrorTMP2->SetParameter(2,79.3);
  fitLinErrorTMP2->SetParameter(3,fitAlpha1_eccl->GetParameter(1));
  fitLinErrorTMP2->SetParameter(4,fitAlpha2_eccl->GetParameter(1));
  fitLinErrorTMP2->SetParameter(5,fit214Po_eccl->GetParameter(1));
  fitLinErrorTMP2->SetParameter(6,8.1);
  fitLinErrorTMP2->SetParameter(7,fitAlpha1_eccl->GetParError(1));
  fitLinErrorTMP2->SetParameter(8,fitAlpha2_eccl->GetParError(1));
  fitLinErrorTMP2->SetParameter(9,fit214Po_eccl->GetParError(1));
  fitLinErrorTMP2->SetParameter(10,724.0);
  fitLinErrorTMP2->SetParameter(11,fitAlpha1_csc->GetParError(1));
  fitLinErrorTMP2->SetParameter(12,fitAlpha2_csc->GetParError(1));
  fitLinErrorTMP2->SetParameter(13,fit214Po_csc->GetParError(1));
  fitLinErrorTMP2->SetParameter(14,-1);

  fitLinErrorTMP2->SetLineWidth(1);
  fitLinErrorTMP2->SetLineColor(kRed);

  TGraphErrors *grTMP = new TGraphErrors(1);
  grTMP->SetPoint(0,79.3,16863.3);
  grTMP->SetPointError(0,8.1,14);
  grTMP->SetMarkerStyle(20);
  grTMP->SetMarkerColor(kRed);
  grTMP->SetLineColor(kRed);

  double sig = 22.5;

  // create efficiency curves
  TF1 *effScint = new TF1("effScint",ErrFunction,0,40000,3);
  effScint->SetParameters(0.5,fitLinChargeScint->Eval(57.5),fitLinChargeScint->Eval(57.5+sig/2.0) - fitLinChargeScint->Eval(57.5-sig/2.0));
  effScint->SetLineStyle(2);
  effScint->SetLineWidth(1);

  TF1 *effScintErr1 = new TF1("effScintErr1",ErrFunction,0,40000,3);
  effScintErr1->SetParameters(0.5,fitLinErrorTMP1->Eval(57.5),fitLinErrorTMP1->Eval(57.5+sig/2.0) - fitLinErrorTMP1->Eval(57.5-sig/2.0));
  effScintErr1->SetLineStyle(1);
  effScintErr1->SetLineWidth(1);

  TF1 *effScintErr2 = new TF1("effScintErr2",ErrFunction,0,40000,3);
  effScintErr2->SetParameters(0.5,fitLinErrorTMP2->Eval(57.5),fitLinErrorTMP2->Eval(57.5+sig/2.0) - fitLinErrorTMP2->Eval(57.5-sig/2.0));
  effScintErr2->SetLineStyle(1);
  effScintErr2->SetLineWidth(1);

/*  TF1 *effCharge = new TF1("effCharge",ErrFunction,0,180,3);
  effCharge->SetParameters(0.5,57.5,15.81);
  effCharge->SetLineWidth(1);
  effCharge->SetLineColor(kBlue);
*/
  TF1 *effChargeConv = new TF1("effChargeConv",ErrFunctionConv,0,180,4);
  effChargeConv->SetParameters(0.5,57.5,15.81,fitAlpha1_eccl->GetParameter(2));
  effChargeConv->SetLineWidth(1);
  effChargeConv->SetLineColor(kBlack);
  effChargeConv->SetLineStyle(2);

  // create copy of scintillation spectrum and scale according to efficiency
  TH1F *hNCL_EQ1_fid_eff = (TH1F*)hNCL_EQ1_fid->Clone("hNCL_EQ1_fid_eff");
  hNCL_EQ1_fid_eff->Divide(effScint);
  hNCL_EQ1_fid_eff->SetLineColor(kRed);

// ************************************************************************************************
// Print out some results
  double IntegralError1;
  double Integral1 = hNCL_EQ1_fid_copy->IntegralAndError(1,100,IntegralError1);

  cout << "Total number of alphas: " << Integral1 << " +- " << IntegralError1 << endl;

  cout << "alpha\tion\tscint\tratio" << endl;
  cout << "Rn222:\t" << fitAlpha1_eccl->GetParameter(1) << "\t" << fitAlpha1_csc->GetParameter(1) << "\t" << fitAlpha1_eccl->GetParameter(1)/fitAlpha1_csc->GetParameter(1) << endl;
  cout << "Po218:\t" << fitAlpha2_eccl->GetParameter(1) << "\t" << fitAlpha2_csc->GetParameter(1) << "\t" << fitAlpha2_eccl->GetParameter(1)/fitAlpha2_csc->GetParameter(1) << endl;
  cout << "Po214:\t" << fit214Po_eccl->GetParameter(1) << "\t" << fit214Po_csc->GetParameter(1) << "\t" << fit214Po_eccl->GetParameter(1)/fit214Po_csc->GetParameter(1) << endl;

  cout << "efficiency paramters: mean = " << fitLinChargeScint->Eval(57.5) << "  sigma = " << fitLinChargeScint->Eval(57.5+sig/2.0) - fitLinChargeScint->Eval(57.5-sig/2.0) << endl;

  cout << "4.27 MeV := " << fitLinScint->Eval(4270) << "  " << fitLinCharge->Eval(4270) << endl;
  cout << "5.59 MeV := " << fitLinScint->Eval(5590) << "  " << fitLinCharge->Eval(5590) << endl;
  cout << "6.115 MeV := " << fitLinScint->Eval(6115) << "  " << fitLinCharge->Eval(6115) << endl;
  cout << "7.833 MeV := " << fitLinScint->Eval(7833) << "  " << fitLinCharge->Eval(7833) << endl;

  // old value 16250
  cout << "Efficiency at 4.27 MeV: " << effScint->Eval(17334) << " + " << effScintErr1->Eval(17334) - effScint->Eval(17334) << " - " << effScint->Eval(17334) - effScintErr2->Eval(17334) << endl;
  cout << "Efficiency at 5.5 MeV: " << effScint->Eval(22473) << " + " << effScintErr1->Eval(22473) - effScint->Eval(22473) << " - " << effScint->Eval(22473) - effScintErr2->Eval(22473) << endl;

  cout << fitLinChargeScint->Eval(57.5) << " +- " << fitLinChargeScint->Eval(57.5+15.81/2.0) - fitLinChargeScint->Eval(57.5-15.81/2.0) << endl;

  cout << "S/I(1) = " << fitIonScintError1->Eval(353.9) << "  S/I(2) = " << fitIonScintError2->Eval(353.9) << endl;
  cout << "change in quenching: " << fitIonScintError1->Eval(353.9)/fitLinIonScint->Eval(353.9) << "  " << fitIonScintError2->Eval(353.9)/fitLinIonScint->Eval(353.9) << endl;

  cout << fitLinChargeScint->Eval(79.3) << "  " << fitLinError1->Eval(79.3) << "  " << fitLinError2->Eval(79.3) << endl;

  cout << "**** Efficiency results *********************************************" << endl;
  cout << "4.27 MeV: " << effChargeConv->Eval(79.3) << " + " << effChargeConv->Eval(79.3+8.1) - effChargeConv->Eval(79.3) << " - " <<  effChargeConv->Eval(79.3) - effChargeConv->Eval(79.3-8.1) << endl;
  cout << "5.59 MeV: " << effChargeConv->Eval(fitAlpha1_eccl->GetParameter(1)) << " + " << effChargeConv->Eval(fitAlpha1_eccl->GetParameter(1)+fitAlpha1_eccl->GetParError(1)) - effChargeConv->Eval(fitAlpha1_eccl->GetParameter(1)) << " - " << effChargeConv->Eval(fitAlpha1_eccl->GetParameter(1)) - effChargeConv->Eval(fitAlpha1_eccl->GetParameter(1)-fitAlpha1_eccl->GetParError(1)) << endl;
  cout << "6.115 MeV: " << effChargeConv->Eval(fitAlpha2_eccl->GetParameter(1)) << " + " << effChargeConv->Eval(fitAlpha2_eccl->GetParameter(1)+fitAlpha2_eccl->GetParError(1)) - effChargeConv->Eval(fitAlpha2_eccl->GetParameter(1)) << " - " << effChargeConv->Eval(fitAlpha2_eccl->GetParameter(1)) - effChargeConv->Eval(fitAlpha2_eccl->GetParameter(1)-fitAlpha2_eccl->GetParError(1)) << endl;
  cout << "7.833 MeV: " << effChargeConv->Eval(fit214Po_eccl->GetParameter(1)) << " + " << effChargeConv->Eval(fit214Po_eccl->GetParameter(1)+fit214Po_eccl->GetParError(1)) - effChargeConv->Eval(fit214Po_eccl->GetParameter(1)) << " - " << effChargeConv->Eval(fit214Po_eccl->GetParameter(1)) - effChargeConv->Eval(fit214Po_eccl->GetParameter(1)-fit214Po_eccl->GetParError(1)) << endl;

  effScint->SetParameter(0,10);
  effScintErr1->SetParameter(0,10);
  effScintErr2->SetParameter(0,10);

// ************************************************************************************************
// Create legends for the various plots

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

  TLegend *l3 = new TLegend(0,0,0.2,0.2);
  l3->AddEntry(hNCL_EQ1_fid,"ncl = 1 (fiducial)");
  l3->AddEntry(hCSC1,"Rn222 (coincidence)");
  l3->AddEntry(hCSC2,"Po218 (coincidence)");
  l3->AddEntry(h214Po_fid,"Po214 (Bi-Po)");
  l3->AddEntry(hCSC_sum,"sum");

  TLegend *l4 = new TLegend(0,0,0.2,0.2);
  l4->AddEntry(hNCL_EQ1_fid_copy,"ncl <= 1 (fiducial)");
  l4->AddEntry(hNCL_EQ1_fid,"ncl == 1 (fiducial)");
  l4->AddEntry(h214Po_fid,"Po214 (fiducial)");

// ************************************************************************************************
// Generate plots

  TCanvas *c1 = new TCanvas("c1","ncl = 1, ncl = 0, ncl <= 0, background fit");
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
  hNCL_EQ1_fid->Draw();
  hNCL_EQ1_fid_copy->Draw("same");
  h214Po_fid->Draw("same");
  l4->Draw("same");

  TCanvas *c6 = new TCanvas("c6","Scintillation spectrum");
  hCSC_sum->Draw();
  hCSC1->Draw("same");
  hCSC2->Draw("same");
  h214Po_fid->Draw("same");
  hNCL_EQ1_fid->Draw("same");
  effScint->Draw("same");
  effScintErr1->Draw("same");
  effScintErr2->Draw("same");
  l3->Draw("same");

  TCanvas *c7 = new TCanvas("c7","Alpha coincidence dt");
  hAlphaDT->Draw();

  TCanvas *c8 = new TCanvas("c8","XY distribution");
  hXY->Draw();

  TCanvas *c9 = new TCanvas("c9","Scintillation spectrum sub");
  hCSC_sub->Draw();

  TCanvas *c10 = new TCanvas("c10","Rn222, Po218, Po214 charge and scintillation");
  c10->Divide(3,2);
  c10->cd(1);
  hCSC1->Draw();
  fitAlpha1_csc->Draw("same");
  c10->cd(2);
  hCSC2->Draw();
  fitAlpha2_csc->Draw("same");
  c10->cd(3);
  h214Po_fid->Draw();
  fit214Po_csc->Draw("same");
  c10->cd(4);
  hECCL1->Draw();
  fitAlpha1_eccl->Draw("same");
  hEPCL1->Draw("same");
  c10->cd(5);
  hECCL2->Draw();
  fitAlpha2_eccl->Draw("same");
  hEPCL2->Draw("same");
  c10->cd(6);
  h214Po_eccl->Draw();
  fit214Po_eccl->Draw("same");
  h214Po_epcl->Draw("same");

  TCanvas *c11 = new TCanvas("c11","Charge collection efficiency");
  //hNCL_EQ1_fid_eff->Draw();
  hNCL_EQ1_fid->Draw();
  effScint->Draw("same");
  effScintErr1->Draw("same");
  effScintErr2->Draw("same");

  TCanvas *c12 = new TCanvas("c12","Mean scintillation vs mean ionization");
  c12->Divide(1,2);
  c12->cd(1);
  fitLinError1->Draw();
  fitLinError2->Draw("same");
  gr->Draw("Psame");
  grTMP->Draw("Psame");
  //fitLinErrorTMP1->Draw("same");
  //fitLinErrorTMP2->Draw("same");
  fitTMP1->Draw("same");
  fitTMP2->Draw("same");
  //grTMP1->Draw("Psame");
  //grTMP2->Draw("Psame");
  //gr2->Draw("Psame");
  //fitLinError1->Draw("same");
  //fitLinError2->Draw("same");
  c12->cd(2);
  effCharge->Draw();
  effChargeConv->Draw("same");

  TCanvas *c13 = new TCanvas("c13","Mean ionization vs true energy");
  grCharge->Draw("AP");

  TCanvas *c14 = new TCanvas("c14","Mean scintillation vs true energy");
  grScint->Draw("AP");

  TCanvas *c15 = new TCanvas("c15","Ion/Scint");
  hCL1->Draw();
  hCL2->Draw("same");
  h214Po_CL->Draw("same");
  fitAlpha1_cl->Draw("same");
  fitAlpha2_cl->Draw("same");
  fit214Po_cl->Draw("same");

  TCanvas *c16 = new TCanvas("c16","ion/scint vs true energy");
  grIonScint->Draw("AP");
  fitIonScintError1->Draw("same");
  fitIonScintError2->Draw("same");
  U238DataPoint->Draw("same");

  TCanvas *c17 = new TCanvas();
  grTMP1->Draw("AP");
  grTMP2->Draw("Psame");

  return;
}

void GetConcidences()
{
  //hCSC1->SetLineColor(kRed);
  //hCSC1->SetFillColor(kRed);
  //hCSC1->SetFillStyle(3004);

  //hCSC2->SetLineColor(kBlue);
  //hCSC2->SetFillColor(kBlue);
  //hCSC2->SetFillStyle(3005);

  int fncl;
  double fX;
  double fY;
  double fZ;
  double fcsc;
  double fcscCorr;
  double ftsc;
  double feccl;
  double fercl;
  int ftrigsec;
  int ftrigsub;
  int ftrigoff;

  t->SetBranchAddress("fncl",&fncl);
  t->SetBranchAddress("fX",&fX);
  t->SetBranchAddress("fY",&fY);
  t->SetBranchAddress("fZ",&fZ);
  t->SetBranchAddress("fcsc",&fcsc);
  t->SetBranchAddress("fcscCorr",&fcscCorr);
  t->SetBranchAddress("ftsc",&ftsc);
  t->SetBranchAddress("feccl",&feccl);
  t->SetBranchAddress("fercl",&fercl);
  t->SetBranchAddress("ftrigsec",&ftrigsec);
  t->SetBranchAddress("ftrigsub",&ftrigsub);
  t->SetBranchAddress("ftrigoff",&ftrigoff);

  // arrays holding the alpha charge cluster object and time information
  double AlphaX[100000];
  double AlphaY[100000];
  double AlphaZ[100000];
  double AlphaTsc[100000];
  double AlphaCSC[100000];
  double AlphaERCL[100000];
  double AlphaECCL[100000];
  int AlphaTrigSec[100000];
  int AlphaTrigSub[100000];
  int AlphaTrigOff[100000];

  // initialize arrays
  for (int i = 0; i < 10000; i++) {
     AlphaX[i] = 0.0;
     AlphaY[i] = 0.0;
     AlphaZ[i] = 0.0;
     AlphaTsc[i] = 0.0;
     AlphaCSC[i] = 0.0;
     AlphaERCL[i] = 0.0;
     AlphaECCL[i] = 0.0;
     AlphaTrigSec[i] = 0;
     AlphaTrigSub[i] = 0;
     AlphaTrigOff[i] = 0;
  }

  int nentries = 0;
  // loop over all alpha events
  for (int i = 0; i < t->GetEntries(); i++) {
     t->GetEntry(i);

     if (fncl != 1) {continue;}

     AlphaX[nentries] = fX;
     AlphaY[nentries] = fY;
     AlphaZ[nentries] = fZ;
     AlphaTsc[nentries] = ftsc;
     AlphaCSC[nentries] = fcscCorr;
     AlphaERCL[nentries] = fercl;
     AlphaECCL[nentries] = feccl;
     AlphaTrigSec[nentries] = ftrigsec;
     AlphaTrigSub[nentries] = ftrigsub;
     AlphaTrigOff[nentries] = ftrigoff;

     nentries++;
  }

  // look for coincidences
  double RCut = 28; // position cut (28 mm)
  double TCut = 300; // time cut (5 min)

  int nCoincidences = 0;
  for (int i = 0; i < nentries-1; i++) {
     for (int k = i+1; k < nentries; k++) {
        if (TMath::Abs(AlphaTrigSec[i] - AlphaTrigSec[k]) > TCut) {continue;}
        double dtAlpha = 0.0;
        double AlphaCSC1 = 0.0;
        double AlphaCSC2 = 0.0;
        double AlphaERCL1 = 0.0;
        double AlphaERCL2 = 0.0;
        double AlphaECCL1 = 0.0;
        double AlphaECCL2 = 0.0;
        double AlphaX1 = 0.0;
        double AlphaX2 = 0.0;
        double AlphaY1 = 0.0;
        double AlphaY2 = 0.0;
        double AlphaZ1 = 0.0;
        double AlphaZ2 = 0.0;
        double tMicro1 = AlphaTsc[i] / 1000 - double(AlphaTrigOff[i]) + double(AlphaTrigSub[i]);
        double tMicro2 = AlphaTsc[k] / 1000 - double(AlphaTrigOff[k]) + double(AlphaTrigSub[k]);

        if (AlphaTrigSec[k] > AlphaTrigSec[i]) {
           dtAlpha = (double(AlphaTrigSec[k]) - double(AlphaTrigSec[i]))*1000000; + (tMicro2 - tMicro1);
           AlphaX1 = AlphaX[i];
           AlphaX2 = AlphaX[k];
           AlphaY1 = AlphaY[i];
           AlphaY2 = AlphaY[k];
           AlphaZ1 = AlphaZ[i];
           AlphaZ2 = AlphaZ[k];
           AlphaCSC1 = AlphaCSC[i];
           AlphaCSC2 = AlphaCSC[k];
           AlphaERCL1 = AlphaERCL[i];
           AlphaERCL2 = AlphaERCL[k];
           AlphaECCL1 = AlphaECCL[i];
           AlphaECCL2 = AlphaECCL[k];
        }
        if (AlphaTrigSec[k] == AlphaTrigSec[i]) {
           //cout << "dt = 0" << endl;
           //continue;
           if (tMicro2 > tMicro1) {
              dtAlpha = tMicro2 - tMicro1;
              AlphaX1 = AlphaX[i];
              AlphaX2 = AlphaX[k];
              AlphaY1 = AlphaY[i];
              AlphaY2 = AlphaY[k];
              AlphaZ1 = AlphaZ[i];
              AlphaZ2 = AlphaZ[k];
              AlphaCSC1 = AlphaCSC[i];
              AlphaCSC2 = AlphaCSC[k];
              AlphaERCL1 = AlphaERCL[i];
              AlphaERCL2 = AlphaERCL[k];
              AlphaECCL1 = AlphaECCL[i];
              AlphaECCL2 = AlphaECCL[k];
           }
           else {
              dtAlpha = tMicro1 - tMicro2;
              AlphaX1 = AlphaX[k];
              AlphaX2 = AlphaX[i];
              AlphaY1 = AlphaY[k];
              AlphaY2 = AlphaY[i];
              AlphaZ1 = AlphaZ[k];
              AlphaZ2 = AlphaZ[i];
              AlphaCSC1 = AlphaCSC[k];
              AlphaCSC2 = AlphaCSC[i];
              AlphaERCL1 = AlphaERCL[k];
              AlphaERCL2 = AlphaERCL[i];
              AlphaECCL1 = AlphaECCL[k];
              AlphaECCL2 = AlphaECCL[i];
           }
        }
        if (AlphaTrigSec[k] < AlphaTrigSec[i]) {
           dtAlpha = (double(AlphaTrigSec[i]) - double(AlphaTrigSec[k]))*1000000 + (tMicro1 - tMicro2);
           AlphaX1 = AlphaX[k];
           AlphaX2 = AlphaX[i];
           AlphaY1 = AlphaY[k];
           AlphaY2 = AlphaY[i];
           AlphaZ1 = AlphaZ[k];
           AlphaZ2 = AlphaZ[i];
           AlphaCSC1 = AlphaCSC[k];
           AlphaCSC2 = AlphaCSC[i];
           AlphaERCL1 = AlphaERCL[k];
           AlphaERCL2 = AlphaERCL[i];
           AlphaECCL1 = AlphaECCL[k];
           AlphaECCL2 = AlphaECCL[i];
        }

        double dR = TMath::Sqrt((AlphaX1 - AlphaX2)*(AlphaX1 - AlphaX2) + (AlphaY1 - AlphaY2)*(AlphaY1 - AlphaY2));

        if (dR > RCut) {continue;}
        //if (AlphaEnergy1 < 4000 || AlphaEnergy2 < 4000) {continue;}

        std::cout << "(" << nCoincidences << ") csc1 = " << AlphaCSC1 << "\tcsc2 = " << AlphaCSC2 << "\tdR = " << dR << "\tdt = " << dtAlpha << std::endl;
        hAlphaDT->Fill(dtAlpha);
        hXY->Fill(AlphaX2,AlphaY2);

        // apply fiducial cut
        if (TMath::Sqrt(AlphaX1*AlphaX1 + AlphaY1*AlphaY1) < 163 && ((AlphaZ1 > -172 && AlphaZ1 < -20) || (AlphaZ1 > 20 && AlphaZ1 < 172))) {hCSC1->Fill(AlphaCSC1); hECCL1->Fill(AlphaERCL1); hEPCL1->Fill(AlphaECCL1); hCL1->Fill(AlphaCSC1/AlphaERCL1);}
        if (TMath::Sqrt(AlphaX2*AlphaX2 + AlphaY2*AlphaY2) < 163 && ((AlphaZ2 > -172 && AlphaZ2 < -20) || (AlphaZ2 > 20 && AlphaZ2 < 172))) {hCSC2->Fill(AlphaCSC2); hECCL2->Fill(AlphaERCL2); hEPCL2->Fill(AlphaECCL2); hCL2->Fill(AlphaCSC2/AlphaERCL2);}

        nCoincidences++;
     }
  }

  TH1F *hCSC_sum = (TH1F*)hCSC1->Clone("hCSC_sum");
  hCSC_sum->Add(hCSC2);
  hCSC_sum->Add(h214Po_fid);
  hCSC_sum->SetLineColor(kGray);
  hCSC_sum->SetFillColor(kGray);
  hCSC_sum->SetFillStyle(1001);

  TH1F *hCSC_sub = (TH1F*)hNCL_EQ1_fid_copy->Clone("hCSC_sub");
  hCSC_sub->Add(hCSC_sum,-1);
  //hCSC_sub->Add(h214Po_fid,-1);

  for (int i = 0; i < 38; i++) {hCSC_sub->SetBinContent(i+1,0);}

  return;
}
