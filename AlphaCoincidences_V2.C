void AlphaCoincidences_V2()
{
  gROOT->SetStyle("Plain");

  TChain *t = new TChain("t");
  TChain *t2 = new TChain("t");

  //t->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1882_physics.root");
  t->Add("output_1333_3232.root");
  t2->Add("../../BiPo/analysis/SecondDataPeriod/Reprocessed_110817/bga_all_1569_1882_physics.root");

  int fncl;
  int fnr;
  int fne;
  double fX;
  double fY;
  double fZ;
  double fcsc;
  double fcscCorr;
  double ftsc;
  int ftrigsec;
  int ftrigsub;
  int ftrigoff;
  double fAlphaEnergy;

  t->SetBranchAddress("fncl",&fncl);
  t->SetBranchAddress("fnr",&fnr);
  t->SetBranchAddress("fne",&fne);
  t->SetBranchAddress("fX",&fX);
  t->SetBranchAddress("fY",&fY);
  t->SetBranchAddress("fZ",&fZ);
  t->SetBranchAddress("fcsc",&fcsc);
  t->SetBranchAddress("fcscCorr",&fcscCorr);
  t->SetBranchAddress("ftsc",&ftsc);
  t->SetBranchAddress("ftrigsec",&ftrigsec);
  t->SetBranchAddress("ftrigsub",&ftrigsub);
  t->SetBranchAddress("ftrigoff",&ftrigoff);
  t->SetBranchAddress("fAlphaEnergy",&fAlphaEnergy);

  TH1F *hAlphaEnergySpectrum = new TH1F("hAlphaEnergySpectrum","Alpha energy spectrum",60,3000,10000);
  TH1F *hAlphaE1 = new TH1F("hAlphaE1","Energy of first alpha in coincidence",60,3000,10000);
  TH1F *hAlphaE2 = new TH1F("hAlphaE2","Energy of second alpha in coincidence",60,3000,10000);
  TH1F *hAlphaE214Po = new TH1F("hAlphaE214Po","Energy of 214Po",60,3000,10000);
  TH1F *hAlphaDT = new TH1F("hAlphaDT","Alpha coincidence dt",100,0,300000000);
  TH1F *hdR = new TH1F("hdR","dR",30,0,300);
  TH1F *hV = new TH1F("hV","Drift Velocity",60,-2.0,2.0);
  TH2F *hdtdr = new TH2F("hdtdr","dt vs dR",50,-20,200,100,0,1000);
  TH1F *hdZ = new TH1F("hdZ","dZ",50,0,200);
  TH3F *h3D = new TH3F("h3D","Rn222 positions",50,-200,200,50,-200,200,50,-200,200);

  TH1F *hCSC = new TH1F("hCSC","Scintillation spectrum",100,0,80000);
  TH1F *hCSC1 = new TH1F("hCSC1","Scintillation spectrum first alpha",100,0,80000);
  TH1F *hCSC2 = new TH1F("hCSC2","Scintillation spectrum second alpha",100,0,80000);

  hCSC->SetMarkerStyle(20);
  hCSC->SetMarkerSize(0.6);

  hAlphaEnergySpectrum->SetMarkerStyle(20);
  hAlphaEnergySpectrum->SetMarkerSize(0.6);

  TH2F *hXY = new TH2F("hXY","Position of coincidence alphas",40,-200,200,40,-200,200); 

  TH1F *h214Po_fid = new TH1F("h214Po_fid","214Po spectrum (fiducial)",100,0,80000);

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
  //double a1 = -11.16;
  //double a2 = 29.46;
  //double b1 = 21112.1;
  //double b2 = 21200.6;
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
     hAlphaE214Po->Fill((csc_corr + 502.884) / 4.06702);
  }

  h214Po_fid->SetLineColor(kGreen);
  h214Po_fid->SetLineStyle(2);
  //h214Po_fid->SetFillColor(kGreen);
  //h214Po_fid->SetFillStyle(3004);

  hAlphaE214Po->SetLineColor(kGreen);
  hAlphaE214Po->SetLineStyle(2);

  //t2->Draw("fcsc2>>h214Po_fid","TMath::Sqrt(fx2*fx2 + fy2*fy2) < 163 && ((fz2 > -172 && fz2 < -20) || (fz2 > 20 && fz2 < 172))");

  hAlphaE1->SetLineColor(kRed);
  hAlphaE1->SetLineStyle(2);
  //hAlphaE1->SetFillColor(kRed);
  //hAlphaE1->SetFillStyle(3004);

  hAlphaE2->SetLineColor(kBlue);
  hAlphaE2->SetLineStyle(2);
  //hAlphaE2->SetFillColor(kBlue);
  //hAlphaE2->SetFillStyle(3005);

  hCSC1->SetLineColor(kRed);
  hCSC1->SetLineStyle(2);
  //hCSC1->SetFillColor(kRed);
  //hCSC1->SetFillStyle(3004);

  hCSC2->SetLineColor(kBlue);
  hCSC2->SetLineStyle(2);
  //hCSC2->SetFillColor(kBlue);
  //hCSC2->SetFillStyle(3005);

  // arrays holding the alpha charge cluster object and time information
  int nr[1000000];
  int ne[1000000];
  double AlphaX[1000000];
  double AlphaY[1000000];
  double AlphaZ[1000000];
  double AlphaTsc[1000000];
  double AlphaCSC[1000000];
  int AlphaTrigSec[1000000];
  int AlphaTrigSub[1000000];
  int AlphaTrigOff[1000000];
  double AlphaEnergy[1000000];

  // initialize arrays
  for (int i = 0; i < 1000000; i++) {
     nr[i] = 0;
     ne[i] = 0;
     AlphaX[i] = 0.0;
     AlphaY[i] = 0.0;
     AlphaZ[i] = 0.0;
     AlphaTsc[i] = 0.0;
     AlphaCSC[i] = 0.0;
     AlphaTrigSec[i] = 0;
     AlphaTrigSub[i] = 0;
     AlphaTrigOff[i] = 0;
     AlphaEnergy[i] = 0.0;
  }

  int nentries = 0;
  // loop over all alpha events
  for (int i = 0; i < t->GetEntries(); i++) {
     t->GetEntry(i);

     if (fncl != 1) {continue;}

     nr[nentries] = fnr;
     ne[nentries] = fne;
     AlphaX[nentries] = fX;
     AlphaY[nentries] = fY;
     AlphaZ[nentries] = fZ;
     AlphaTsc[nentries] = ftsc;
     AlphaCSC[nentries] = fcscCorr;
     AlphaTrigSec[nentries] = ftrigsec;
     AlphaTrigSub[nentries] = ftrigsub;
     AlphaTrigOff[nentries] = ftrigoff;
     AlphaEnergy[nentries] = fAlphaEnergy;

     //hAlphaSpectrum->Fill(fAlphaEnergy);

     if (TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -172 && fZ < -20) || (fZ > 20 && fZ < 172))) {hCSC->Fill(fcscCorr); hAlphaEnergySpectrum->Fill((fcscCorr + 502.884) / 4.06702);}

     nentries++;
  }

  // look for coincidences
  double RCut = 15; // position cut (28 mm)
  double TCut = 300; // time cut (5 min)

  int nCoincidences = 0;
  for (int i = 0; i < nentries-1; i++) {
     for (int k = i+1; k < nentries; k++) {
        if (TMath::Abs(AlphaTrigSec[i] - AlphaTrigSec[k]) > TCut) {continue;}
        int nr1 = 0;
        int nr2 = 0;
        int ne1 = 0;
        int ne2 = 0;
        double dtAlpha = 0.0;
        double AlphaEnergy1 = 0.0;
        double AlphaEnergy2 = 0.0;
        double AlphaCSC1 = 0.0;
        double AlphaCSC2 = 0.0;
        double AlphaX1 = 0.0;
        double AlphaX2 = 0.0;
        double AlphaY1 = 0.0;
        double AlphaY2 = 0.0;
        double AlphaZ1 = 0.0;
        double AlphaZ2 = 0.0;
        double tMicro1 = AlphaTsc[i] / 1000.0 - double(AlphaTrigOff[i]) + double(AlphaTrigSub[i]);
        double tMicro2 = AlphaTsc[k] / 1000.0 - double(AlphaTrigOff[k]) + double(AlphaTrigSub[k]);

        if (AlphaTrigSec[k] > AlphaTrigSec[i]) {
           dtAlpha = (double(AlphaTrigSec[k]) - double(AlphaTrigSec[i]))*1000000.0 + (tMicro2 - tMicro1);
           nr1 = nr[i];
           nr2 = nr[k];
           ne1 = ne[i];
           ne2 = ne[k];
           AlphaX1 = AlphaX[i];
           AlphaX2 = AlphaX[k];
           AlphaY1 = AlphaY[i];
           AlphaY2 = AlphaY[k];
           AlphaZ1 = AlphaZ[i];
           AlphaZ2 = AlphaZ[k];
           AlphaEnergy1 = AlphaEnergy[i];
           AlphaEnergy2 = AlphaEnergy[k];
           AlphaCSC1 = AlphaCSC[i];
           AlphaCSC2 = AlphaCSC[k];
        }
        if (AlphaTrigSec[k] == AlphaTrigSec[i]) {
           //cout << "dt = 0" << endl;
           //continue;
           if (tMicro2 > tMicro1) {
              dtAlpha = tMicro2 - tMicro1;
              nr1 = nr[i];
              nr2 = nr[k];
              ne1 = ne[i];
              ne2 = ne[k];
              AlphaX1 = AlphaX[i];
              AlphaX2 = AlphaX[k];
              AlphaY1 = AlphaY[i];
              AlphaY2 = AlphaY[k];
              AlphaZ1 = AlphaZ[i];
              AlphaZ2 = AlphaZ[k];
              AlphaEnergy1 = AlphaEnergy[i];
              AlphaEnergy2 = AlphaEnergy[k];
              AlphaCSC1 = AlphaCSC[i];
              AlphaCSC2 = AlphaCSC[k];
           }
           else {
              dtAlpha = tMicro1 - tMicro2;
              nr1 = nr[k];
              nr2 = nr[i];
              ne1 = ne[k];
              ne2 = ne[i];
              AlphaX1 = AlphaX[k];
              AlphaX2 = AlphaX[i];
              AlphaY1 = AlphaY[k];
              AlphaY2 = AlphaY[i];
              AlphaZ1 = AlphaZ[k];
              AlphaZ2 = AlphaZ[i];
              AlphaEnergy1 = AlphaEnergy[k];
              AlphaEnergy2 = AlphaEnergy[i];
              AlphaCSC1 = AlphaCSC[k];
              AlphaCSC2 = AlphaCSC[i];
           }
        }
        if (AlphaTrigSec[k] < AlphaTrigSec[i]) {
           dtAlpha = (double(AlphaTrigSec[i]) - double(AlphaTrigSec[k]))*1000000.0 + (tMicro1 - tMicro2);
           nr1 = nr[k];
           nr2 = nr[i];
           ne1 = ne[k];
           ne2 = ne[i];
           AlphaX1 = AlphaX[k];
           AlphaX2 = AlphaX[i];
           AlphaY1 = AlphaY[k];
           AlphaY2 = AlphaY[i];
           AlphaZ1 = AlphaZ[k];
           AlphaZ2 = AlphaZ[i];
           AlphaEnergy1 = AlphaEnergy[k];
           AlphaEnergy2 = AlphaEnergy[i];
           AlphaCSC1 = AlphaCSC[k];
           AlphaCSC2 = AlphaCSC[i];
        }
       
        if (nr1 == nr2 && ne1 == ne2 && AlphaCSC1 == AlphaCSC2) {continue;}

        double dR = TMath::Sqrt((AlphaX1 - AlphaX2)*(AlphaX1 - AlphaX2) + (AlphaY1 - AlphaY2)*(AlphaY1 - AlphaY2));

        if (TMath::Sqrt(AlphaX1*AlphaX1 + AlphaY1*AlphaY1) < 163 && ((AlphaZ1 > -172 && AlphaZ1 < -20) || (AlphaZ1 > 20 && AlphaZ1 < 172)) && TMath::Sqrt(AlphaX2*AlphaX2 + AlphaY2*AlphaY2) < 163 && ((AlphaZ2 > -172 && AlphaZ2 < -20) || (AlphaZ2 > 20 && AlphaZ2 < 172))) {hdR->Fill(dR);}

        if (dR > RCut) {continue;}
        if (AlphaZ1 < 0.0 && AlphaZ2 > 0.0 || AlphaZ1 > 0.0 && AlphaZ2 < 0.0) {continue;}
        //if (TMath::Abs(AlphaZ1) > 70) {continue;}
        //if (AlphaEnergy1 < 4000 || AlphaEnergy2 < 4000) {continue;}

        std::cout << "(" << nCoincidences << ") csc1 = " << AlphaCSC1 << "\tcsc2 = " << AlphaCSC2 << "\tdR = " << dR << "\tdt = " << dtAlpha << "\tdZ = " << TMath::Abs(AlphaZ1) - TMath::Abs(AlphaZ2) << "\tdt = " << dtAlpha/1000000.0 << "\tfnr1 = " << nr1 << "\tfnr2 = " << nr2 << "\tfne1 = " << ne1 << "\tfne2 = " << ne2 << std::endl;
        hAlphaDT->Fill(dtAlpha);
        //hAlpha1Coinc->Fill(AlphaEnergy1);
        //hAlpha2Coinc->Fill(AlphaEnergy2);
        hXY->Fill(AlphaX2,AlphaY2);

        // apply fiducial cut
        //if (TMath::Sqrt(AlphaX1*AlphaX1 + AlphaY1*AlphaY1) < 163 && ((AlphaZ1 > -172 && AlphaZ1 < -20) || (AlphaZ1 > 20 && AlphaZ1 < 172))) {hCSC1->Fill(AlphaCSC1);}
        //if (TMath::Sqrt(AlphaX2*AlphaX2 + AlphaY2*AlphaY2) < 163 && ((AlphaZ2 > -172 && AlphaZ2 < -20) || (AlphaZ2 > 20 && AlphaZ2 < 172))) {hCSC2->Fill(AlphaCSC2);}
        if (TMath::Sqrt(AlphaX1*AlphaX1 + AlphaY1*AlphaY1) < 163 && ((AlphaZ1 > -172 && AlphaZ1 < -20) || (AlphaZ1 > 20 && AlphaZ1 < 172)) && TMath::Sqrt(AlphaX2*AlphaX2 + AlphaY2*AlphaY2) < 163 && ((AlphaZ2 > -172 && AlphaZ2 < -20) || (AlphaZ2 > 20 && AlphaZ2 < 172))) {
           hCSC1->Fill(AlphaCSC1);
           hCSC2->Fill(AlphaCSC2);
           hAlphaE1->Fill((AlphaCSC1 + 502.884) / 4.06702);
           hAlphaE2->Fill((AlphaCSC2 + 502.884) / 4.06702);
           if (dR == 0.0) {
              hV->Fill((TMath::Abs(AlphaZ1) - TMath::Abs(AlphaZ2)) /dtAlpha * 1000000.0);
              hdtdr->Fill(TMath::Abs(AlphaZ1) - TMath::Abs(AlphaZ2), dtAlpha/1000000.0);
              if (dtAlpha/1000000.0 > (AlphaZ1 - AlphaZ2) * 1.1) {hdZ->Fill(AlphaZ1 - AlphaZ2);}
              h3D->Fill(AlphaX1,AlphaY1,AlphaZ1);
           }
        }

        nCoincidences++;
     }
  }

  //hCSC1->Scale(2.1278);
  //hCSC2->Scale(1.3325);
  hCSC1->Scale(3.55);
  hCSC2->Scale(2.45);

  hAlphaE1->Scale(3.55);
  hAlphaE2->Scale(2.45);

  TH1F *hCSC_sum = (TH1F*)hCSC1->Clone("hCSC_sum");
  hCSC_sum->Add(hCSC2);
  hCSC_sum->Add(h214Po_fid);
  hCSC_sum->SetLineColor(kGray);
  hCSC_sum->SetLineStyle(1);
  hCSC_sum->SetFillColor(kGray);
  hCSC_sum->SetFillStyle(1001);

  TH1F *hE_sum = (TH1F*)hAlphaE1->Clone("hE_sum");
  hE_sum->Add(hAlphaE2);
  hE_sum->Add(hAlphaE214Po);
  hE_sum->SetLineColor(kGray);
  hE_sum->SetLineStyle(1);
  hE_sum->SetFillColor(kGray);
  hE_sum->SetFillStyle(1001);

  TH1F *hCSC_sub = (TH1F*)hCSC->Clone("hCSC_sub");
  hCSC_sub->Add(hCSC_sum,-1);
  //hCSC_sub->Add(h214Po_fid,-1);

  for (int i = 0; i < 38; i++) {hCSC_sub->SetBinContent(i+1,0);}

  double dt = 31.4025;
  double integralErr;
  double integral = hCSC1->IntegralAndError(0,80,integralErr);
  cout << "Rn222: " << integral << " +- " << integralErr << "  (" << integral / dt << " +- " << integralErr / dt << " d-1)" << endl;

  TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
  l->AddEntry(hCSC,"Scintillation spectrum");
  l->AddEntry(hCSC1,"Rn222 (#alpha-#alpha coincidence)");
  l->AddEntry(hCSC2,"Po218 (#alpha-#alpha coincidence)");
  l->AddEntry(h214Po_fid,"Po214 (Bi-Po)");

  TCanvas *c1 = new TCanvas("c1","Scintillation spectrum");
  hCSC->Draw("EP");
  hCSC_sum->Draw();
  hCSC1->Draw("same");
  hCSC2->Draw("same");
  h214Po_fid->Draw("same");
  hCSC->Draw("EPsame");
  l->Draw("same");

  TCanvas *c2 = new TCanvas("c2","Energy spectrum");
  hAlphaEnergySpectrum->Draw("EP");
  hE_sum->Draw();
  hAlphaE1->Draw("same");
  hAlphaE2->Draw("same");
  hAlphaE214Po->Draw("same");
  hAlphaEnergySpectrum->Draw("EPsame");
  l->Draw("same");

  TCanvas *c3 = new TCanvas("c3","Alpha coincidence dt");
  hAlphaDT->Draw();

  TCanvas *c4 = new TCanvas("c4","XY distribution");
  hdR->Draw();

  TCanvas *c5 = new TCanvas("c5","Scintillation spectrum sub");
  hV->Draw();

  TCanvas *c6 = new TCanvas("c6","dt vs dR");
  hdZ->Draw();
  
  TCanvas *c7 = new TCanvas("c7","dt vs dz");
  hdtdr->Draw();
  
  TCanvas *c8 = new TCanvas("c8","Rn222 positions");
  h3D->Draw();

  return;
}

