void ChargeCollectionEfficiency2()
{
  TChain *t = new TChain("t");
  t->Add("../analysis/SecondDataPeriod/output_ncl_lt1_1569_1782.root");

  TH1F *hNCL_EQ1 = new TH1F("hNCL_EQ1","ncl = 1",100,0,40000);
  TH1F *hNCL_LT1 = new TH1F("hNCL_LT1","ncl <= 1",100,0,40000);
  TH1F *hCSC_EQ1_C = new TH1F("hCSC_EQ1_C","cathode/anode",100,0,40000);
  TH1F *hCSC_EQ1_B = new TH1F("hCSC_EQ1_B","bulk",100,0,40000);
  TH1F *hCSC_EQ1_S = new TH1F("hCSC_EQ1_S","all",100,0,40000);

  hCSC_EQ1_C->SetLineColor(kRed);
  hCSC_EQ1_S->SetLineColor(kGray);
  hCSC_EQ1_S->SetFillColor(kGray);

  t->Draw("fcsc>>hNCL_EQ1","fncl == 1");
  t->Draw("fcsc>>hNCL_LT1");
  t->Draw("fcsc>>hCSC_EQ1_S","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163");
  t->Draw("fcsc>>hCSC_EQ1_B","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -182 && fZ < -10) || (fZ > 10 && fZ < 182))");
  t->Draw("fcsc>>hCSC_EQ1_C","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && (fZ < -182 || fZ > 182 || (fZ > -10 && fZ < 10))");

  TGraphErrors *grEff = new TGraphErrors(100);
  for (int i = 0; i < 100; i++) {
     double BinContNCL_EQ1 = hNCL_EQ1->GetBinContent(i+1);
     double BinContNCL_LT1 = hNCL_LT1->GetBinContent(i+1);
     grEff->SetPoint(i,(i+0.5)*400,BinContNCL_EQ1 / BinContNCL_LT1);

     double err = TMath::Sqrt(1/(BinContNCL_LT1*BinContNCL_LT1) * TMath::Sqrt(BinContNCL_EQ1) * TMath::Sqrt(BinContNCL_EQ1) + (BinContNCL_EQ1 / (BinContNCL_LT1*BinContNCL_LT1)) * (BinContNCL_EQ1 / (BinContNCL_LT1*BinContNCL_LT1)) * TMath::Sqrt(BinContNCL_LT1) * TMath::Sqrt(BinContNCL_LT1));
     grEff->SetPointError(i,0,err);
  }

  TF1 *fitEff = new TF1("fitEff","0.5*(TMath::Erf((x-[0])/[1]) + 1.0)",19000,35000);
  fitEff->SetParameters(20000,1000);

  grEff->Fit("fitEff","rn");

  TF1 *fitEff2 = new TF1("fitEff2","0.5*(TMath::Erf((x-[0])/[1]) + 1.0)",0,35000);
  double *params = fitEff->GetParameters();
  fitEff2->SetParameters(params);

  TLegend *l1 = new TLegend(0.0,0.0,0.2,0.2);
  l1->AddEntry(hCSC_EQ1_B,"bulk");
  l1->AddEntry(hCSC_EQ1_C,"anode/cathode");
  l1->AddEntry(hCSC_EQ1_S,"all");

  TCanvas *c1 = new TCanvas("csc anode/cathode and bulk");
  hCSC_EQ1_S->Draw();
  hCSC_EQ1_B->Draw("same");
  hCSC_EQ1_C->Draw("same");
  l1->Draw("same");

  TCanvas *c2 = new TCanvas("c1sc/c2sc vs c1sc");
  t->Draw("TMath::Log10(fc1sc/fc2sc):fc1sc","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && (fZ < -187 || fZ > 187 || (fZ > -5 && fZ < 5))");
  t->Draw("TMath::Log10(fc1sc/fc2sc):fc1sc","fncl == 1 && TMath::Sqrt(fX*fX + fY*fY) < 163 && ((fZ > -187 && fZ < -5) || (fZ > 5 && fZ < 187))","same");

  TCanvas *c3 = new TCanvas();
  grEff->Draw("AP");
  fitEff2->Draw("same");

  return;
}
