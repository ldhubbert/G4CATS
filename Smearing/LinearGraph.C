#include <TMath.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TStyle.h>
void maincode()
{
	//Code From 100MeV Smear
	//Gaussian and Filling Histogram
	TFile a("~/Vincent/G4CATS/Out/B4_100MeV.root");

	TF1 *f1 = new TF1("f1", "gaus(0)", -3, 3);
	f1->SetParameter(0, (1/((0.07*TMath::Sqrt(100))*(TMath::Sqrt(2*TMath::Pi())))));
	f1->SetParameter(1, 0);
	f1->SetParameter(2, 0.07*TMath::Sqrt(100));

	TTreeReader r1("B4", &a);

	TTreeReaderValue<Double_t> Ecore(r1, "Ecore");
	TTreeReaderValue<Double_t> Eann1(r1, "Eann1");
	TTreeReaderValue<Double_t> Eann2(r1, "Eann2");
	TTreeReaderValue<Double_t> Eann3(r1, "Eann3");
	TTreeReaderValue<Double_t> Eann4(r1, "Eann4");
	TTreeReaderValue<Double_t> Eann5(r1, "Eann5");
	TTreeReaderValue<Double_t> Eann6(r1, "Eann6");

	TH1F *h1 = new TH1F("", "", 200, 85, 105);

	while (r1.Next())
	{
		h1->Fill(*Ecore + *Eann1 + *Eann2 + *Eann3 + *Eann4 + *Eann5 + *Eann6 + f1->GetRandom());
	}
	
	//100MeV FWHM Code
	Double_t BinWithMostCounts = h1->GetMaximumBin();
	Double_t MaxYValue = h1->GetBinContent(BinWithMostCounts);
	Double_t HalfMaxYValue = MaxYValue/2;

	Double_t FWHMLeftXValue = 0;
	int binA = 0;
	h1->GetBinWithContent(HalfMaxYValue, binA, 0, BinWithMostCounts, 0);

	if (binA == 0)
	{
		Double_t Lower2 = h1->FindFirstBinAbove(HalfMaxYValue, 1, 0, BinWithMostCounts);
		Double_t Lower2Contents = h1->GetBinContent(Lower2);

		Double_t Lower1 = Lower2 - 1;
		Double_t Lower1Contents = h1->GetBinContent(Lower1);

		Double_t LowerBinFraction = Lower1Contents/Lower2Contents;
		Double_t CenterLower1 = h1->GetBinCenter(Lower1);
		Double_t CenterLower2 = h1->GetBinCenter(Lower2);
		Double_t LowerBinWidth = CenterLower2 - CenterLower1;
		FWHMLeftXValue = CenterLower1 + (LowerBinWidth)*(LowerBinFraction);
	}

	else if (binA != 0)
	{
		FWHMLeftXValue = h1->GetBinCenter(binA);
	}

	Double_t FWHMRightXValue = 0;
	int binB = 0;
	h1->GetBinWithContent(HalfMaxYValue, binB, BinWithMostCounts, 200, 0);
	if (binB == 0)
	{
		Double_t Upper1 = h1->FindLastBinAbove(HalfMaxYValue, 1, BinWithMostCounts, 200);
		Double_t Upper1Contents = h1->GetBinContent(Upper1);

		Double_t Upper2 = Upper1 + 1;
		Double_t Upper2Contents = h1->GetBinContent(Upper2);

		Double_t UpperBinFraction = Upper2Contents/Upper1Contents;
		Double_t CenterUpper1 = h1->GetBinCenter(Upper1);
		Double_t CenterUpper2 = h1->GetBinCenter(Upper2);
		Double_t UpperBinWidth = CenterUpper2 - CenterUpper1;
		FWHMRightXValue = CenterUpper2 + (UpperBinWidth)*(UpperBinFraction);
	}

	else if (binB != 0)
	{
		FWHMRightXValue = h1->GetBinCenter(binB);
	}

	Double_t FinalWidth = FWHMRightXValue - FWHMLeftXValue;
	Double_t FWHM = (FinalWidth/100)*100;
	cout << FWHM << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////

	//Code From 200MeV Smear
	//Gaussian and Filling Histogram
	TFile b("~/Vincent/G4CATS/Out/B4_200MeV.root");

	TF1 *f2 = new TF1("f2", "gaus(0)", -5, 5);
	f2->SetParameter(0, (1/((0.09*TMath::Sqrt(200))*(TMath::Sqrt(2*TMath::Pi())))));
	f2->SetParameter(1, 0);
	f2->SetParameter(2, 0.09*TMath::Sqrt(200));

	TTreeReader r2("B4", &b);

	TTreeReaderValue<Double_t> Ecore_2(r2, "Ecore");
	TTreeReaderValue<Double_t> Eann1_2(r2, "Eann1");
	TTreeReaderValue<Double_t> Eann2_2(r2, "Eann2");
	TTreeReaderValue<Double_t> Eann3_2(r2, "Eann3");
	TTreeReaderValue<Double_t> Eann4_2(r2, "Eann4");
	TTreeReaderValue<Double_t> Eann5_2(r2, "Eann5");
	TTreeReaderValue<Double_t> Eann6_2(r2, "Eann6");

	TH1F *h2 = new TH1F("", "", 300, 180, 210);

	while (r2.Next())
	{
		h2->Fill(*Ecore_2 + *Eann1_2 + *Eann2_2 + *Eann3_2 + *Eann4_2 + *Eann5_2 + *Eann6_2 + f2->GetRandom());
	}
	
	//200MeV FWHM Code
	Double_t BinWithMostCounts_2 = h2->GetMaximumBin();
	Double_t MaxYValue_2 = h2->GetBinContent(BinWithMostCounts_2);
	Double_t HalfMaxYValue_2 = MaxYValue_2/2;

	Double_t FWHMLeftXValue_2 = 0;
	int binA_2 = 0;
	h2->GetBinWithContent(HalfMaxYValue_2, binA_2, 0, BinWithMostCounts_2, 0);

	if (binA_2 == 0)
	{
		Double_t Lower2_2 = h2->FindFirstBinAbove(HalfMaxYValue_2, 1, 0, BinWithMostCounts_2);
		Double_t Lower2Contents_2 = h2->GetBinContent(Lower2_2);

		Double_t Lower1_2 = Lower2_2 - 1;
		Double_t Lower1Contents_2 = h2->GetBinContent(Lower1_2);

		Double_t LowerBinFraction_2 = Lower1Contents_2/Lower2Contents_2;
		Double_t CenterLower1_2 = h2->GetBinCenter(Lower1_2);
		Double_t CenterLower2_2 = h2->GetBinCenter(Lower2_2);
		Double_t LowerBinWidth_2 = CenterLower2_2 - CenterLower1_2;
		FWHMLeftXValue_2 = CenterLower1_2 + (LowerBinWidth_2)*(LowerBinFraction_2);
	}

	else if (binA_2 != 0)
	{
		FWHMLeftXValue_2 = h2->GetBinCenter(binA_2);
	}

	Double_t FWHMRightXValue_2 = 0;
	int binB_2 = 0;
	h2->GetBinWithContent(HalfMaxYValue_2, binB_2, BinWithMostCounts_2, 300, 0);
	if (binB_2 == 0)
	{
		Double_t Upper1_2 = h2->FindLastBinAbove(HalfMaxYValue_2, 1, BinWithMostCounts_2, 300);
		Double_t Upper1Contents_2 = h2->GetBinContent(Upper1_2);

		Double_t Upper2_2 = Upper1_2 + 1;
		Double_t Upper2Contents_2 = h2->GetBinContent(Upper2_2);

		Double_t UpperBinFraction_2 = Upper2Contents_2/Upper1Contents_2;
		Double_t CenterUpper1_2 = h2->GetBinCenter(Upper1_2);
		Double_t CenterUpper2_2 = h2->GetBinCenter(Upper2_2);
		Double_t UpperBinWidth_2 = CenterUpper2_2 - CenterUpper1_2;
		FWHMRightXValue_2 = CenterUpper2_2 + (UpperBinWidth_2)*(UpperBinFraction_2);
	}

	else if (binB_2 != 0)
	{
		FWHMRightXValue_2 = h2->GetBinCenter(binB_2);
	}

	Double_t FinalWidth_2 = FWHMRightXValue_2 - FWHMLeftXValue_2;
	Double_t FWHM_2 = (FinalWidth_2/200)*100;
	cout << FWHM_2 << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////

	//Code From 300MeV Smear
	//Gaussian and Filling Histogram
	TFile c("~/Vincent/G4CATS/Out/B4_300MeV.root");

	TF1 *f3 = new TF1("f3", "gaus(0)", -5, 5);
	f3->SetParameter(0, (1/((0.11*TMath::Sqrt(300))*(TMath::Sqrt(2*TMath::Pi())))));
	f3->SetParameter(1, 0);
	f3->SetParameter(2, 0.11*TMath::Sqrt(300));

	TTreeReader r3("B4", &c);

	TTreeReaderValue<Double_t> Ecore_3(r3, "Ecore");
	TTreeReaderValue<Double_t> Eann1_3(r3, "Eann1");
	TTreeReaderValue<Double_t> Eann2_3(r3, "Eann2");
	TTreeReaderValue<Double_t> Eann3_3(r3, "Eann3");
	TTreeReaderValue<Double_t> Eann4_3(r3, "Eann4");
	TTreeReaderValue<Double_t> Eann5_3(r3, "Eann5");
	TTreeReaderValue<Double_t> Eann6_3(r3, "Eann6");

	TH1F *h3 = new TH1F("", "", 400, 270, 310);

	while (r3.Next())
	{
		h3->Fill(*Ecore_3 + *Eann1_3 + *Eann2_3 + *Eann3_3 + *Eann4_3 + *Eann5_3 + *Eann6_3 + f3->GetRandom());
	}
	
	//300MeV FWHM Code
	Double_t BinWithMostCounts_3 = h3->GetMaximumBin();
	Double_t MaxYValue_3 = h3->GetBinContent(BinWithMostCounts_3);
	Double_t HalfMaxYValue_3 = MaxYValue_3/2;

	Double_t FWHMLeftXValue_3 = 0;
	int binA_3 = 0;
	h3->GetBinWithContent(HalfMaxYValue_3, binA_3, 0, BinWithMostCounts_3, 0);

	if (binA_3 == 0)
	{
		Double_t Lower2_3 = h3->FindFirstBinAbove(HalfMaxYValue_3, 1, 0, BinWithMostCounts_3);
		Double_t Lower2Contents_3 = h3->GetBinContent(Lower2_3);

		Double_t Lower1_3 = Lower2_3 - 1;
		Double_t Lower1Contents_3 = h1->GetBinContent(Lower1_3);

		Double_t LowerBinFraction_3 = Lower1Contents_3/Lower2Contents_3;
		Double_t CenterLower1_3 = h3->GetBinCenter(Lower1_3);
		Double_t CenterLower2_3 = h3->GetBinCenter(Lower2_3);
		Double_t LowerBinWidth_3 = CenterLower2_3 - CenterLower1_3;
		FWHMLeftXValue_3 = CenterLower1_3 + (LowerBinWidth_3)*(LowerBinFraction_3);
	}

	else if (binA_3 != 0)
	{
		FWHMLeftXValue_3 = h3->GetBinCenter(binA_3);
	}

	Double_t FWHMRightXValue_3 = 0;
	int binB_3 = 0;
	h3->GetBinWithContent(HalfMaxYValue_3, binB_3, BinWithMostCounts_3, 400, 0);
	if (binB_3 == 0)
	{
		Double_t Upper1_3 = h3->FindLastBinAbove(HalfMaxYValue_3, 1, BinWithMostCounts_3, 400);
		Double_t Upper1Contents_3 = h3->GetBinContent(Upper1_3);

		Double_t Upper2_3 = Upper1_3 + 1;
		Double_t Upper2Contents_3 = h3->GetBinContent(Upper2_3);

		Double_t UpperBinFraction_3 = Upper2Contents_3/Upper1Contents_3;
		Double_t CenterUpper1_3 = h3->GetBinCenter(Upper1_3);
		Double_t CenterUpper2_3 = h3->GetBinCenter(Upper2_3);
		Double_t UpperBinWidth_3 = CenterUpper2_3 - CenterUpper1_3;
		FWHMRightXValue_3 = CenterUpper2_3 + (UpperBinWidth_3)*(UpperBinFraction_3);
	}

	else if (binB_3 != 0)
	{
		FWHMRightXValue_3 = h3->GetBinCenter(binB_3);
	}

	Double_t FinalWidth_3 = FWHMRightXValue_3 - FWHMLeftXValue_3;
	Double_t FWHM_3 = (FinalWidth_3/300)*100;
	cout << FWHM_3 << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Code From 400MeV Smear
	//Gaussian and Filling Histogram
	TFile d("~/Vincent/G4CATS/Out/B4_400MeV.root");

	TF1 *f4 = new TF1("f4", "gaus(0)", -10, 10);
	f4->SetParameter(0, (1/((0.13*TMath::Sqrt(400))*(TMath::Sqrt(2*TMath::Pi())))));
	f4->SetParameter(1, 0);
	f4->SetParameter(2, 0.13*TMath::Sqrt(400));

	TTreeReader r4("B4", &d);

	TTreeReaderValue<Double_t> Ecore_4(r4, "Ecore");
	TTreeReaderValue<Double_t> Eann1_4(r4, "Eann1");
	TTreeReaderValue<Double_t> Eann2_4(r4, "Eann2");
	TTreeReaderValue<Double_t> Eann3_4(r4, "Eann3");
	TTreeReaderValue<Double_t> Eann4_4(r4, "Eann4");
	TTreeReaderValue<Double_t> Eann5_4(r4, "Eann5");
	TTreeReaderValue<Double_t> Eann6_4(r4, "Eann6");

	TH1F *h4 = new TH1F("", "", 500, 365, 415);

	while (r4.Next())
	{
		h4->Fill(*Ecore_4 + *Eann1_4 + *Eann2_4 + *Eann3_4 + *Eann4_4 + *Eann5_4 + *Eann6_4 + f4->GetRandom());
	}
	
	//400MeV FWHM Code
	Double_t BinWithMostCounts_4 = h4->GetMaximumBin();
	Double_t MaxYValue_4 = h4->GetBinContent(BinWithMostCounts_4);
	Double_t HalfMaxYValue_4 = MaxYValue_4/2;

	Double_t FWHMLeftXValue_4 = 0;
	int binA_4 = 0;
	h4->GetBinWithContent(HalfMaxYValue_4, binA_4, 0, BinWithMostCounts_4, 0);

	if (binA_4 == 0)
	{
		Double_t Lower2_4 = h4->FindFirstBinAbove(HalfMaxYValue_4, 1, 0, BinWithMostCounts_4);
		Double_t Lower2Contents_4 = h4->GetBinContent(Lower2_4);

		Double_t Lower1_4 = Lower2_4 - 1;
		Double_t Lower1Contents_4 = h4->GetBinContent(Lower1_4);

		Double_t LowerBinFraction_4 = Lower1Contents_4/Lower2Contents_4;
		Double_t CenterLower1_4 = h4->GetBinCenter(Lower1_4);
		Double_t CenterLower2_4 = h4->GetBinCenter(Lower2_4);
		Double_t LowerBinWidth_4 = CenterLower2_4 - CenterLower1_4;
		FWHMLeftXValue_4 = CenterLower1_4 + (LowerBinWidth_4)*(LowerBinFraction_4);
	}

	else if (binA_4 != 0)
	{
		FWHMLeftXValue_4 = h4->GetBinCenter(binA_4);
	}

	Double_t FWHMRightXValue_4 = 0;
	int binB_4 = 0;
	h4->GetBinWithContent(HalfMaxYValue_4, binB_4, BinWithMostCounts_4, 500, 0);
	if (binB_4 == 0)
	{
		Double_t Upper1_4 = h4->FindLastBinAbove(HalfMaxYValue_4, 1, BinWithMostCounts_4, 500);
		Double_t Upper1Contents_4 = h4->GetBinContent(Upper1_4);

		Double_t Upper2_4 = Upper1_4 + 1;
		Double_t Upper2Contents_4 = h4->GetBinContent(Upper2_4);

		Double_t UpperBinFraction_4 = Upper2Contents_4/Upper1Contents_4;
		Double_t CenterUpper1_4 = h4->GetBinCenter(Upper1_4);
		Double_t CenterUpper2_4 = h4->GetBinCenter(Upper2_4);
		Double_t UpperBinWidth_4 = CenterUpper2_4 - CenterUpper1_4;
		FWHMRightXValue_4 = CenterUpper2_4 + (UpperBinWidth_4)*(UpperBinFraction_4);
	}

	else if (binB_4 != 0)
	{
		FWHMRightXValue_4 = h4->GetBinCenter(binB_4);
	}

	Double_t FinalWidth_4 = FWHMRightXValue_4 - FWHMLeftXValue_4;
	Double_t FWHM_4 = (FinalWidth_4/400)*100;
	cout << FWHM_4 << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Code From 500MeV Smear
	//Gaussian and Filling Histogram
	TFile e("~/Vincent/G4CATS/Out/B4_500MeV.root");

	TF1 *f5 = new TF1("f1", "gaus(0)", -12, 12);
	f5->SetParameter(0, (1/((0.15*TMath::Sqrt(500))*(TMath::Sqrt(2*TMath::Pi())))));
	f5->SetParameter(1, 0);
	f5->SetParameter(2, 0.15*TMath::Sqrt(500));

	TTreeReader r5("B4", &e);

	TTreeReaderValue<Double_t> Ecore_5(r5, "Ecore");
	TTreeReaderValue<Double_t> Eann1_5(r5, "Eann1");
	TTreeReaderValue<Double_t> Eann2_5(r5, "Eann2");
	TTreeReaderValue<Double_t> Eann3_5(r5, "Eann3");
	TTreeReaderValue<Double_t> Eann4_5(r5, "Eann4");
	TTreeReaderValue<Double_t> Eann5_5(r5, "Eann5");
	TTreeReaderValue<Double_t> Eann6_5(r5, "Eann6");

	TH1F *h5 = new TH1F("", "", 600, 455, 515);

	while (r5.Next())
	{
		h5->Fill(*Ecore_5 + *Eann1_5 + *Eann2_5 + *Eann3_5 + *Eann4_5 + *Eann5_5 + *Eann6_5 + f5->GetRandom());
	}
	
	//500MeV FWHM Code
	Double_t BinWithMostCounts_5 = h5->GetMaximumBin();
	Double_t MaxYValue_5 = h5->GetBinContent(BinWithMostCounts_5);
	Double_t HalfMaxYValue_5 = MaxYValue_5/2;

	Double_t FWHMLeftXValue_5 = 0;
	int binA_5 = 0;
	h5->GetBinWithContent(HalfMaxYValue_5, binA_5, 0, BinWithMostCounts_5, 0);

	if (binA_5 == 0)
	{
		Double_t Lower2_5 = h5->FindFirstBinAbove(HalfMaxYValue_5, 1, 0, BinWithMostCounts_5);
		Double_t Lower2Contents_5 = h5->GetBinContent(Lower2_5);

		Double_t Lower1_5 = Lower2_5 - 1;
		Double_t Lower1Contents_5 = h5->GetBinContent(Lower1_5);

		Double_t LowerBinFraction_5 = Lower1Contents_5/Lower2Contents_5;
		Double_t CenterLower1_5 = h5->GetBinCenter(Lower1_5);
		Double_t CenterLower2_5 = h5->GetBinCenter(Lower2_5);
		Double_t LowerBinWidth_5 = CenterLower2_5 - CenterLower1_5;
		FWHMLeftXValue_5 = CenterLower1_5 + (LowerBinWidth_5)*(LowerBinFraction_5);
	}

	else if (binA_5 != 0)
	{
		FWHMLeftXValue_5 = h5->GetBinCenter(binA_5);
	}

	Double_t FWHMRightXValue_5 = 0;
	int binB_5 = 0;
	h5->GetBinWithContent(HalfMaxYValue_5, binB_5, BinWithMostCounts_5, 600, 0);
	if (binB_5 == 0)
	{
		Double_t Upper1_5 = h5->FindLastBinAbove(HalfMaxYValue_5, 1, BinWithMostCounts_5, 600);
		Double_t Upper1Contents_5 = h5->GetBinContent(Upper1_5);

		Double_t Upper2_5 = Upper1_5 + 1;
		Double_t Upper2Contents_5 = h5->GetBinContent(Upper2_5);

		Double_t UpperBinFraction_5 = Upper2Contents_5/Upper1Contents_5;
		Double_t CenterUpper1_5 = h5->GetBinCenter(Upper1_5);
		Double_t CenterUpper2_5 = h5->GetBinCenter(Upper2_5);
		Double_t UpperBinWidth_5 = CenterUpper2_5 - CenterUpper1_5;
		FWHMRightXValue_5 = CenterUpper2_5 + (UpperBinWidth_5)*(UpperBinFraction_5);
	}

	else if (binB_5 != 0)
	{
		FWHMRightXValue_5 = h5->GetBinCenter(binB_5);
	}

	Double_t FinalWidth_5 = FWHMRightXValue_5 - FWHMLeftXValue_5;
	Double_t FWHM_5 = (FinalWidth_5/500)*100;
	cout << FWHM_5 << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//END
	TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
	c1->Divide(1,1);
	c1->cd(1);

	Int_t array_size = 5;
	Double_t x_values[array_size];
	Double_t y_values[array_size];
	Double_t z_values[array_size];

	x_values[0] = 100;
	x_values[1] = 200;
	x_values[2] = 300;
	x_values[3] = 400;
	x_values[4] = 500;

	y_values[0] = FWHM;
	y_values[1] = FWHM_2;
	y_values[2] = FWHM_3;
	y_values[3] = FWHM_4;
	y_values[4] = FWHM_5;

	z_values[0] = 7;
	z_values[1] = 9;
	z_values[2] = 11;
	z_values[3] = 13;
	z_values[4] = 15;

	TGraph2D *g1 = new TGraph2D(array_size, x_values, y_values, z_values);
	g1->SetTitle("Smearing factor for FWHM of various Energy Levels");
	g1->GetXaxis()->SetTitle("Energy (in MeV)");
	g1->GetYaxis()->SetTitle("FWHM (in % of incident energy)");
	g1->GetZaxis()->SetTitle("Percentage");
	gStyle->SetPalette(1);

	g1->Draw();




}
