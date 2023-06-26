#include <TMath.h>
{
  //Open file filled by Geant4 simulation 
  TFile f("~/Vincent/G4CATS/Out/B4_100MeV.root");

  //Create a canvas and divide it into 2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(1,2);

  //Draw smeared histogram in the pad 2
  c1->cd(2);

  //Smearing the data:
  //STEP ONE
  //To smear the data, first we need to make a Gaussian distribution to simulate the detector efficiency.
  //The mean of this Gaussian will be 0MeV, and the standard deviation will be 0.0075, since we want to smear the results by 0.75%.
  //The more we want to smear the results by, the wider the Gaussian curve will be, attributing more error to the histogram results.
  //gaus(0) refers to a Gaussian distribution with parameters as commented below

  TF1 *f1 = new TF1("f1", "gaus", 0, 2);
  //Fraction being raised to power
  f1->SetParameter(0, (1/((0.0075)*(TMath::Sqrt(2*TMath::Pi())))));
  //f1->SetParameter(0, 1);
  //Mean
  f1->SetParameter(1, 1);
  //Standard Deviation
  f1->SetParameter(2, 0.0075);

  //Looking for the branch, "B4", in file f (the 100MeV output file)
  TTreeReader r1("B4", &f);

  //Reading each branch below as defined in B4
  //<Double_t> tells TTreeReaderValue what type of data it's extracting
  TTreeReaderValue<Double_t> Ecore(r1, "Ecore");
  TTreeReaderValue<Double_t> Eann1(r1, "Eann1");
  TTreeReaderValue<Double_t> Eann2(r1, "Eann2");
  TTreeReaderValue<Double_t> Eann3(r1, "Eann3");
  TTreeReaderValue<Double_t> Eann4(r1, "Eann4");
  TTreeReaderValue<Double_t> Eann5(r1, "Eann5");
  TTreeReaderValue<Double_t> Eann6(r1, "Eann6");

  //Create histogram
  TH1F *h1 = new TH1F("Histogram 2", "h1", 200, 85, 105);

  //The while loop goes through each branch and reads entries.
  //So for the first time around the loop, entry 1 is read from each branch (Ecore up to Eann6)
  while (r1.Next())
  {
	h1->Fill((*Ecore + *Eann1 + *Eann2 + *Eann3 + *Eann4 + *Eann5 + *Eann6)*(f1->GetRandom()));
  } 

  TLine *line = new TLine(100, 0, 100, 3650);

  h1->Draw();
  h1->SetTitle("0.75% Gaussian Smear on 100MeV Beam");

  line->Draw();

  //Compare to normal 100MeV Histogram (taken from code in Histogram folder)
  c1->cd(1);
  TH1F *h2 = new TH1F("Histogram 1", "", 200, 85, 105);
  TTree *B4 = (TTree*)f.Get("B4");
  B4->Draw("Ecore+Eann1+Eann2+Eann3+Eann4+Eann5+Eann6>>Histogram 1");
  h2->GetXaxis()->SetTitle("Energy (MeV)");
  h2->GetYaxis()->SetTitle("Counts");
  h2->SetTitle("100MeV Beam Before Smearing");
  h2->Draw();

  //Start of FWHM Section
  Double_t BinWithMostCounts = h1->GetMaximumBin();
  Double_t MaxYValue = h1->GetBinContent(BinWithMostCounts);
  Double_t HalfMaxYValue = MaxYValue/2;
 
  //Finding HalfMaxYValue bin on left-hand side of peak
  Double_t FWHMLeftXValue = 0;
  int binA = 0;
  //This searches for the bin in the range of bin 0 and BinWithMostCounts
  h1->GetBinWithContent(HalfMaxYValue, binA, 0, BinWithMostCounts, 0);
  //The if condition pertains if no bin was found with content of HalfMaxYValue
  if (binA == 0)
  	{
	//Finding the first bin (on the left-hand side of the max peak) that has contents above HalfMaxYValue (Lower2)
	//Not sure why but second parameter has to be 1 when using this search ability
	Double_t Lower2 = h1->FindFirstBinAbove(HalfMaxYValue, 1, 0, BinWithMostCounts);
	Double_t Lower2Contents = h1->GetBinContent(Lower2);

	//Finding the bin below HalfMaxYValue on the left-hand side of the max peak (Lower1)
	Double_t Lower1 = Lower2 - 1;
	Double_t Lower1Contents = h1->GetBinContent(Lower1);

	//Finding an approximation for the bin x-value with HalfMaxYValue on the left-hand side of the max peak
	Double_t LowerBinFraction = Lower1Contents/Lower2Contents;
	Double_t CenterLower1 = h1->GetBinCenter(Lower1);
	Double_t CenterLower2 = h1->GetBinCenter(Lower2);
	Double_t LowerBinWidth = CenterLower2 - CenterLower1;
	FWHMLeftXValue = CenterLower1 + (LowerBinWidth)*(LowerBinFraction);
	}

  //The else if condition pertains if a bin was found with content of HalfMaxYValue
  else if (binA != 0)
  	{
	FWHMLeftXValue = h1->GetBinCenter(binA);
	}

  //Finding HalfMaxYValue bin on right-hand side of peak
  Double_t FWHMRightXValue = 0;
  int binB = 0;
  //This searches for the bin in the range of BinWithMostCounts and bin 200, which is the last bin
  h1->GetBinWithContent(HalfMaxYValue, binB, BinWithMostCounts, 200, 0);
  //The if condition pertains if no bin was found with content of HalfMaxYValue
  if (binB == 0)
  	{
	//Finding the last bin (on the right-hand side of the max peak) that has contents above HalfMaxYValue (Upper1)
	//Not sure why but second parameter has to be 1 when using this search ability
	Double_t Upper1 = h1->FindLastBinAbove(HalfMaxYValue, 1, BinWithMostCounts, 200);
	Double_t Upper1Contents = h1->GetBinContent(Upper1);
	 
	//Finding the bin below HalfMaxYValue on the right-hand side of the max peak (Upper2)
	Double_t Upper2 = Upper1 + 1;
	Double_t Upper2Contents = h1->GetBinContent(Upper2);

	//Finding an approximation for the bin x-value with HalfMaxYValue on the right-hand side of the max peak
 	Double_t UpperBinFraction = Upper2Contents/Upper1Contents;
	Double_t CenterUpper1 = h1->GetBinCenter(Upper1);
	Double_t CenterUpper2 = h1->GetBinCenter(Upper2);
	Double_t UpperBinWidth = CenterUpper2 - CenterUpper1;
	FWHMRightXValue = CenterUpper2 + (UpperBinWidth)*(UpperBinFraction);
	}

  //The else if condition pertains if a bin was found with content of HalfMaxYValue
  else if (binB != 0)
  	{
	FWHMRightXValue = h1->GetBinCenter(binB);
	}

  //Finding the FWHM
  Double_t FinalWidth = FWHMRightXValue - FWHMLeftXValue;
  Double_t FWHM = (FinalWidth/100)*100;
  cout << "FWHM:" << endl;
  cout << FWHM << endl;

  //Display FWHM on bottom pad 
  c1->cd(2);
  TString FWHM_string;
  FWHM_string = Form("FWHM: %lf", FWHM);
  TPaveLabel *a = new TPaveLabel(86,2800,90,3300, FWHM_string);
  a->Draw();
}