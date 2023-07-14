//#include "TF1.h"
//#include "TMath.h"

Double_t MyMethod(Double_t *x, Double_t *par)
{

	Double_t exponent = (-1* (((x[0] - par[2])/par[3]) - par[0]) * (((x[0] - par[2])/par[3]) - par[0])) / (2 * par[1] * par[1]);
	Double_t erf_top_limit = par[4] * (x[0] - par[2]) / (par[3] * TMath::Sqrt(2));
	Double_t function = par[5] * (1 / (par[3] * par[1] * TMath::Sqrt(2 * TMath::Pi()))) * TMath::Exp(exponent) * (1 + TMath::Erf(erf_top_limit));

	return function;

}

void maincode()
{
	//Normal Histogram Stuff
	//gStyle->SetOptTitle(1);
	//gStyle->SetOptStat(0);
	TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
	c1->Divide(1,1);
	c1->cd(1);

	TString filename = "~/Vincent/G4CATS/Out/B4_200MeV.root";
	TFile *f = TFile::Open(filename);

	TH1F *h1 = new TH1F("h1", "", 300, 180, 210);

	TTreeReader r1("B4", f);
	TTreeReaderValue<Double_t> Ecore(r1, "Ecore");
	TTreeReaderValue<Double_t> Eann1(r1, "Eann1");
	TTreeReaderValue<Double_t> Eann2(r1, "Eann2");
	TTreeReaderValue<Double_t> Eann3(r1, "Eann3");
	TTreeReaderValue<Double_t> Eann4(r1, "Eann4");
	TTreeReaderValue<Double_t> Eann5(r1, "Eann5");
	TTreeReaderValue<Double_t> Eann6(r1, "Eann6");

	while (r1.Next())
	{
		h1->Fill(*Ecore + *Eann1 + *Eann2 + *Eann3 + *Eann4 + *Eann5 + *Eann6);
	}

  	h1->GetXaxis()->SetTitle("Energy (MeV)");
  	h1->GetYaxis()->SetTitle("Counts");
  	h1->SetTitle("Energy Recorded by G4 CATS sim -- 200MeV Photon Beam");
	h1->Draw();

  	//Start of FWHM Section
  	//Creating a Skewed Gaussian Function
  	Double_t BinWithMostCounts = h1->GetMaximumBin();
  	Double_t MaxYValue = h1->GetBinContent(BinWithMostCounts);
  	Double_t CenterPeak = h1->GetBinCenter(BinWithMostCounts);
	Double_t StdDev = h1->GetStdDev();

	TF1 *graph_function = new TF1("graph_function", MyMethod, 180, 210, 6);

	//0 is the mean; 1 is the standard deviation; 2 is E; 3 is w; 4 is a (the shape parameter); 5 is the height
	graph_function->SetParameter(0, 0);
	graph_function->SetParameter(1, StdDev);
	graph_function->SetParameter(2, CenterPeak);
	graph_function->SetParameter(3, 1);
	graph_function->SetParameter(4, -8);
	graph_function->SetParameter(5, 1);

	//Scaling
	Double_t max_of_function = graph_function->GetMaximum();
	Double_t scale_factor = MaxYValue/max_of_function;
	graph_function->SetParameter(5, scale_factor);

	h1->Fit("graph_function");
	
	Double_t parameter0 = graph_function->GetParameter(0);
	cout << parameter0 << endl;
	Double_t parameter1 = graph_function->GetParameter(1);
	cout << parameter1 << endl;
	Double_t parameter2 = graph_function->GetParameter(2);
	cout << parameter2 << endl;
	Double_t parameter3 = graph_function->GetParameter(3);
	cout << parameter3 << endl;
	Double_t parameter4 = graph_function->GetParameter(4);
	cout << parameter4 << endl;
	Double_t parameter5 = graph_function->GetParameter(5);
	cout << parameter5 << endl;
}
