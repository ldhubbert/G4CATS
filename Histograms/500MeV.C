void maincode()
{
	TCanvas *c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
	c1->Divide(1,1);
	c1->cd(1);

	TString filename = "~/Vincent/G4CATS/Out/B4_500MeV.root";
	TFile *f = TFile::Open(filename);

	TH1F *h1 = new TH1F("Histogram Statistics", "", 600, 455, 515);

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
	h1->GetXaxis()->CenterTitle();
	h1->GetYaxis()->CenterTitle();
	h1->SetTitle("500MeV Incident Beam");
	h1->Draw();
}

