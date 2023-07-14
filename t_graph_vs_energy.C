#include <TMath.h>
#include <TF1.h>
#include <TAxis.h>

Double_t t_graph(Double_t *x)
{
	//x[0] represents energy
	
	//Two parameters
	Double_t mass_proton = 938.27208816;
	//Changing 30 degrees into radians
	Double_t theta = 30 * (TMath::Pi() / 180);

	Double_t t = (-2 * x[0] * x[0] * mass_proton * (1 - TMath::Cos(theta))) / (mass_proton + (x[0] * (1 - TMath::Cos(theta))));
	return t;
}

void maincode()
{
	TCanvas *c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);

	int array_size = 15;
	Double_t x_values[array_size];
	Double_t y_values[array_size];

	for (int i = 100; i <= 1500; i = i + 100)
	{
		Double_t double_i = static_cast<Double_t>(i);
		int array_numbering = (i / 100) - 1;
		x_values[array_numbering] = double_i;
		y_values[array_numbering] = t_graph(&double_i)/1000000;
	}
	
	cout << y_values[0] << endl;

	TGraph *g1 = new TGraph(array_size, x_values, y_values);
	g1->SetTitle("t vs. Incident Energy for a 30-Degree Scattering Angle");
	g1->GetXaxis()->SetTitle("Energy (MeV)");
	g1->GetXaxis()->CenterTitle();
	g1->GetYaxis()->SetTitle("t (GeV^2)");
	g1->GetYaxis()->CenterTitle();
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(21);
	g1->Draw("AP");
}
