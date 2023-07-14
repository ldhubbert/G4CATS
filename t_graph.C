#include <TMath.h>
#include <TF1.h>
#include <TAxis.h>

//This code has the rearranged equation for t. Nevertheless, this code should produce the same information as that of the t_graph_og_equation.C file
Double_t t_graph(Double_t *x)
{
	//x[0] represents theta
	
	//Two parameters
	Double_t initial_energy = 1000;
	Double_t mass_proton = 938.27208816;

	Double_t t = (-2 * initial_energy * initial_energy * mass_proton * (1 - TMath::Cos(x[0]))) / (mass_proton + (initial_energy * (1 - TMath::Cos(x[0]))));
	return t;
}

void maincode()
{
	TCanvas *c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);

	int array_size = 8;
	//Here, the x_values correspond to theta
	Double_t x_values[array_size];
	//Here, the y_values correspond to the t values
	Double_t y_values[array_size];

	for (int i = 20; i <= 160; i = i + 20)
	{
		//This is to pass a double into t
		//Also, change from degrees to radians for the TMath::Cos statement to work
		Double_t double_i = static_cast<Double_t>(i) * (TMath::Pi() / 180);
		//So array starts numbering at 0
		int array_numbering = (i/20) - 1;
		x_values[array_numbering] = i;
		y_values[array_numbering] = t_graph(&double_i)/1000000;

		cout << x_values[array_numbering] << endl;
	}

	cout << "///" << endl;
	cout << x_values[0] << endl;
	cout << x_values[1] << endl;

	TGraph *g1 = new TGraph(array_size, x_values, y_values);
	g1->SetTitle("t vs. Theta for a 1000MeV Incident Photon Beam");
	g1->GetYaxis()->SetTitle("t (in GeV^2)");
	g1->GetXaxis()->SetTitle("Theta (in degrees)");
	g1->GetYaxis()->CenterTitle();
	g1->GetXaxis()->CenterTitle();
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(21);
	g1->Draw("AP");
}
