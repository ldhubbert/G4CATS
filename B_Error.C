//Define the function that we want the computer to fit to the unsmeared data
Double_t fit_function(Double_t *x, Double_t *par)
{
	Double_t function = TMath::Exp(-1 * par[0] * x[0]);
	return function;
}

void maincode()
{
	TCanvas *c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
	c1->Divide(2, 1);
	c1->GetFrame()->SetBorderSize(1);


	//Creating Gaussian function
	//Standard deviation will be defined for each y-value in the "for" loop
	TF1 *gaussian = new TF1("gaussian", "gaus", -10, 10);
	gaussian->SetParameter(0, 1);
	gaussian->SetParameter(1, 1);

	//Defining arrays to be filled with data from the "for" loop
	int array_size = 11;
	Double_t x_values[array_size];
	Double_t y_values[array_size];
	Double_t y_value_errors[array_size];
	Double_t y_smeared_errors[array_size];
	Double_t gaussian_sigma[array_size];
	Double_t smeared_data_points[array_size];

	//In order to place results from the "for" loop into the arrays, x has to be an integer
	//Since we actually want to go from x = 0 to x = 1 in 0.1 steps, divide x by 10 in the "for" loop
	for (int x = 0; x <= 10; x++)
	{
		//Since x is an integer, it has to be momentarily changed to a double for this calculation to work
		Double_t realx = static_cast<Double_t>(x)/10;
		x_values[x] = realx;
		//y_values array are filled with the y-values of this exponential decay function
		Double_t function = TMath::Exp(-1*realx);
		y_values[x] = function;
		y_value_errors[x] = 0.05*y_values[x];

		//This gives a different gaussian sigma for each y-value
		gaussian_sigma[x] = 0.05*y_values[x];
		//Parameter 2 is redefined each time through the loop, which we want
		gaussian->SetParameter(2, gaussian_sigma[x]);

		smeared_data_points[x] = y_values[x] * gaussian->GetRandom();
		y_smeared_errors[x] = 0.05*smeared_data_points[x];

	}
	
	//We want error bars in the x-direction to be zero, so fill the x_value_errors array with zeroes
	Double_t x_value_errors[array_size];
	for (int i = 0; i <= 10; i++)
	{
		x_value_errors[i] = 0;
	}

	//Smeared graph
	c1->cd(2);
	TGraphErrors *smeared = new TGraphErrors(11, x_values, smeared_data_points, x_value_errors, y_smeared_errors);
	smeared->SetTitle("Smeared Data");
	//control y-axis range
	smeared->GetYaxis()->SetRangeUser(0.2, 1.3);
	smeared->GetXaxis()->SetTitle("t");
	smeared->GetYaxis()->SetTitle("Differential Cross Section");
	//Title offset prevents axis titles from getting cut-off
	smeared->GetXaxis()->SetTitleOffset(1.0);
	smeared->GetYaxis()->SetTitleOffset(1.4);
	smeared->GetXaxis()->CenterTitle();
	smeared->GetYaxis()->CenterTitle();
	smeared->SetMarkerColor(4);
	smeared->SetMarkerStyle(21);
	smeared->Draw("AP");

	//Unsmeared graph
	c1->cd(1);
	TGraphErrors *unsmeared = new TGraphErrors(11, x_values, y_values, x_value_errors, y_value_errors);
	unsmeared->SetTitle("Unsmeared Data");
	//control y-axis range (make it the same as for the smeared graph)
	unsmeared->GetYaxis()->SetRangeUser(0.2, 1.3);
	unsmeared->GetXaxis()->SetTitle("t");
	unsmeared->GetYaxis()->SetTitle("Differential Cross Section");
	//Title offset prevents axis titles from getting cut-off
	unsmeared->GetXaxis()->SetTitleOffset(1.0);
	unsmeared->GetYaxis()->SetTitleOffset(1.4);
	unsmeared->GetXaxis()->CenterTitle();
	unsmeared->GetYaxis()->CenterTitle();
	unsmeared->SetMarkerColor(4);
	unsmeared->SetMarkerStyle(21);
	unsmeared->Draw("AP");

	//Fit exponential decay curve to smeared data
	//Fifth parameter means there is one parameter (par[0] in the Double_t method at the top of the code)
	//The computer will give us the best par[0] value for this fit
	TF1 *f1 = new TF1("f1", fit_function, 0, 1, 1);
	smeared->Fit("f1");

	//Print parameter on screen:
	gStyle->SetOptFit(111);
}
