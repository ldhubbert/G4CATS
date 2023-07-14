const Double_t kPI = TMath::Pi();
Double_t MyMethod( Double_t*, Double_t*);

void PlotStuff()
{
	TCanvas* c1 = new TCanvas( "c1", "Canvas", 20, 20, 1000, 1000);
	c1->Divide( 1, 2);
	c1->cd( 1);

	TString filename = "~/Vincent/G4CATS/Out/B4_200MeV.root";
	TFile *f = TFile::Open(filename);

	TH1F* h1 = new TH1F( "h1", "Histogram Statistics", 300, 180, 210);

	TTreeReader r1( "B4", f);
	TTreeReaderValue<Double_t> Ecore( r1, "Ecore");
	TTreeReaderValue<Double_t> Eann1( r1, "Eann1");
	TTreeReaderValue<Double_t> Eann2( r1, "Eann2");
	TTreeReaderValue<Double_t> Eann3( r1, "Eann3");
	TTreeReaderValue<Double_t> Eann4( r1, "Eann4");
	TTreeReaderValue<Double_t> Eann5( r1, "Eann5");
	TTreeReaderValue<Double_t> Eann6( r1, "Eann6");

	while (r1.Next()) h1->Fill(*Ecore + *Eann1 + *Eann2 + *Eann3 + *Eann4 + *Eann5 + *Eann6);

  	h1->GetXaxis()->SetTitle( "Energy (MeV)");
  	h1->GetYaxis()->SetTitle( "Counts");
  	h1->SetTitle( "Energy Recorded by G4 CATS sim -- 200MeV Photon Beam");
	h1->Draw();

  	Double_t BinWithMostCounts = h1->GetMaximumBin();
  	Double_t MaxYValue = h1->GetBinContent( BinWithMostCounts);
  	Double_t CenterPeak = h1->GetBinCenter( BinWithMostCounts);
	Double_t StdDev = h1->GetStdDev();

//	cout << CenterPeak;
//	cout << "  " << StdDev;
//	cout << endl;

	c1->cd(2);

	TF1 *f1 = new TF1( "f1", MyMethod, 180, 210, 5);

	f1->SetParameter( 0, StdDev);				// standard deviation
	f1->SetParameter( 1, 0);					// average of data
	f1->SetParameter( 2, CenterPeak);		// peak location
	f1->SetParameter( 3, 1);					// scale parameter
	f1->SetParameter( 4, -10);					// shape parameter
 
//	cout << f1->Eval( 180);
//	cout << "  " << f1->Eval( 210);
//	cout << endl;

  	f1->SetTitle( "Skewed Gaussian");
	f1->Draw();

}

Double_t Sqr( Double_t x)
{
	return( x*x);
}

Double_t MyMethod( Double_t *y, Double_t *par)
{
	Double_t exponent, lim, erfterm, function;
	Double_t sig, mu, xi, om, al;
	Double_t x;

	x = y[0];
	sig = par[0];	// sigma -> standard deviation
	mu = par[1];	// mu    -> average of data
	xi = par[2];	// xi    -> location parameter
	om = par[3];	// omega -> scale parameter
	al = par[4];	// alpha -> shape parameter

	exponent = -Sqr( (x-xi)/om - mu)/2/Sqr( sig);
	lim = al*(x-xi)/om/sqrt( 2);
	erfterm = 1 + TMath::Erf( lim);

	function = exp( exponent)*erfterm/2/om/sig/sqrt( 2*kPI);

	return( function);

}

/*
void TestFunc( Double_t x = 100)
{
	Double_t exponent, lim, erfterm, function;
	Double_t sig, mu, xi, om, al;

	sig = 3;		// sigma -> standard deviation
	mu = 0;		// mu    -> average of data
	xi = 200;   // xi    -> location parameter
	om = 1;     // omega -> scale parameter
	al = -10;   // alpha -> shape parameter

	exponent = -Sqr( (x-xi)/om - mu)/2/Sqr( sig);
	lim = al*(x-xi)/om/sqrt( 2);
	erfterm = 1 + TMath::Erf( lim);

	function = exp( exponent)*erfterm/2/om/sig/sqrt( 2*kPI);

	cout << exponent;
	cout << " " << lim;
	cout << " " << term2;
	cout << " " << function;
	cout << endl;

}
*/
