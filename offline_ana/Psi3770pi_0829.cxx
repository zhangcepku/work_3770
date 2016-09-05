#include "Psi3770pi_0426.cxx"

void Psi3770pi_0829(){
	double cut2 = 30;
	TH1F* hd6 = Psi3770pi_0426(6);
	TH1F* hd4 = Psi3770pi_0426(4);
	TH1F* hd2 = Psi3770pi_0426(2);
//	TH1F* hd1 = Psi3770pi_0426(2);

	TCanvas *c = new TCanvas(); c->cd();
	hd6->SetFillColor(kGreen-10);hd6->SetAxisRange(0,140,"Y");
	hd6->Draw();
	hd4->SetFillColor(kGreen);hd4->Draw("same");
	hd2->SetFillColor(kGreen+3);hd2->Draw("same");
	//hd1->SetLineColor(kBlue-2);hd1->Draw("same");

	TLegend *legend = new TLegend(0.6,0.65,0.7,0.8);

    legend->AddEntry(hd6,"cut=6");
    
    legend->AddEntry(hd4,"cut=4");

    legend->AddEntry(hd2,"cut=2");

    //legend->AddEntry(hd5,"data-cut=6");
	legend->Draw("same");


}