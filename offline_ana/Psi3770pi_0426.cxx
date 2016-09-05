
#include <stdio.h>
# include <stdlib.h>
#include <iostream>
#include <string>
//#include "Psi3770pi_data.cxx"

#include <string>
//#include <Math/Point3D.h>



#include "chi2.cxx"
#include "chi3.cxx"
#include "chi7.cxx"
#include "chi6.cxx"
using namespace std;


TH1F* Psi3770pi_0426(double cut2){

	Color_t lcolor1 = kRed;
	Color_t lcolor2 = kGreen;
	Color_t lcolor3 = kBlue;
	Color_t lcolor4 = kBlack;
	Color_t lcolor5 = kYellow;
	Color_t lcolor6 = kOrange;
    Color_t lcolor7 = kMagenta+2;
    Color_t lcolor8 = kTeal-6;
    Color_t lcolor9 = kCyan;
    Color_t lcolor10 = kAzure+7;
    Color_t lcolor11 = kSpring-7;
    Color_t lcolor12 = kBlue+3;

	Color_t lcolor13 = kViolet-7;

	Color_t lcolor14 = kViolet-6;
	Color_t lcolor15 = kViolet-5;
	Color_t lcolor16 = kViolet-4;
	Color_t lcolor17 = kViolet-9;

  Color_t lcolor18 = kGreen+1;
  Color_t lcolor19 = kGreen+2;
  Color_t lcolor20 = kGreen+3;
  Color_t lcolor21 = kGreen+4;

  
    //TCanvas *c = new TCanvas("mass2","mass2");
  	//TH1F *hh = Psi3770pi_data(1,427,4);
  	/*
    TH1F* hh = new TH1F();
    hh->SetLineWidth(2);hh->SetLineColor(lcolor6);

    hh->SetTitle(
  		"the_decay");
      //"mass-pi0");
  		//"mass-3770");
  		//"gamma1");
		//"gamma2");
		//"gamma3");
		//"chisq-4c");
  	*/
/*
    TH1F *hqq1 = Psi3770pi_newmc(1,0);//0-qq
    TH1F *hqq2 = Psi3770pi_newmc(2,0);
    TH1F *hqq3 = Psi3770pi_newmc(3,0);
*/
    //TH1F *hee = Psi3770pi_newmc(4,1);

    //TH1F *hd1 = Psi3770pi_data(1,120,7,0.8);


/*
    TH1F *hd3 = Psi3770pi_data(1,427,4,0.8);
    TH1F *hd4 = Psi3770pi_data(1,427,4,0.7);
    TH1F *hd5 = Psi3770pi_data(1,427,4,0.6);
  */  
    //TH1F *ht2 = Psi3770pi_mc(5,2);

    //hqq->SetLineColor(lcolor3);
    //hee->SetLineColor(lcolor2);
/*
    TFile *f1 = new TFile("../root/sigmc/sigmc042501.root");
    TH1F *h1 = chi3(f1,1,0.8);
    
    TFile *f2 = new TFile("../root/sigmc/sigmc042502.root");
    TH1F *h2 = chi3(f2,1,0.8);
    
    TFile *f3 = new TFile("../root/sigmc/sigmc042501.root");
    TH1F *h3 = chi3(f3,1,0.8);
*/
    /*
    int o = 4;int k = 6;
    
    TFile *fd1 = new TFile("./rootnew/data1_0521.root");
    TH1F *hd1 = chi3(fd1,o,0.8,k);
    
    TFile *fd2 = new TFile("./rootnew/data2_0521.root");
    TH1F *hd2 = chi3(fd2,o,0.8,k);
/*
    TFile *fd3 = new TFile("./rootnew/data1_0521.root");
    TH1F *hd3 = chi3(fd3,o,0.8,10);
    
    TFile *fd4 = new TFile("./rootnew/data2_0521.root");
    TH1F *hd4 = chi3(fd4,o,0.8,10);

    TFile *fd5 = new TFile("./rootnew/data1_0521.root");
    TH1F *hd5 = chi3(fd5,o,0.8,6);
    
    TFile *fd6 = new TFile("./rootnew/data2_0521.root");
    TH1F *hd6 = chi3(fd6,o,0.8,6);
*/

  //  hd1->Add(hd2);
    //hd3->Add(hd4);
    //hd5->Add(hd6);
    /*
    TFile *fs = new TFile("./rootnew/sigmc0521.root");
    TH1F *hs = chi3(fs,o,0.8,k);

    
/*
    TFile *f1 = new TFile("./rootnew/ee1.root");
    TH1F *h1 = chi3(f1,o,0.8,k);
    
    TFile *f2 = new TFile("./rootnew/ee2.root");
    TH1F *h2 = chi3(f2,o,0.8,k);
    
    TFile *f3 = new TFile("./rootnew/ee3.root");
    TH1F *h3 = chi3(f3,o,0.8,k);
    TFile *f4 = new TFile("./rootnew/ee4.root");
    TH1F *h4 = chi3(f4,o,0.8,k);
    TFile *f5 = new TFile("./rootnew/ee5.root");
    TH1F *h5 = chi3(f5,o,0.8,k);
    
    h1->Add(h2);h1->Add(h3);h1->Add(h4);h1->Add(h5);
*/

    /*
    TFile *f4 = new TFile("../rootnew/");
    TH1F *h4 = chi3(f4,7,0.8);
*/
  /*  
    TFile *fee = new TFile("../root/eemc/mc_ee_0504.root");

    TH1F *hee1 = chi3(fee,4,1.0);hee1->SetLineColor(lcolor2);
    
    TH1F *hee2 = chi3(fee,4,0.9);hee2->SetLineColor(lcolor18);
    TH1F *hee3 = chi3(fee,4,0.8);hee3->SetLineColor(lcolor19);
    TH1F *hee4 = chi3(fee,4,0.7);hee4->SetLineColor(lcolor20);
    TH1F *hee5 = chi3(fee,4,0.6);hee5->SetLineColor(lcolor21);
    */  

    //TH1F *hh = new TH1F("p","p",80,1.4,2.0);

    //TH1F *hh = new TH1F("Invariant_Mass_2","Invariant_Mass_2",100,0.08,0.18);
    //TH1F *hh = new TH1F("the_decay","the_decay",80,0.0,1.0);//6
    //h1->Add(h2);h1->Add(h3);//hh->Add(h3);

    //hh->SetLineColor(lcolor6);

  	//hqq->SetNormFactor(1133);
    //hee->SetNormFactor(316);
    //ht1->Add(ht2,50);
  	//h2->SetLineColor(lcolor6);
    //hee->DrawNormalized();
    //hqq->Add(hee);
    //hqq->Add(hee,97.6);
    //hqq->SetNormFactor(1449);
    //hd1->SetLineColor(lcolor1);

int f = 3;
int op = 4;


//chi3,chi6专门为相关root版本设计。chi2能统计wire(针对0825记录了X1245)

// 0825 two cone 0.12
if(f==0){
    double cut1 = 0.8;
    double cut2 = 2;

        TFile *fd1 = new TFile("./root/3770data_rec7_0825.root");
    TH1F *hd1 = chi3(fd1,op,cut1,cut2);
        TFile *fd2 = new TFile("./root/3770data_rec1_0825.root");
    TH1F *hd2 = chi3(fd2,op,cut1,cut2);
        TFile *fd3 = new TFile("./root/3770data_rec2_0825.root");
    TH1F *hd3 = chi3(fd3,op,cut1,cut2);
        TFile *fd4 = new TFile("./root/3770data_rec3_0825.root");
    TH1F *hd4 = chi3(fd4,op,cut1,cut2);
        TFile *fd5 = new TFile("./root/3770data_rec4_0825.root");
    TH1F *hd5 = chi3(fd5,op,cut1,cut2);
        TFile *fd6 = new TFile("./root/3770data_rec5_0825.root");
    TH1F *hd6 = chi3(fd6,op,cut1,cut2);
        TFile *fd7 = new TFile("./root/3770data_rec6_0825.root");
    TH1F *hd7 = chi3(fd7,op,cut1,cut2);
    /*
        TFile *fd1 = new TFile("./root/3770data_rec7_0825.root");
    TH1F *hd1 = chi2(fd1,op,cut1,cut2);
    
    TFile *fd2 = new TFile("./root/3770data_rec1_0825.root");
    TH1F *hd2 = chi2(fd2,op,cut1,cut2);

        TFile *fd3 = new TFile("./root/3770data_rec2_0825.root");
    TH1F *hd3 = chi2(fd3,op,cut1,cut2);
        TFile *fd4 = new TFile("./root/3770data_rec3_0825.root");
    TH1F *hd4 = chi2(fd4,op,cut1,cut2);
        TFile *fd5 = new TFile("./root/3770data_rec4_0825.root");
    TH1F *hd5 = chi2(fd5,op,cut1,cut2);
        TFile *fd6 = new TFile("./root/3770data_rec5_0825.root");
    TH1F *hd6 = chi2(fd6,op,cut1,cut2);
        TFile *fd7 = new TFile("./root/3770data_rec6_0825.root");
    TH1F *hd7 = chi2(fd7,op,cut1,cut2);
    */
}

// 

if(f==1){
    double cut1 = 0.8;
    

    TFile *fd1 = new TFile("./root/3770data_rec1_0829.root");
    TH1F *hd1 = chi6(fd1,op,cut1,cut2);
    
    TFile *fd2 = new TFile("./root/3770data_rec2_0829.root");
    TH1F *hd2 = chi6(fd2,op,cut1,cut2);

    TFile *fd3 = new TFile("./root/3770data_rec3_0829.root");
    TH1F *hd3 = chi6(fd3,op,cut1,cut2);

    TFile *fd4 = new TFile("./root/3770data_rec4_0829.root");
    TH1F *hd4 = chi6(fd4,op,cut1,cut2);
    
    TFile *fd5 = new TFile("./root/3770data_rec5_0829.root");
    TH1F *hd5 = chi6(fd5,op,cut1,cut2);

    TFile *fd6 = new TFile("./root/3770data_rec6_0829.root");
    TH1F *hd6 = chi6(fd6,op,cut1,cut2);


}

if(f==2){

    double cut1 = 0.8;


    TFile *fd1 = new TFile("./root/3770data_rec1_0822.root");
    TH1F *hd1 = chi3(fd1,op,cut1,cut2);
    
    TFile *fd2 = new TFile("./root/3770data_rec2_0822.root");
    TH1F *hd2 = chi3(fd2,op,cut1,cut2);
}

if(f==3){

    double cut1 = 0.8;
    //double cut2 = 2;


    TFile *fd1 = new TFile("./root/3770data_rec1_0830.root");
    TH1F *hd1 = chi6(fd1,op,cut1,cut2);
    
    TFile *fd2 = new TFile("./root/3770data_rec2_0830.root");
    TH1F *hd2 = chi6(fd2,op,cut1,cut2);

        TFile *fd3 = new TFile("./root/3770data_rec3_0830.root");
    TH1F *hd3 = chi6(fd3,op,cut1,cut2);

        TFile *fd5 = new TFile("./root/3770data_rec5_0830.root");
    TH1F *hd5 = chi6(fd5,op,cut1,cut2);
            TFile *fd4 = new TFile("./root/3770data_rec4_0830.root");
    TH1F *hd4 = chi6(fd4,op,cut1,cut2);
            TFile *fd6 = new TFile("./root/3770data_rec6_0830.root");
    TH1F *hd6 = chi6(fd6,op,cut1,cut2);
}


if(f==4){

    double cut1 = 0.8;
    double cut2 = 2;


    TFile *fd1 = new TFile("./root/3770data_rec1_0826.root");
    TH1F *hd1 = chi6(fd1,op,cut1,cut2);
    
    TFile *fd2 = new TFile("./root/3770data_rec2_0826.root");
    TH1F *hd2 = chi6(fd2,op,cut1,cut2);

    TFile *fd3 = new TFile("./root/3770data_rec3_0826.root");
    TH1F *hd3 = chi6(fd3,op,cut1,cut2);
}

if(f==5){

    double cut1 = 0.8;
    //double cut2 = 2;


    TFile *fd1 = new TFile("./root/sigmc1_0818.root");
    TH1F *hd1 = chi3(fd1,op,cut1,cut2);
  
}
    
    TFile *fd0 = new TFile("./root/qqbar6_0826.root");
    TH1F *hd0 = chi3(fd0,op,cut1,cut2);




    //TLegend *legend = new TLegend(0.6,0.65,0.7,0.8);

    //legend->AddEntry(hs,"signal");
    /*
    legend->AddEntry(hd1,"data-cut=20");

    legend->AddEntry(hd3,"data-cut=10");

    legend->AddEntry(hd5,"data-cut=6");
    //legend->AddEntry(h1,"ee");
    */
    TCanvas *c1 = new TCanvas();c1->cd();
    if(f!=5)hd1->Add(hd2);
    if(f==0){hd1->Add(hd3);hd1->Add(hd4);hd1->Add(hd5);hd1->Add(hd6);hd1->Add(hd7);}
    if(f==4){hd1->Add(hd3);}
    if(f==1){hd1->Add(hd3);hd1->Add(hd4);hd1->Add(hd5);hd1->Add(hd6);}
    if(f==3){hd1->Add(hd3);hd1->Add(hd5);hd1->Add(hd4);hd1->Add(hd6);}
    hd1->SetTitle("mdc hit count");

    //hs->SetLineColor(kRed);hd1->SetLineColor(kBlack);
    //hd2->SetLineColor(kRed);
    //h3->SetLineColor(kRed);
    hd1->SetLineColor(kBlack);//hd3->SetLineColor(lcolor14);hd5->SetLineColor(lcolor15); //hd5->SetFillColor(lcolor3);
    //hd1->SetHistLineStyle();
    //hd1->SaveAs("result.root");
    hd0->SetLineColor(kBlue);
    hd0->SetNormFactor(10*hd0->GetEntries());
    //hd1->Draw("E");
    //hd0->Draw("same");
    return hd1;
    
    //hd3->Draw("same");
    //hd5->Draw("same");
    //legend->Draw();
/*
    hs->DrawNormalized();
    hd1->DrawNormalized("same");
    h1->DrawNormalized("same");
*/
   // h1->SetLineColor(kBlack);
    //h1->Draw("same");
    //hh->SetTitle("signal- the_decay");
    //hd->DrawNormalized();
    //hd->SetTitle("p_pi0");
    
    //h1->SetLineColor(lcolor1);
    
    //h1->DrawNormalized("P");
    //hd1->DrawNormalized("same");
    //h1->SetLineColor(lcolor1);
    //hee1->Draw("same");
    //h1->SetTitle("signal");
    //h1->DrawNormalized();
    //hd1->DrawNormalized("same");
    //hd1->DrawNormalized("same");
 /*
    TCanvas *c2 = new TCanvas();c2->cd();
    hs->DrawNormalized();
    hd1->DrawNormalized("same");
    //hqq2->Draw();
    //hee2->SetLineColor(lcolor14);
    //hee2->Draw("same");
   /*
    TCanvas *c3 = new TCanvas();c3->cd();
    //hqq3->Draw();
    //hee3->SetLineColor(lcolor15);
    //hee3->Draw("same");
    //TCanvas *c4 = new TCanvas();c4->cd();
    hd3->Draw();
    TCanvas *c4 = new TCanvas();c4->cd();
    //hqq3->Draw();
    //hee3->SetLineColor(lcolor15);
    //hee3->Draw("same");
    //TCanvas *c4 = new TCanvas();c4->cd();
    hd5->Draw();
    //hee4->SetLineColor(lcolor16);
    //hee4->Draw("same");
    //TCanvas *c5 = new TCanvas();c5->cd();
    
    //hee5->SetLineColor(lcolor17);
    //hee5->Draw("same");


/*
    legend->AddEntry(hee4,"cut=0.7");
    legend->AddEntry(hee5,"cut=0.6");
*/

  //legend->AddEntry(h1,"sig");
  //legend->AddEntry(hd1,"data");legend->AddEntry(hee1,"ee");  
    
    
    
    //hee5->Draw("same");

    //legend->Draw();
    

    /*
    hee->SetTitle("p_pi0");
    hee->SetLineColor(lcolor1);
    hee->DrawNormalized();
    */
    //hd->DrawNormalized("same");
//    hqq->DrawNormalized("same");
    //TCanvas *c2 = new TCanvas();c2->cd();
    //hqq->Draw();

    //hee->DrawNormalized();
    //hqq->DrawNormalized("same");
    
    //hqq->DrawNormalized("same");
  	
  	
  }