#这是bak 在跑3770 rec.file之前就拷贝过来的chi3,应该是做？用的
#include <stdio.h>
# include <stdlib.h>
#include <iostream>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include <string>


//处理root 文件的底层函数
//TH1F *chi3(TFile *f1,int op,double cut){
TH1F *chi3(TFile *f1,int op,double cut, int cut2){

	TTree *t1 = (TTree*)f1->Get("TreeAna");
	


	// ---------------- gamma part -----------------
	Double_t 
    the_decay,p,
    masst;
    Double_t Chisq1, the_decay;
    double mass2;

    int runid,evtid;
    int N = 0;


/*
    TLorentzVector *gamma1 = new TLorentzVector();
    TLorentzVector *gamma2 = new TLorentzVector();
    TLorentzVector *gamma3 = new TLorentzVector();
*/




    t1->SetBranchAddress("mass_2",&mass2);

    
    
 /*   
    t1->SetBranchAddress("gamma1_unfitted",&gamma1);
    t1->SetBranchAddress("gamma2_unfitted",&gamma2);
    //t1->SetBranchAddress("gamma2",&gamma2);
    t1->SetBranchAddress("gamma3_unfitted",&gamma3);

    //t1->SetBranchAddress("gamma",&gamma);
    
*/
    t1->SetBranchAddress("Chisq_low",&Chisq1);
    t1->SetBranchAddress("m_thedecay",&the_decay);

    TH1F *h2 = new TH1F("Invariant_Mass_2","Invariant_Mass_2",100,0.10,0.14);




    //-----------------------------------------------
	
	
	Long64_t nentries1 = t1->GetEntries();

    
    for (Long64_t i=0;i<nentries1;i++) {
        
        t1->GetEntry(i);



       



 

        if (Chisq1<=50 && mass2>=0.10 && mass2 <=0.14){
        //if(Chisq1<=50){

            if(the_decay<=0.8){
            h2->Fill(mass2);
            }
            
            //对于data 暂时没有the_decay,故自己写，否则要解除上面的注释h4

    
        }    
    }

	if(op==0) return m1;

	if(op==1) return g1;
    if(op==2) return g2;
    if(op==3) return g3;
    if(op==4) return h2;
    if(op==5) return h3;
    if(op==6) return h4;
    if(op==7) return h6;

}