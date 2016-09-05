#include <stdio.h>
# include <stdlib.h>
#include <iostream>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include <string>
#include "chi3.cxx"
using namespace std;

//为了处理批量的data文件而构造的函数。返回与chi是一样的,只不过多出了批量add

TH1F* Psi3770pi_data(int n1,int n2,int opt, double cut){


 /*
    int n1=1;
    int n2=120;
    int opt = 4;
    double cut = 0.8;
*/
    int n = n2-n1+1;
    ostringstream s[n];
    int N = 0;
    
    //TFile *f[n];
    //TFile *f = new TFile[n];

    if(opt ==0)TH1F *hh = new TH1F("chi4c","chi4c",100,0.01,60);
    if(opt ==1)TH1F *hh = new TH1F("g1","g1",150,1.4,2.2);
    if(opt ==2)TH1F *hh = new TH1F("g2","g2",150,0.4,2.2);
    if(opt ==3)TH1F *hh = new TH1F("g3","g3",100,0.0,1.2);
    if(opt ==4)TH1F *hh = new TH1F("m2","m2",50,0.10,0.16);
    if(opt ==5)TH1F *hh = new TH1F("mt","mt",150,3.2,4.0);
    if(opt ==6)TH1F *hh = new TH1F("the_decay","the_decay",80,0.0,1.0);
//    if(opt==7)TH1F *hh = new TH1F("p","p",80,1.4,2.0);

    for(int i=n1; i <= n2;i++){
        int p = i-n1;        
        //s[i] << "../p3770Alg/dataroot/0414" << p << ".root"<<endl;

        
        s[p]<<"./psip/psip09_"<<i<<".root";

        //cout << s[p].str() << endl;

        const char* r = s[p].str().c_str();//"test.root";
          
        cout<<r<<endl;
        
        
        TFile *f = new TFile(r);

        /*
        cout<<" "<<endl;
        cout<<"is open "<<f->IsOpen()<<" zombine "<<f->IsZombie()<<endl;
        cout<<" "<<endl;
        */
        
        

        if (f->IsOpen()) {
            
            if ( (TTree*)f->Get("TreeAna") ){
                TH1F *h = chi3(f,opt,cut,1); 
                //h->Draw();
                hh->Add(h);       
            }
        }

        f->Close();

    }


    
//draw --------------------------------------------------------



    //hh->Draw();
    return hh;
   
//    legend->Draw();
    



}
