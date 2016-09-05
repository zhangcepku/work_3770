#include <stdio.h>
# include <stdlib.h>
#include <iostream>
//#include "TMVA/Factory.h"
//#include "TMVA/Tools.h"
#include <string>
#include "shape0823.cxx"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

#ifndef __CINT__
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/AxisAngle.h"
#include "Math/Quaternion.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"
#include "Math/BoostX.h"
#include "Math/BoostY.h"
#include "Math/BoostZ.h"
#include "Math/Transform3D.h"
#include "Math/Plane3D.h"
#include "Math/VectorUtil.h"

using namespace ROOT::Math;

#endif

//处理root 文件的底层函数
//TH1F *chi3(TFile *f1,int op,double cut){
TH1F *chi7(TFile *f1,int op,double cut, int cut2){



     double x0,y0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;


    TTree *t1 = (TTree*)f1->Get("TreeAna");

    // ---------------- gamma part -----------------
    Double_t 
    the_decay,p,
    masst;
    Double_t Chisq1, the_decay;

    int runid,evtid;
    int N = 0;

    int Nhit;
    int count;
    Double_t phiE;
    Double_t phiW;
    Double_t phi[10000];


    TLorentzVector *gamma1 = new TLorentzVector();
    TLorentzVector *gamma2 = new TLorentzVector();
    TLorentzVector *gamma3 = new TLorentzVector();


    TLorentzVector gamma,gamma_t1,gamma_t2,gamma_t3;

/*
    Double_t X_1[1000];
    Double_t Y_1[1000];
    Double_t Z_1[1000];
    
    Double_t X_2[1000];
    Double_t Y_2[1000];
    Double_t Z_2[1000];

    Double_t X_4;
    Double_t Y_4;
    Double_t Z_4;

    Double_t X_5;
    Double_t Y_5;
    Double_t Z_5;
*/
    Int_t Nhit = 0;
//    Int_t Nhits[3];
    Int_t count = 0;
    Int_t count_cone = 0;
    Int_t count_phi = 0;

    //Nhit 是统计所有的hit数，8以外的。Nhits是bhabha中根据动量统计出的trk的hits数

/*
    Double_t X0[3];
    Double_t Y0[3];
    Double_t Z0[3];

    Double_t PX0[3];
    Double_t PY0[3];
    Double_t PZ0[3];
*/

/*
    t1->SetBranchAddress("X_1",&X_1);
    t1->SetBranchAddress("Y_1",&Y_1);
    t1->SetBranchAddress("Z_1",&Z_1);
    t1->SetBranchAddress("X_2",&X_2);
    t1->SetBranchAddress("Y_2",&Y_2);
    t1->SetBranchAddress("Z_2",&Z_2);

    t1->SetBranchAddress("X_3",&X_3);
    t1->SetBranchAddress("Y_3",&Y_3);
    t1->SetBranchAddress("Z_3",&Z_3);
    t1->SetBranchAddress("X_4",&X_4);
    t1->SetBranchAddress("Y_4",&Y_4);
    t1->SetBranchAddress("Z_4",&Z_4);
*/    
    t1->SetBranchAddress("count",&count);
//
    t1->SetBranchAddress("Nhit",&Nhit);

/*
    t1->SetBranchAddress("X0",&X0);
    t1->SetBranchAddress("PX0",&PX0);
    t1->SetBranchAddress("Y0",&Y0);
    t1->SetBranchAddress("PY0",&PY0);
    t1->SetBranchAddress("Z0",&Z0);
    t1->SetBranchAddress("PZ0",&PZ0);
    t1->SetBranchAddress("Nhits",&Nhits);
*/


    //t1->SetBranchAddress("mass_2",&mass2);
    t1->SetBranchAddress("mass_t",&masst);
    t1->SetBranchAddress("count_cone",&count_cone);
    //-----------------------------------.rec only
    
    t1->SetBranchAddress("phiE",&phiE);
    t1->SetBranchAddress("phiW",&phiW);
    t1->SetBranchAddress("phi",&phi);
    
    
    
    t1->SetBranchAddress("gamma1_unfitted",&gamma1);
    t1->SetBranchAddress("gamma2_unfitted",&gamma2);
    //t1->SetBranchAddress("gamma2",&gamma2);
    t1->SetBranchAddress("gamma3_unfitted",&gamma3);

    t1->SetBranchAddress("runid",&runid);
    t1->SetBranchAddress("evtid",&evtid);

    //t1->SetBranchAddress("gamma",&gamma);
    

    t1->SetBranchAddress("Chisq_low",&Chisq1);
    //t1->SetBranchAddress("m_thedecay",&the_decay);

    TH1F *h2 = new TH1F("Invariant_Mass_2","Invariant_Mass_2",50,0.10,0.16);

    TH1F *h3 = new TH1F("Invariant_Mass_t","Invariant_Mass_t",150,3.2,4.0);

    TH1F *h4 = new TH1F("the_decay","the_decay",80,0.0,1.0);//6

    TH1F *h5 = new TH1F("p","p",80,1.4,2.0);
    TH1F *h6 = new TH1F("count","count",70,0.0,70);
    
    
    TH1F *g1 = new TH1F("gamma1","gamma1_e",150,1.4,2.2);
    TH1F *g2 = new TH1F("gamma2","gamma2_e",150,0.4,2.2);
    TH1F *g3 = new TH1F("gamma3","gamma3_e",100,0,1.2);



    TH1F *m1 = new TH1F("Chisq_4C","Chisq_4C",100,0.01,60);

    //-----------------------------------------------
    
    
    Long64_t nentries1 = t1->GetEntries();

//cout<<"done initial"<<endl;
    
    for (Long64_t i=0;i<nentries1;i++) {
        
        t1->GetEntry(i);

/*
        x4 = PX0[0];
        y4 = PY0[0];
        z4 = PZ0[0];
        
        x5 = PX0[1];
        y5 = PY0[1];
        z5 = PZ0[1];
*/
/*        
        x4 = X_4;
        y4 = X_4;
        z4 = X_4;
        
        x5 = X_5;
        y5 = X_5;
        z5 = X_5;
*/
//        XYZVector *g4 = new XYZVector(x4,y4,z4);
//        XYZVector *g5 = new XYZVector(x5,y5,z5);
/*
        double dgree = 0.12;
        double l1 = sqrt(g4->mag2());
        double l2 = sqrt(g5->mag2());
        double r1 = l1*TMath::Tan(dgree);
        double r2 = l2*TMath::Tan(dgree);
*/

        double r = 6.0;
        double  l = 100;


/*
        x4 = l*x4/sqrt(g4->mag2());
        y4 = l*y4/sqrt(g4->mag2());
        z4 = l*z4/sqrt(g4->mag2());
        x5 = l*x5/sqrt(g5->mag2());
        y5 = l*y5/sqrt(g5->mag2());        
        z5 = l*z5/sqrt(g5->mag2());
*/
        double L = 230;
        double  LL = L;

        //make cone
int m = 0;
for(int j =0;j<count;j++){if(phiE<=phi[j]<=phiW)m++;}


        gamma_t3 = *gamma1+*gamma2;
        gamma_t2 = *gamma1+*gamma3;
        gamma_t1 = *gamma2+*gamma3;

        double masspi = 0.1349766;
        int flag =0;
        double e1,e2,e3;
        double mass2 = 0;

        e1 = gamma_t1.M() - masspi;if(e1<0)e1=-e1;
        e2 = gamma_t2.M() - masspi;if(e2<0)e2=-e2;
        e3 = gamma_t3.M() - masspi;if(e3<0)e3=-e3;

        if(e1<e2 && e1<e3){
            flag = 1;
            //cout<<"1"<<endl;
            gamma = *gamma2 + *gamma3;
            the_decay = (gamma2->Energy() - gamma3->Energy())/gamma.P();
            mass2 = gamma.M();
            //cout<<mass2<<endl;
        
        }
        //else{cout<<e1<<" "<<e2<<" "<<e3<<" "<<endl;}
        
        if(e2<e1 && e2<e3){
            flag =2;
            //cout<<"2"<<endl;
            gamma = gamma_t2;
            the_decay = (gamma1->Energy() - gamma3->Energy())/gamma.P();
            //if(the_decay<0.8)cout<<"------------------------------------------"<<endl;
            mass2 = gamma.M();
        }

        if(e3<e1 && e3<e2){
            flag =3;
           // cout<<"3"<<endl;
            gamma = gamma_t3;
            the_decay = (gamma2->Energy() - gamma1->Energy())/gamma.P();
            mass2 = gamma.M();
        }

        if (Chisq1<=50 && mass2>=0.1 && mass2 <=0.16){

            //gamma = *gamma2 + *gamma3;
            //the_decay = (gamma2->Energy() - gamma3->Energy())/gamma.P();
            //p = gamma.P();
            if(the_decay<0)the_decay = - the_decay;

            //if(the_decay<1.0){
            if (the_decay<cut && m<=cut2){
            //if (the_decay<cut){

            //if (the_decay<0.8&& p>1.8){




                h4->Fill(the_decay);
                h6->Fill(m);

                //h5->Fill(p);

                m1->Fill(Chisq1);
                h2->Fill(mass2);
                h3->Fill(masst);
            //h4->Fill(the_decay);
                if(flag ==1){
                    g1->Fill(gamma1->Energy());
                    g2->Fill(gamma2->Energy());
                    g3->Fill(gamma3->Energy());
                }
                if(flag ==2){
                    g1->Fill(gamma2->Energy());
                    g2->Fill(gamma1->Energy());
                    g3->Fill(gamma3->Energy());
                }
                if(flag ==3){
                    g1->Fill(gamma3->Energy());
                    g2->Fill(gamma1->Energy());
                    g3->Fill(gamma2->Energy());
                }
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
    //if(op==6) return h4;
    if(op==7) return h6;

}