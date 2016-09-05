#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
//#include <Math/Point3D.h>


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

void shape0822(){

	 gSystem->Load("libGenVector");
     #ifdef __CINT__
     gSystem->Load("libMathCore");
     using namespace ROOT::Math;
     #endif


     double x0,y0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x3_,y3_,z3_,x4,y4,z4,x5,y5,z5;

	TFile *fd1 = new TFile("./root/3770data_rec7_0825.root");
	TTree *t1 = (TTree*)fd1->Get("TreeAna");

	Int_t runid =0;
	Int_t evtid = 0;
    Double_t X_1[1000];
    Double_t Y_1[1000];
    Double_t Z_1[1000];
    
    Double_t X_2[1000];
    Double_t Y_2[1000];
    Double_t Z_2[1000];

    Int_t Nhit = 0;
    Int_t count = 0;

    Double_t X_4;
    Double_t Y_4;
    Double_t Z_4;
    
    Double_t X_5;
    Double_t Y_5;
    Double_t Z_5;

    t1->SetBranchAddress("X_1",&X_1);
    t1->SetBranchAddress("Y_1",&Y_1);
    t1->SetBranchAddress("Z_1",&Z_1);
    t1->SetBranchAddress("X_2",&X_2);
    t1->SetBranchAddress("Y_2",&Y_2);
    t1->SetBranchAddress("Z_2",&Z_2);

    t1->SetBranchAddress("X_4",&X_4);
    t1->SetBranchAddress("Y_4",&Y_4);
    t1->SetBranchAddress("Z_4",&Z_4);
    t1->SetBranchAddress("X_5",&X_5);
    t1->SetBranchAddress("Y_5",&Y_5);
    t1->SetBranchAddress("Z_5",&Z_5);
    t1->SetBranchAddress("count",&count);
    t1->SetBranchAddress("Nhit",&Nhit);

	t1->SetBranchAddress("runid",&runid);
	t1->SetBranchAddress("evtid",&evtid);

    //test zone

	Long64_t nentries1 = 2;
	nentries1 = t1->GetEntries();
    
    int f = 6;
    for (Long64_t i=f;i<f+1;i++) {

        t1->GetEntry(i);
        cout<<"runid "<<runid<<endl;
        cout<<"evtid "<<evtid<<endl;

        cout<<X_4<<" "<<Y_4<<" "<<Z_4<<" "<<X_5<<" "<<Y_5<<" "<<Z_5<<endl;
      double L = 230;
    	double  LL = 1*L;

		//TSystem::Load("libGeom");
		TGeoManager *gGeoManager = new TGeoManager("world","for mdc hit count");
		TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
		TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);

		TGeoVolume *all = gGeoManager->MakeBox("all", med, LL, LL, LL);
		gGeoManager->SetTopVolume(all);

		TGeoVolume *tube[1000];
		TGeoCombiTrans *combi_tube[1000];
		double theta_tube[1000];
		double phi_tube[1000];

		//make cone for 1st evt
		
//		double l =100;double r =50;
    double degree = 0.12;
    XYZVector *g1 = new XYZVector(X_4,Y_4,Z_4);
    XYZVector *g2 = new XYZVector(X_5,Y_5,Z_5);
//  
//    x3_ = X_4; y3_ = Y_4; z3_ = Z_4;

    double r1 = 1.0/sqrt(g1->mag2());//*TMath::Sin(dgree); double l1 = gamma_tag1->mag();
    double r2 = 1.0/sqrt(g2->mag2());//*TMath::Sin(dgree); double l2 = gamma_tag2->mag();


    XYZVector *g1_ = new XYZVector(r1*g1->x(),r1*g1->y(),r1*g1->z());
    XYZVector *g2_ = new XYZVector(r2*g2->x(),r2*g2->y(),r2*g2->z());

    XYZVector g3 = *g1_ + *g2_;

    double r_ = 1.0/sqrt(g3.mag2());

    XYZVector *g4 = new XYZVector(g3.x()*r_,g3.y()*r_,g3.z()*r_);
              
    double l = g4->Dot(*g1);//*sqrt(x4*x4+y4*y4+z4*z4);
    double r = TMath::Tan(TMath::ACos(g1_->Dot(*g4))+0.06)*l;

    x3 = l*g3.x()*r_;
    y3 = l*g3.y()*r_;
    z3 = l*g3.z()*r_;


if(x3>0&&y3>0&&z3>0){
              double theta_cone = 3.14/2 + TMath::ATan(x3/z3);
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x<0
if(x3<0&&y3>0&&z3>0){
              double theta_cone = 3.14/2 - TMath::ATan((-1*x3)/z3);
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//y<0
if(x3>0&&y3<0&&z3>0){
              double theta_cone = 3.14/2 + TMath::ATan(x3/z3);
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/(-1*y3));
              phi_cone = phi_cone*360/(2*3.14);
}
//z<0
if(x3>0&&y3>0&&z3<0){
              double theta_cone = -3.14/2 - TMath::ATan(x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x>0
if(x3>0&&y3<0&&z3<0){

              double theta_cone = -3.14/2 - TMath::ATan(x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//y>0
if(x3<0&&y3>0&&z3<0){
              double theta_cone = -3.14/2 + TMath::ATan(-x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//z>0
if(x3<0&&y3<0&&z3>0){
              double theta_cone = 3.14/2 - TMath::ATan(-x3/(1*z3));
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x,y,z<0
if(x3<0&&y3<0&&z3<0){
              double theta_cone = -3.14/2 + TMath::ATan(-x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
              double phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
		//double theta_cone = 98.3528;double phi_cone = 89.6927;
		//double theta_cone = 118.967;double phi_cone = 92.1633;
		//double theta_cone = -91.0076;double phi_cone = 112.478;
		//double theta_cone = -88.1158;double phi_cone = 112.388;
		TGeoVolume *cone = gGeoManager->MakeCone("cone",med,l,0,0,0,r);
		TGeoCombiTrans *combi_cone = new TGeoCombiTrans(z3,x3,y3, 
                                   new TGeoRotation("rot_cone",theta_cone,phi_cone,0));
		all->AddNode(cone,0,combi_cone);

		
    	for(int j =0;j<count;j++){

    		//L = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    
			tube[j] = gGeoManager->MakeTube("tube",med,0,0.00125,L);
			tube[j]->SetLineColor(kRed);

    		x1 = X_1[j];
    		y1 = Y_1[j];
    		z1 = Z_1[j];
    		x2 = X_2[j];
    		y2 = Y_2[j];
    		z2 = Z_2[j];

    		x0 = (x1+x2)/2;
            y0 = (y1+y2)/2;

			//////////////////////// tube
			if((x1-x2)>=0&&(y1-y2)>=0){
			      theta_tube[j] = 3.14/2 + TMath::ATan((x1-x2)/(z1-z2));
			      theta_tube[j] = theta_tube[j]*360/(2*3.14);
			      phi_tube[j] = TMath::ATan(L/(y1-y2));
			      phi_tube[j] = phi_tube[j]*360/(2*3.14);
			}
			if((x1-x2)<=0&&(y1-y2)>=0){
			      theta_tube[j] = 3.14/2 - TMath::ATan(-(x1-x2)/(z1-z2));
			      theta_tube[j] = theta_tube[j]*360/(2*3.14);
			      phi_tube[j] = TMath::ATan(L/(y1-y2));
			      phi_tube[j] = phi_tube[j]*360/(2*3.14);
			}
			if((x1-x2)>=0&&(y1-y2)<=0){
			      theta_tube[j] = 3.14/2 + TMath::ATan((x1-x2)/(z1-z2));
			      theta_tube[j] = theta_tube[j]*360/(2*3.14);
			      phi_tube[j] = 3.14 - TMath::ATan(L/-(y1-y2));
			      phi_tube[j] = phi_tube[j]*360/(2*3.14);
			}
			if((x1-x2)<=0&&(y1-y2)<=0){
			      theta_tube[j] = 3.14/2 - TMath::ATan(-(x1-x2)/(z1-z2));
			      theta_tube[j] = theta_tube[j]*360/(2*3.14);
			      phi_tube[j] = 3.14 - TMath::ATan(L/-(y1-y2));
			      phi_tube[j] = phi_tube[j]*360/(2*3.14);
			}
       		combi_tube[j] = new TGeoCombiTrans(0,x0,y0, 
       				new TGeoRotation("rot_tube",theta_tube[j],phi_tube[j],0));
       		all->AddNode(tube[j], j+1, combi_tube[j]);
    	
    	}

		gGeoManager->CloseGeometry();
		gGeoManager->SetTopVisible();
		all->Draw();

    }



}