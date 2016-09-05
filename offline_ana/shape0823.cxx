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
int shape0823(double x0,double y0,
  double x1,double y1,double z1,
  double x2,double y2,double z2,
  double x3,double y3,double z3,
  double l,double r ){


  	   gSystem->Load("libGenVector");
       #ifdef __CINT__
       gSystem->Load("libMathCore");
       using namespace ROOT::Math;
       #endif

  double L = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  double theta_tube,phi_tube,theta_cone,phi_cone;

//////////////////////// tube

if((x1-x2)>=0&&(y1-y2)>=0){
               theta_tube = 3.14/2 + TMath::ATan((x1-x2)/(z1-z2));
              theta_tube = theta_tube*360/(2*3.14);
               phi_tube = TMath::ATan(L/(y1-y2));
              phi_tube = phi_tube*360/(2*3.14);
}
if((x1-x2)<=0&&(y1-y2)>=0){
               theta_tube = 3.14/2 - TMath::ATan(-(x1-x2)/(z1-z2));
              theta_tube = theta_tube*360/(2*3.14);
               phi_tube = TMath::ATan(L/(y1-y2));
              phi_tube = phi_tube*360/(2*3.14);
}
if((x1-x2)>=0&&(y1-y2)<=0){
               theta_tube = 3.14/2 + TMath::ATan((x1-x2)/(z1-z2));
              theta_tube = theta_tube*360/(2*3.14);
               phi_tube = 3.14 - TMath::ATan(L/-(y1-y2));
              phi_tube = phi_tube*360/(2*3.14);
}
if((x1-x2)<=0&&(y1-y2)<=0){
               theta_tube = 3.14/2 - TMath::ATan(-(x1-x2)/(z1-z2));
              theta_tube = theta_tube*360/(2*3.14);
               phi_tube = 3.14 - TMath::ATan(L/-(y1-y2));
              phi_tube = phi_tube*360/(2*3.14);
}


//////////////////////// cone
//x,y,z>0
if(x3>0&&y3>0&&z3>0){
               theta_cone = 3.14/2 + TMath::ATan(x3/z3);
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x<0
if(x3<0&&y3>0&&z3>0){
               theta_cone = 3.14/2 - TMath::ATan((-1*x3)/z3);
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//y<0
if(x3>0&&y3<0&&z3>0){
               theta_cone = 3.14/2 + TMath::ATan(x3/z3);
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/(-1*y3));
              phi_cone = phi_cone*360/(2*3.14);
}
//z<0
if(x3>0&&y3>0&&z3<0){
               theta_cone = -3.14/2 - TMath::ATan(x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x>0
if(x3>0&&y3<0&&z3<0){

               theta_cone = -3.14/2 - TMath::ATan(x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//y>0
if(x3<0&&y3>0&&z3<0){
               theta_cone = -3.14/2 + TMath::ATan(-x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = TMath::ATan(sqrt(z3*z3+x3*x3)/y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//z>0
if(x3<0&&y3<0&&z3>0){
               theta_cone = 3.14/2 - TMath::ATan(-x3/(1*z3));
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
//x,y,z<0
if(x3<0&&y3<0&&z3<0){
               theta_cone = -3.14/2 + TMath::ATan(-x3/(-1*z3));
              theta_cone = theta_cone*360/(2*3.14);
               phi_cone = 3.14 - TMath::ATan(sqrt(z3*z3+x3*x3)/-y3);
              phi_cone = phi_cone*360/(2*3.14);
}
              ///////////////////////////////////////////////////

              TGeoManager *gGeoManager = new TGeoManager("world","for mdc hit count");

              

              TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
              TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);

              double LL = 5*L;
              TGeoVolume *all = gGeoManager->MakeBox("all", med, LL, LL, LL);
              //TGeoVolume *detector = gGeoManager->MakeTube("detector",med,40,100,L);
              TGeoVolume *tube = gGeoManager->MakeTube("tube",med,0,0.00125,L);
              TGeoVolume *cone = gGeoManager->MakeCone("cone",med,l,0,0,0,r);
              tube->SetLineColor(kRed);
              cone->SetLineColor(kGreen);
              
              
              
              ///////////////////////////////////////////////////
//              cout<<"theta_cone = "<<theta_cone<<";phi_cone = "<<phi_cone<<";"<<endl;

              TGeoCombiTrans *combi_cone = new TGeoCombiTrans(z3,x3,y3, 
                                   new TGeoRotation("rot_cone",theta_cone,phi_cone,0));
              TGeoCombiTrans *combi_tube = new TGeoCombiTrans(0,x0,y0, 
                                   new TGeoRotation("rot_tube",theta_tube,phi_tube,0));
              
              TGeoCombiTrans *combi_tube2 = new TGeoCombiTrans(0,0,0, 
                                   new TGeoRotation("rot_tube",0,0,0));
              TGeoCombiTrans *combi_tube3 = new TGeoCombiTrans(0,0,0, 
                                   new TGeoRotation("rot_tube",90,90,0));
              TGeoCombiTrans *combi_box1 = new TGeoCombiTrans(60,30,30, 
                                   new TGeoRotation("rot_tube",0,0,0));

              gGeoManager->SetTopVolume(all);
              //all->AddNode(detector,1,new TGeoRotation("rot_detector",90,90,0));
              all->AddNode(tube, 1, combi_tube);
              all->AddNode(cone,2,combi_cone);
              //all->AddNode(tube2,3,combi_tube2);
              //all->AddNode(tube3,4,combi_tube3);
              //all->AddNode(mark,5,combi_box1);

              gGeoManager->CloseGeometry();

              //////////////////////////////////////////////////

              gGeoManager->SetTopVisible();
              //all->Draw();
              all->CheckOverlaps();

              TIter next(gGeoManager->GetListOfOverlaps());
              TGeoOverlap* ov = new TGeoOverlap();
              int mm=0;
              while((ov=(TGeoOverlap*)next())){mm++;cout<<ov->GetOverlap()<<endl;}
             // while((ov=(TGeoOverlap*)next())){mm++;cout<<ov->GetOverlap()<<endl;}
              if(mm==1){return 1;}
  
  return 0;
}