#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"


#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/IService.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ITHistSvc.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"


#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "RawDataProviderSvc/RawDataProviderSvc.h"
#include "RawDataProviderSvc/MdcRawDataProvider.h"
#include "RawDataProviderSvc/RawDataProviderSvc.h"


#include "VertexFit/IVertexDbSvc.h"
 
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IJobOptionsSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EmcRawEvent/EmcDigi.h"
#include "EmcRecEventModel/RecEmcHit.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "TofRecEvent/RecTofTrack.h"


#include "EvTimeEvent/RecEsTime.h"

#include "Identifier/Identifier.h"
#include "McTruth/McParticle.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/MucMcHit.h"  
#include "McTruth/McEvent.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"

#include "MucRawEvent/MucDigi.h"
#include "MucRecEvent/RecMucTrack.h"
#include "MucRecEvent/MucRecHit.h"
#include "MucRecEvent/MucRecHitContainer.h"

#include "MdcGeomSvc/MdcGeomSvc.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "TMath.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"

#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "RootCnvSvc/RootInterface.h"
#include "ParticleID/ParticleID.h"
#include "EventNavigator/EventNavigator.h"
#include "PartPropSvc/PartPropSvc.h"


//my include
#include "MdcRecEvent/RecMdcHit.h"
#include "MdcData/MdcHit.h"
#include "MdcData/MdcHitUse.h"
#include "MdcData/MdcRecoHitOnTrack.h"
#include "MdcPrintSvc/IMdcPrintSvc.h"

#include "BesDChain/CDDecayList.h"
#include "BesDChain/CDPhotonList.h"
#include "BesDChain/CDChargedPionList.h"
#include "BesDChain/CDChargedKaonList.h"


#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "Psi3770piAlg/Psi3770pi.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <ctime>
#include "TRandom.h"


#include "MdcPrintSvc/MdcPrintSvc.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataSvc.h"
#include "EventModel/EventHeader.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Identifier/MdcID.h"
#include "RawEvent/RawDataUtil.h"
#include "MdcRawEvent/MdcDigi.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "EventModel/EventModel.h"

#include <iomanip>
#include <iostream>

using namespace std;
using namespace Event;

//const double pi = 3.1415927;
const double mpi = 0.13957;
const double mpi0 = 0.1349766;
const double mk = 0.493677;
const double meta = 0.547853;
const double me = 0.000510998910;
const double mmu = 0.105658367;
const double mp = 0.93827203; 
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double E3770=3.686;//3.686;//3.773;
//const double E2175=3.686;

const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;

using namespace std;
/*
    //total event
int Ncut1 = 0;    //Pass good charged tracks cut
int Ncut2 = 0;    //Pass PID

int Ncut3 = 0;    //Pass Good Photon
int Ncut4 = 0;    //Pass 4C 
int Ncut5 = 0;    //Pass 5C

int Ncut6 = 0;    //Pass E/p
int Ncut7 = 0;    //Pass MUC
*/
/////////////////////////////////////////////////////////////////////////////
Psi3770pi::Psi3770pi(const std::string& name, ISvcLocator* pSvcLocator) :   //do not change this function
  Algorithm(name, pSvcLocator) {


    //Declare the properties  
    declareProperty("OutputFileName",  m_OutputFileName = "out.root");
    declareProperty("PsiType",m_psiType = 1);
    declareProperty("Vr0cut", m_vr0cut=1.0);
    declareProperty("Vz0cut", m_vz0cut=10.0);
    declareProperty("Ccoscut", m_ccoscut=0.93);
    declareProperty("IsoAngleCut", m_isoAngleCut=10.0);
    declareProperty("CutFlow", m_cutFlow = 1);
    declareProperty("ReadSig", m_readSig = 1);
    declareProperty("MCTruth", m_mcTruth = 0);
                
      declareProperty("UseVxfitCut", m_useVxfitCut = 1);
        declareProperty("UseKmfitCut", m_useKmfitCut = 1);     
    }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


StatusCode Psi3770pi::initialize(){

  Ncut0 = 0;
  Ncut1 = 0;
  Ncut2 = 0;
  Ncut3 = 0;
  Ncut4 = 0;
  Ncut5 = 0;
  Ncut6 = 0;
  Ncut7 = 0;
  Nhit = 0;
  phi1 = 0;
  phi2 = 0;
  phi3 = 0;
  count =0;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in initialize()" << endmsg;
  StatusCode status;

  cout<<"initialize"<<endl;


  //***********Initialize the output structure**************

  TString s_OutputFileName(m_OutputFileName);
  //s_OutputFileName.ReplaceAll("[\"","");
  //s_OutputFileName.ReplaceAll("\"]","");
  saveFile = new TFile(s_OutputFileName, "recreate");

  //For cut flow, DO NOT ACTIVATE it except for analysing signal MC
  if(m_cutFlow==1){
    TCutFlow = new TTree("TCutFlow", "cut flow");
    TCutFlow->Branch("cutflow", &cutflow, "cutflow/I");
  }

  //For Analysis
  if(m_readSig==1)
  {
    TreeAna = new TTree("TreeAna", "analysis");

    //To store mc truth.
        if(m_mcTruth==1)
        {
          p_f0_mctruth = new TLorentzVector();

          TreeAna->Branch("p_f0_mctruth",&p_f0_mctruth,32000, 0);
          TreeAna->Branch("indexmc",&m_idxmc,"indexmc/I");
          TreeAna->Branch("pdgid", m_pdgid,"pdgid[indexmc]/I");
          TreeAna->Branch("motheridx",m_motheridx,"motheridx[indexmc]/I");
      }

        //After 4c.
        pip   = new TLorentzVector();
        pim   = new TLorentzVector();
        kp   = new TLorentzVector();
        km   = new TLorentzVector();
        gamma1 = new TLorentzVector();
        gamma2 = new TLorentzVector();
        gamma3 = new TLorentzVector();
        


        //Before 4c.
        pip_unfitted   = new TLorentzVector();
        pim_unfitted   = new TLorentzVector();
        kp_unfitted   = new TLorentzVector();
        km_unfitted   = new TLorentzVector();
        gamma1_unfitted = new TLorentzVector();
        gamma2_unfitted = new TLorentzVector();
        gamma3_unfitted = new TLorentzVector();
        gamma4_unfitted = new TLorentzVector();


               
        //Header
        TreeAna->Branch("runid", &runid, "runid/I");
        TreeAna->Branch("evtid", &evtid, "evtid/I");


        TreeAna->Branch("n_charged",&n_charged,"n_charged/I");
        TreeAna->Branch("n_gamma",&n_gamma,"n_gamma/I");

        //TreeAna->Branch("n_gamma_e",&n_gamma_e,"n_gamma_e/I");

        TreeAna->Branch("n_endgamma",&n_endgamma,"n_endgamma/I"); 
        TreeAna->Branch("m_dang",m_dang,"m_dang[n_gamma]/D");
        TreeAna->Branch("m_Rxy",m_Rxy,"m_Rxy[n_charged]/D");
        TreeAna->Branch("m_tdc",m_tdc,"m_tdc[n_gamma]/D");
        TreeAna->Branch("m_z0",m_z0,"m_z0[n_charged]/D");
          
        TreeAna->Branch("Chisq_low", &Chisq_low, "Chisq_low/D");  
        TreeAna->Branch("vxchisq", &vxchisq, "vxchisq/D");
        //Particles

              //After 4c
        TreeAna->Branch("gamma1",&gamma1,32000, 0);
        TreeAna->Branch("gamma2",&gamma2,32000, 0); 
        TreeAna->Branch("gamma3",&gamma3,32000, 0);
        TreeAna->Branch("e_rest",&e_rest,"e_rest/D");
        TreeAna->Branch("e_tot",&e_tot,"e_tot/D");
        //TreeAna->Branch("mdchit",&mdchit,"mdchit/D");


        TreeAna->Branch("m_thedecay",&m_thedecay,"m_thedecay/D");



        TreeAna->Branch("mass_t",&mass_t,"mass_t/D");

        TreeAna->Branch("mass_2",&mass_2,"mass_2/D");

        TreeAna->Branch("phi1",&phi1,"phi1/D");
        TreeAna->Branch("phi2",&phi2,"phi2/D");
        TreeAna->Branch("phi3",&phi3,"phi3/D");
        TreeAna->Branch("count",&count,"count/I");

        TreeAna->Branch("kp",&kp,32000,0);
        TreeAna->Branch("km",&km,32000,0);
        //Before 4c
        TreeAna->Branch("gamma1_unfitted",&gamma1_unfitted,32000, 0);
        TreeAna->Branch("gamma2_unfitted",&gamma2_unfitted,32000, 0);
        TreeAna->Branch("gamma3_unfitted",&gamma3_unfitted,32000, 0);
        //TreeAna->Branch("gamma4_unfitted",&gamma4_unfitted,32000, 0);
        TreeAna->Branch("kp_unfitted",&kp_unfitted,32000,0);
        TreeAna->Branch("km_unfitted",&km_unfitted,32000,0);

        TreeAna->Branch("Nhit",&Nhit,32000,0);
        TreeAna->Branch("phiE",&phiE,"phiE[10000]/D");
        TreeAna->Branch("phiW",&phiE,"phiW[10000]/D");
        TreeAna->Branch("phi",&phiE,"phi[10000]/D");


  }



  //********************************************************
       
  static const bool CREATEIFNOTTHERE(true);
  StatusCode PartPropStatus = Gaudi::svcLocator()->service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
  if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
    std::cerr << "Could not initialize Particle Properties Service" << std::endl;
    return StatusCode::FAILURE;
  }     
  m_particleTable = p_PartPropSvc->PDT();

  //cout<<"ncut0++ in initialize before MdcPrintSvc "<<Ncut0<<endl;


  IMdcPrintSvc* imdcPrintSvc;
  StatusCode sc = service ("MdcPrintSvc", imdcPrintSvc);
  m_mdcPrintSvc = dynamic_cast<MdcPrintSvc*> (imdcPrintSvc);
  if ( sc.isFailure() ){
      log << MSG::FATAL << "Could not load MdcPrintSvc!" << endreq;
      return StatusCode::FAILURE;
  }

  IRawDataProviderSvc* irawDataProviderSvc;
    sc = service ("RawDataProviderSvc", irawDataProviderSvc);
    m_rawDataProviderSvc = dynamic_cast<RawDataProviderSvc*> (irawDataProviderSvc);
    if ( sc.isFailure() ){
      log << MSG::FATAL << "Could not load RawDataProviderSvc!" << endreq;
      return StatusCode::FAILURE;
    }

  



  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  //cout<<"ncut0++ in initialize "<<Ncut0<<endl;
  return StatusCode::SUCCESS;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Psi3770pi::execute() {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");  
  int runNo=eventHeader->runNumber();
  int event=eventHeader->eventNumber();
  
  log << MSG::DEBUG <<"run, evtnum = "
    << runNo << " , "
    << event <<endreq;

  Ncut0++;            //Total Event counter
  
  if(m_cutFlow ==1){
    cutflow=0;
    TCutFlow->Fill();
  }

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  //  log << MSG::INFO << "get event tag OK" << endreq;
  log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
    << evtRecEvent->totalCharged() << " , "
    << evtRecEvent->totalNeutral() << " , "
    << evtRecEvent->totalTracks() <<endreq;

  if(Ncut0%100 == 0)
  {
    cout<<"Processing "<<Ncut0<<"th event..."<<endl;
    cout<<"Event Level Cut flow: total events: "<<Ncut0<<" charged track selection: "<<Ncut1<<" | photon selection: "<<Ncut2<<"  PID: "<<Ncut3<<" | vertex fit "<<Ncut4<<" 4c fit "<<Ncut5<<endl;
  }

  //******************MCTruth************************

    if(m_mcTruth==1&&runNo<0)
    {

        int m_numParticle = 0; 
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

  if (!mcParticleCol){    
      
      std::cout << "Could not retrieve McParticelCol" << std::endl;
      return StatusCode::FAILURE;
  }    
  else{    
      int rootIndex = 0; 
      m_pdgid[m_numParticle] = 0; 
      m_motheridx[m_numParticle] = -1;
      m_numParticle ++;

      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++){
      //if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        
        int mcidx;
        
        if((*iter_mc)->primaryParticle())mcidx = 0;
        else mcidx = ((*iter_mc)->mother()).trackIndex() + 1;
        
        int pdgid = (*iter_mc)->particleProperty();
        m_pdgid[m_numParticle] = pdgid;
        m_motheridx[m_numParticle] = mcidx;
        m_numParticle ++;
      }
  }
  m_idxmc = m_numParticle;


    }
    



  
  //******************Hits part - new **************
  int N=0; 
  unsigned _l[10000],_w[10000];

 //   Gaudi::svcLocator()->service("MdcPrintSvc",m_mdcPrintSvc);
 
  //m_mdcPrintSvc->printDigi(); 

    StatusCode scc;
    IService* ser;
    scc = Gaudi::svcLocator()->getService("MdcGeomSvc",ser);
    if (!scc.isSuccess())std::cout <<" MdcDetector::Could not open Geometry Service"<<std::endl;
    MdcGeomSvc* mdcsvc = dynamic_cast<MdcGeomSvc*> (ser);
    if(!mdcsvc) std::cout <<"MdcDetector::Could not open Geometry Service"<<std::endl;


/*
  const std::string& name = " ";
  ISvcLocator* svcloc = NULL ;
    MdcGeomSvc fGeom(name, svcloc);
    fGeom.initialize();
    cout <<"dump"<<endl;
    fGeom.Dump();
    cout<<"Hi, event routine is running"<<endl;
    cout<<"1st wire id:"<<fGeom.Wire(0)->Id()<<endl;
    cout<<"2nd wire lyr id:"<<fGeom.Wire(1)->Lyr()->Id()<<endl;
    cout<<"6860th wire lyr id:"<<fGeom.Wire(6859)->Lyr()->Id()<<endl;
    fGeom.finalize();
*/



    SmartDataPtr<Event::MdcMcHitCol> mcMdcMcHitCol(eventSvc(),"/Event/MC/MdcMcHitCol");

   // if (!mcMdcMcHitCol) {
  //    cout<<"not find MdcMCHit, it is data!"<<endl;
      log << MSG::INFO << "Could not find Mdc_Mc_Hit, start at print digi" << endreq;
    
  
      MdcDigiVec mdcDigiVec = m_rawDataProviderSvc->getMdcDigiVec();
      MdcDigiCol::iterator iter = mdcDigiVec.begin();

      for (int iDigi=0;iter!= mdcDigiVec.end(); iter++,iDigi++ ) {
        _l[iDigi] = MdcID::layer((*iter)->identify());
        _w[iDigi] = MdcID::wire((*iter)->identify());
        N++;

      }

  //  }
  /*
    else{
//     if(m_debug>2) 
      cout<<"it is mc!"<<endl;
      //std::cout<<"nTruthHit="<< mcMdcMcHitCol->size()<<" Truth hits :   "<<std::endl;
      Event::MdcMcHitCol::iterator iter_mchit = mcMdcMcHitCol->begin();
     
      for (int i = 0; iter_mchit != mcMdcMcHitCol->end(); iter_mchit++,i++) {
          const Identifier id= (*iter_mchit)->identify();
          _l[i]   = MdcID::layer(id);
          _w[i]   = MdcID::wire(id);
          N++;
      }
    }
*/

/*
    _nSWire = mdcsvc->getWireSize();
    _nLayer = mdcsvc->getLayerSize();
    _nSlay  = mdcsvc->getSuperLayerSize();
*/
  //cout<<__LINE__<<endl;

    Nhit=N;
    cout<<"Nhit "<<Nhit<<" N "<<N<<endl;

/*
    int ihit=0;

    for(int j=0; j<N; j++){

      //cout<<__LINE__<<endl;

        const MdcGeoWire *geowir=mdcsvc->Wire(_l[j],_w[j]);

        //if(((8<=_l[j]&&_l[j]<20)||(36<=_l[j]&&_l[j]<43))){
          if(8<=_l[j]&&_l[j]<43){
    //   if(8<=_l[j]){

          HepPoint3D eastP = geowir->Backward()/10.0;
          HepPoint3D westP = geowir->Forward() /10.0;
          
          phiE[ihit]=eastP.phi();
          phiW[ihit]=westP.phi();
          phi[ihit]=(eastP.phi()+westP.phi())/2.0;

          cout<<"Layer: "<<_l[j]<<endl;
          cout<<"phiE: "<<phiE[ihit]<<endl;
          cout<<"phiW: "<<phiW[ihit]<<endl;
          cout<<"phi: "<<phi[ihit]<<endl;

          ihit++;
      }
    }
*/


  //*************Global Event Parameters************
  
  //do not change below
  double ecms;
  int psi;
  double ESpread;

  if(m_psiType == 1){
    ecms=E3770;    
  //  psi=443;
    ESpread=0.0013;
  }


  HepLorentzVector cms(0.011*ecms, 0., 0., ecms);    //for 4C fit

    //get the track collection 
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

  //************Global Vectors**********************

  //do not change below
  Vint iGood, iGamma;
    //Vint ipip,ipim,ikp,ikm;

  int nGood = 0; 
  int nCharge = 0;
  int nGamma = 0;

  iGood.clear();
  iGamma.clear();

  //*****************Primary Vertex*****************

  Hep3Vector xorigin(0,0,0);
    HepSymMatrix VtxErr(3,0);
  IVertexDbSvc*  vtxsvc;

  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    //cout<<"valid for vertex"<<endl;
    double* dbv = vtxsvc->PrimaryVertex(); 
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
                
        VtxErr[0][0] = vv[0]*vv[0];
      VtxErr[1][1] = vv[1]*vv[1];
      VtxErr[2][2] = vv[2]*vv[2];
  }



//*****************************************************
//Good Charged Track Selection
//****************************************************

  n_charged =0;
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
  HepPoint3D OP(0,0,0);



/*
  for(int i = 0; i< evtRecEvent->totalTracks(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    cout<<(*itTrk)->trackId()<<endl;
  }
  cout<<"-------------------------------------"<<evtRecEvent->totalCharged()<<endl;
*/
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    if(i >= evtRecTrkCol->size()) break;

    //cout<<"-----------track----------- "<<i<<" event "<<Ncut0<<endl;
    
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;

    if(!(*itTrk)->isMdcTrackValid()){
      cout<<"<<<<<<<<<<<<<<<<<<invalid<<<<<<<<<<<<<<<<"<<endl;
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

      cout<< " theta "<<mdcTrk->theta()<<" phi "<<mdcTrk->phi()<<" phi0helix "<<mdcTrk->helix(1)<<" first layer "<<mdcTrk->firstLayer()<<" last layer "<<mdcTrk->lastLayer()<<" nlayer "<<mdcTrk->nlayer()<<endl;
      continue;
    }

    if(!(*itTrk)->isMdcKalTrackValid()) continue;

    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();



    double pch=mdcTrk->p();
    double ptch=mdcTrk->pxy();
    double x0=mdcTrk->x();
    double y0=mdcTrk->y();
    double z0=mdcTrk->z();
    double phi0=mdcTrk->helix(1);
    double xv=xorigin.x();
    double yv=xorigin.y();
    double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);


    //cout<<"-----------end--------------"<<endl;
    //raw vertax
    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
  //  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
    VFHelix helixip(point0,a,Ea); 
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
    double  Rvphi0=vecipa[1];
    double ccos=cos(mdcTrk->theta());
    //
    //Good Charged Track Cut is set here
    if(fabs(ccos) > m_ccoscut) continue;
    if(fabs(Rvz0) >= m_vz0cut ) continue;
    if(fabs(Rvxy0) >= m_vr0cut) continue;
  
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    //if(mdcKalTrk->charge() == 0) continue;
              
        m_Rxy[n_charged] = Rvxy0;
        m_z0[n_charged]  = Rvz0;
        n_charged ++;
    iGood.push_back((*itTrk)->trackId());
    nCharge += mdcKalTrk->charge();
  }

  nGood = iGood.size(); 
  log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
  
  //cout<<"nGood "<<nGood<<", nCharge"<<nCharge<<endl;

  if((nGood != 0) || (nCharge != 0)) return StatusCode::SUCCESS;

  Ncut1++;  

  if(m_cutFlow ==1)
  {
    cutflow=1;
    TCutFlow->Fill();
  }

//**************************************
//  Good Photon selection
//*************************************

        n_gamma=0;
        n_endgamma=0;
        e_tot = 0;

  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++){
  if(i >= evtRecTrkCol->size()) break;

//Emcshower: electron and photon    information about incident angle and deposited energy

    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;

    RecEmcShower *emcTrk = (*itTrk)->emcShower();

    //cout<<(*itTrk)->trackId()<<endl;;
    //cout<<"_____________"<<endl;
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

        double eraw = emcTrk->energy();
    // good photon cut will be set here
    
    double the = emcpos.theta();
    double e_threshold = 10.0; //Mev

    if(fabs(cos(the)) < 0.8)   e_threshold = 0.025;
    //else if((fabs(cos(the)) > 0.86) && (fabs(cos(the)) < 0.92)) e_threshold = 0.050;

    if(eraw < e_threshold) continue;
        //if(fabs(dang) < m_isoAngleCut) continue;
    if(emcTrk->time()>14  || emcTrk->time()<0) continue;

    
    if(e_threshold ==0.050 ) n_endgamma++;
                //m_dang[n_gamma] =dang;
                m_tdc[n_gamma]= ( emcTrk->time());
                n_gamma++;
    

    iGamma.push_back((*itTrk)->trackId());
    e_tot = e_tot + eraw;
    //cout<<Ncut0<<" eraw "<<eraw<<endl;

        
  }

  nGamma = iGamma.size();

  log << MSG::DEBUG << "num Good Photon " << nGamma  << " , " <<evtRecEvent->totalNeutral()<<endreq;
  
  //cout<<"nGamma = "<< nGamma <<endl;

  if(nGamma < 3) return StatusCode::SUCCESS;

  Ncut2++; 

  if(m_cutFlow ==1)
  {
    cutflow=2;
    TCutFlow->Fill();
  } 

// Assign 4-momentum to each photon
  
    HepLorentzVector ptotal(0,0,0,0);
  Vp4 pGamma;
  pGamma.clear();

  for(int i = 0; i < nGamma; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGamma[i]; 
    RecEmcShower* emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();


    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z()); 
    Hep3Vector gammaDirection = emcpos - xorigin;
    double phi = gammaDirection.phi();
    double the = gammaDirection.theta();


    HepLorentzVector ptrk;
    ptrk.setPx(eraw*sin(the)*cos(phi));
    ptrk.setPy(eraw*sin(the)*sin(phi));
    ptrk.setPz(eraw*cos(the));
    ptrk.setE(eraw);
    ptotal=ptotal+ptrk;
    pGamma.push_back(ptrk);
  }


//**************************************
// Kinematic Fit
//**************************************



    //cout<<"----------end of event before fit------"<< Ncut0 <<endl;

    nevt=0;
    
    Chisq_low=50;
    //Chisq_low=200;
    
    mass_t = 0;
    mass_2 = 0.0;
    e_rest = 0.0;


        TLorentzVector gamma;
        TLorentzVector gamma_t;

        KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();

        for(int i1 = 0; (i1 < nGamma-2); i1++){
          
          RecEmcShower* gammaTrk1 = (*(evtRecTrkCol->begin()+iGamma[i1]))->emcShower();
          
        for(int i2 = i1+1; i2 < (nGamma-1); i2++){
            
            RecEmcShower* gammaTrk2 = (*(evtRecTrkCol->begin()+iGamma[i2]))->emcShower();

        for(int i3 = i2+1; i3 < nGamma; i3++){

                RecEmcShower* gammaTrk3 = (*(evtRecTrkCol->begin()+iGamma[i3]))->emcShower();
           
                kmfit->init();

          kmfit->AddTrack(0, 0.0, gammaTrk1);
          kmfit->AddTrack(1, 0.0, gammaTrk2);         
          kmfit->AddTrack(2, 0.0, gammaTrk3);


          kmfit->AddFourMomentum(0, cms);


          bool okKmfit = kmfit->Fit();
          if((m_useKmfitCut==1)&&(!okKmfit)) continue;

          kmchisq_4c = kmfit->chisq();

          //cout<<"okfit"<<endl;
          if((kmchisq_4c>Chisq_low)||(kmchisq_4c<0.01))continue;

          Chisq_low = kmchisq_4c;

          HepLorentzVector gamma1_ = kmfit->pfit(0);
          HepLorentzVector gamma2_ = kmfit->pfit(1);
          HepLorentzVector gamma3_ = kmfit->pfit(2);


          Hep3Vector emcpos1(gammaTrk1->x(), gammaTrk1->y(), gammaTrk1->z()); 
          Hep3Vector gammaDirection1 = emcpos1 - xorigin;
          phi1 = gammaDirection1.phi();
        
          Hep3Vector emcpos2(gammaTrk2->x(), gammaTrk2->y(), gammaTrk2->z()); 
          Hep3Vector gammaDirection2 = emcpos2 - xorigin;
          phi2 = gammaDirection2.phi();
        
          Hep3Vector emcpos3(gammaTrk3->x(), gammaTrk3->y(), gammaTrk3->z()); 
          Hep3Vector gammaDirection3 = emcpos3 - xorigin;
          phi3 = gammaDirection3.phi();

          
          int m = 0;
          double dphi1,dphi2,dphi3,tag1,tag2,tag1p,tag1m,tag2p,tag2m;
          dphi1 = phi1 - phi2;if(dphi1<0){m=1;dphi1=-dphi1;}
          dphi2 = phi1 - phi3;if(dphi2<0){m=1;dphi2=-dphi2;}
          dphi3 = phi2 - phi3;if(dphi3<0){m=1;dphi3=-dphi3;}
          
          if(dphi1<dphi2 && dphi1<dphi3){if(m==1){tag1=phi1;tag2=phi2;}else{tag1=phi2;tag2=phi1;}}
          if(dphi2<dphi1 && dphi2<dphi3){if(m==1){tag1=phi1;tag2=phi3;}else{tag1=phi3;tag2=phi1;}}
          if(dphi3<dphi2 && dphi3<dphi2){if(m==1){tag1=phi2;tag2=phi3;}else{tag1=phi3;tag2=phi2;}}

          tag1p = tag1 + 0.06;
          tag1m = tag1 - 0.06;


          tag2p = tag2 + 0.06;
          tag2m = tag2 - 0.06;

          tag1 = tag1 - 0.06;
          tag2 = tag2 + 0.6;//0.06 = 3 dgree;



          int ihit=0;
          count = 0;
          for(int j=0; j<N; j++){

            const MdcGeoWire *geowir=mdcsvc->Wire(_l[j],_w[j]);
            if(8<=_l[j]&&_l[j]<43){

              HepPoint3D eastP = geowir->Backward()/10.0;
              HepPoint3D westP = geowir->Forward() /10.0;
          
              phiE[ihit]=eastP.phi();
              phiW[ihit]=westP.phi();
              phi[ihit]=(eastP.phi()+westP.phi())/2.0;

              //if(tag1<=phi[ihit]){
              //if(tag1<=phi[ihit]&&phi[ihit]<=tag2){
              if(  (tag1m<=phi[ihit]&&phi[ihit]<=tag1p) ||  (tag2m<=phi[ihit]&&phi[ihit]<=tag2p) ) {
                count++;
              }
              //cout<<"phi "<<phi[ihit]<<" tag1m "<<tag1m<<" tag1p "<<tag1p<<" tag2m "<<tag2m<<" tag2p "<<tag2p<<endl;

              ihit++;
            }
          }
          if(count!=0){
          cout<<"ihit "<<ihit<<endl;
          cout<<"count "<<count<<endl;
          }


          //After 4c
          gamma1->SetPxPyPzE(gamma1_.px(), gamma1_.py(), gamma1_.pz(), gamma1_.e());
          gamma2->SetPxPyPzE(gamma2_.px(), gamma2_.py(), gamma2_.pz(), gamma2_.e());
          gamma3->SetPxPyPzE(gamma3_.px(), gamma3_.py(), gamma3_.pz(), gamma3_.e());


          //Before 4c
          //Photon:after  Assign 4-momentum to each photon.

          gamma1_unfitted->SetPxPyPzE(pGamma[i1].px(), pGamma[i1].py(), pGamma[i1].pz(), pGamma[i1].e());
          gamma2_unfitted->SetPxPyPzE(pGamma[i2].px(), pGamma[i2].py(), pGamma[i2].pz(), pGamma[i2].e());
          gamma3_unfitted->SetPxPyPzE(pGamma[i3].px(), pGamma[i3].py(), pGamma[i3].pz(), pGamma[i3].e());
          //gamma4_unfitted->SetPxPyPzE(pGamma[i4].px(), pGamma[i4].py(), pGamma[i4].pz(), pGamma[i4].e());


          gamma_t = *gamma2_unfitted + *gamma3_unfitted + *gamma1_unfitted;
          gamma = *gamma2_unfitted + *gamma3_unfitted;

          mass_2 = gamma.M();//sqrt(gamma.E()*gamma.E()-gamma.Pz()*gamma.Pz()-gamma.Py()*gamma.Py()-gamma.Px()*gamma.Px());
          mass_t = gamma_t.M();//sqrt(gamma_t.E()*gamma_t.E()-gamma_t.Pz()*gamma_t.Pz()-gamma_t.Py()*gamma_t.Py()-gamma_t.Px()*gamma_t.Px());


          e_rest = e_tot - gamma_t.E();

          m_thedecay = (gamma2_unfitted->E() - gamma3_unfitted->E())/gamma.P();
          if(m_thedecay<0)m_thedecay = - m_thedecay;


          nevt++;

            }//third photon 
         }//second photon
      }//first photon

  runid=runNo;
  evtid=event;

  //cout<<__LINE__<<endl;


  //cout<< "nevt: "<<nevt<<endl;
  if(nevt==0) return StatusCode::SUCCESS;
  
  Ncut5++;

  if(m_cutFlow ==1)
  {
    cutflow=5;
    TCutFlow->Fill();
  } 

  TreeAna->Fill();
  //cout<<"----------end of event ------"<< Ncut0 <<endl;

//*************************************************************************
  return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Psi3770pi::finalize()
{
  cout << "In finalize()..." << endl;

  saveFile->cd();
  if(m_cutFlow==1) TCutFlow->Write();
  if(m_readSig==1) TreeAna->Write();
  saveFile->Close();

  cout << "Total Event Number: " << Ncut0 << endl;
  cout << "After charged track cut: " << Ncut1 << endl;
  cout << "After Good Photon: " << Ncut2 << endl;
  cout << "After PID: " << Ncut3 << endl;
  cout << "After vertex fit :" << Ncut4 << endl;
    cout << "After 4C:" << Ncut5 << endl;
    cout << "After E/P:" << Ncut6 << endl;
    cout << "After MUC:" << Ncut7 << endl;

  return StatusCode::SUCCESS;
}
