#ifndef Physics_Analysis_Psi3770pi_H
#define Physics_Analysis_Psi3770pi_H 
#include "MdcGeomSvc/MdcGeomSvc.h"

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"    //No NTuple!
//#include "VertexFit/ReadBeamParFromDb.h"
#include "PartPropSvc/PartPropSvc.h"  

#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "RootCnvSvc/RootCnvSvc.h"

#include "TROOT.h"
#include "TBenchmark.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IService.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataSvc.h"
#include "RawDataProviderSvc/RawDataProviderSvc.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcPrintSvc/IMdcPrintSvc.h"

#include <vector>
class Psi3770pi : public Algorithm {

	public:
		Psi3770pi(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  

	private:
		IMdcPrintSvc*              m_mdcPrintSvc; 
		RawDataProviderSvc*        	m_rawDataProviderSvc; 

		//MdcGeomSvc* m_mdcsvc;
		//IMdcGeomSvc* mdcsvc;


		IPartPropSvc *p_PartPropSvc;
		HepPDT::ParticleDataTable* m_particleTable;

		void corgen(HepMatrix &, HepVector &, int );
		void corset(HepSymMatrix &, HepMatrix &, int );
		void calibration(RecMdcKalTrack * , HepVector &, int );
		//Declare r0, z0 and cos cut for charged tracks
		Double_t m_vr0cut;
		Double_t m_vz0cut;
		Double_t m_ccoscut;
		Double_t m_ptcut;
		Double_t m_pcut;

		//Declare energy, dphi, dthe cuts for fake gamma's
		Double_t m_isoAngleCut;

		//Declare flag for cut flow
		Int_t m_cutFlow;

		//Declare flag for analysis
		Int_t m_readSig;

		//Declare flag for MCTruth
		Int_t m_mcTruth;

		//Declare type of psi. For Jpsi psitype=1, for psip psitype=0
		Int_t m_psiType;

		Int_t m_useVxfitCut;
		Int_t m_useKmfitCut;
		//Declare name of output file
		std::string m_OutputFileName;
		TFile *saveFile;

		//Define TCutFlow here
		TTree *TCutFlow;
		Int_t cutflow;

		//Define TreeAna here
		TTree *TreeAna;

		//For header
		Int_t runid;
		Int_t evtid;
		Int_t nevt;
		Int_t Ncut0;
		Int_t Ncut1;
		Int_t Ncut2;
		Int_t Ncut3;
		Int_t Ncut4;
		Int_t Ncut5;
		Int_t Ncut6;
		Int_t Ncut7;
		Int_t n_gamma;
		Int_t n_endgamma;
		Int_t n_charged;

		Double_t m_dang[500];
		Double_t m_Rxy[500];
		Double_t m_tdc[500];
		Double_t m_z0[500];


		//For storing 4-mom of particles. Change names below according to your channel
		//After 4c
		TLorentzVector *gamma1;
		TLorentzVector *gamma2;
		TLorentzVector *gamma3;
		    TLorentzVector *pip;
		TLorentzVector *pim;
		TLorentzVector *kp;
		TLorentzVector *km;

		//TLorentzVector gamma_t;
		//TLorentzVector gamma;
		//Before 4c
		TLorentzVector *gamma1_unfitted;
		TLorentzVector *gamma2_unfitted;
		TLorentzVector *gamma3_unfitted;
		TLorentzVector *gamma4_unfitted;
		TLorentzVector *gamma1tmp;
		TLorentzVector *gamma2tmp;
		TLorentzVector *gamma3tmp;

		TLorentzVector *pip_unfitted;
		TLorentzVector *pim_unfitted;
		TLorentzVector *kp_unfitted;
		TLorentzVector *km_unfitted;

		TLorentzVector *ep_unfitted;
		TLorentzVector *em_unfitted;
		//For further cut. Add your own cut here
		TLorentzVector *etap_mc;
		TLorentzVector *pi0_mc;
		TLorentzVector *eta_mc;
		TLorentzVector *gamma2_mc;
		TLorentzVector *ep_mc;
		TLorentzVector *em_mc;
		TLorentzVector *gamma1_mc;
		TLorentzVector *gamma_mc;
		TLorentzVector *Rgamma_mc;
		TLorentzVector *gamma3_mc;
		//		TLorentzVector *gamma;

		Double_t V_xy;
		Double_t vxchisq;
		Double_t kmchisq_4c;
		Double_t Chisq_low;

		Double_t mass_2, mass_22;
		Double_t mass_t, mass_tt;		

		Double_t EOP_p;
		Double_t Edpop;
		Double_t Edpom;
		Double_t EOP_m;
		Double_t Prob_ep;
		Double_t Prob_em;
		Double_t Prob_pip;
		Double_t Prob_pim;
		Double_t P_ep;
		Double_t P_ep2;
		Double_t P_em;
		Double_t P_em2;
		Double_t PD_p;
		Double_t PD_m;
		//                Double_t kmchisq_eta[500];
		//gamma conversion
		Double_t m_xconv1;
		Double_t m_yconv1;
		Double_t m_zconv1;
		Double_t m_rconv1;
		Double_t m_xconv2;
		Double_t m_yconv2;
		Double_t m_zconv2;
		Double_t m_rconv2;
		Double_t m_xiep;
		Double_t m_deltaxy;

		Double_t m_deltaz1;
		Double_t m_deltaz2;

		Double_t m_lep;
		Double_t m_psipair;
		//                Double_t m_dgamma;
		Double_t MEE;
		Double_t m_vx_x;
		Double_t m_vx_y;
		Double_t m_vx_r;
		Double_t m_thetaeg1;
		Double_t m_thetaeg2; 
		Double_t m_cthep;
		Double_t m_ptrkp;
		Double_t m_ptrkm;
		Double_t m_mgamma;
		Double_t m_egamma;
		Double_t m_theta;
		Double_t m_cosTheta;
		Double_t m_phi;

		Double_t m_thedecay;
		Double_t e_rest;
		Double_t e_tot;

		Double_t m_rp;
		Double_t m_re;
		Double_t m_deltaeq;
		Double_t m_case;
		Double_t phi1;
		Double_t phi2;
		Double_t phi3;
             		//Define TMCTruth Here

		TTree *TMCTruth;
//MC truth is stored in TreeAna.
		Int_t m_pdgid[500];
		Int_t m_idxmc;
		Int_t m_motheridx[500];

        TLorentzVector *p_f0_mctruth;

		//special for nhit
		int Nhit;
		Double_t phiE[10000];
		Double_t phiW[10000];
		Double_t phi[10000];


		int	_nSWire; 
		int _nLayer; 
		int _nSlay;  
		int count;
    	

};

#endif 
