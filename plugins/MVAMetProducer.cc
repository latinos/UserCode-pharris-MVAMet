

// system include files
#include <memory>
#include <cmath>
#include <algorithm>
#include <TLorentzVector.h>

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/CommonMETData.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "pharris/MVAMet/plugins/MVAMetProducer.h"

using namespace edm;
using namespace std;
using namespace reco;

MVAMetProducer::MVAMetProducer(const edm::ParameterSet& iConfig) {
  produces<reco::PFMETCollection>();
  fCorrJetName    = iConfig.getParameter<edm::InputTag>("CorrJetName");
  fUnCorrJetName  = iConfig.getParameter<edm::InputTag>("JetName");
  fPFCandName     = iConfig.getParameter<edm::InputTag>("PFCandidateName");
  fVertexName     = iConfig.getParameter<edm::InputTag>("VertexName");
  fRhoName        = iConfig.getParameter<edm::InputTag>("RhoName");
  fJetPtMin       = iConfig.getParameter<double>       ("JetPtMin");
  fDZMin          = iConfig.getParameter<double>       ("dZMin");

  fPUJetIdAlgo    = new PileupJetIdAlgo(iConfig);
  fMVAMet         = new MVAMet(fDZMin);
  fMVAMet         ->Initialize(iConfig,
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmet_52.root"))),        //U
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmetphi_52.root"))),     //U Phi
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmetu1cov_52.root"))),   //U1 Cov
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmetu2cov_52.root"))) //U2 Cov
			       );

}
MVAMetProducer::~MVAMetProducer() { 
  delete fMVAMet;
}
void MVAMetProducer::beginJob() { }

void MVAMetProducer::endJob() { } 

void MVAMetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {  
  //Uncorrected Jets
  Handle<PFJetCollection>       lHUCJets;
  iEvent.getByLabel(fUnCorrJetName, lHUCJets);
  PFJetCollection               lUCJets = *lHUCJets;

  //Corrected Jets
  Handle<PFJetCollection>       lHCJets;
  iEvent.getByLabel(fCorrJetName  , lHCJets);
  PFJetCollection               lCJets = *lHCJets;

  //Get pfCandidates
  Handle<PFCandidateCollection> lHCands;
  iEvent.getByLabel(fPFCandName   , lHCands);
  PFCandidateCollection         lCands = *lHCands;

  // vertices    
  Handle<reco::VertexCollection> lHVertices;
  iEvent.getByLabel(fVertexName      , lHVertices); 
  VertexCollection lVertices = *lHVertices;
  Vertex *lPV = 0; if(lVertices.size() > 0) lPV = &lVertices[0]; 

  //Get PF Met
  Handle<reco::PFMETCollection> lHPFMet;
  iEvent.getByLabel("pfMet"      , lHPFMet); 
  PFMETCollection lPFMET = *lHPFMet;

  //Get Rho
  Handle<double>                                        lHRho;
  iEvent.getByLabel(fRhoName                          , lHRho);
  double lRho = *lHRho;
  
  //Make Generic Objects
  std::vector<LorentzVector >                                     lVisible;
  std::vector<          std::pair<LorentzVector,double> >         lPFInfo;  makeCandidates(lPFInfo, lCands,lPV);
  std::vector<MetUtilities::JetInfo>                              lJetInfo; makeJets      (lJetInfo,lUCJets,lCJets,lVertices,lRho);
  std::vector<Vector>                                             lVtxInfo; makeVertices  (lVtxInfo,lVertices);
 
  //Dummy visible stuff
  double lDummyPt0  = 10; 
  double lDummyPhi0 = 0;
  double lDummyEta0 = 1;
  double lDummyPt1  = 10;
  double lDummyPhi1 = 2.; 
  double lDummyEta1 = -1;
  TLorentzVector lD0; lD0.SetPtEtaPhiM(lDummyPt0,lDummyEta0,lDummyPhi0,0.);
  TLorentzVector lD1; lD1.SetPtEtaPhiM(lDummyPt1,lDummyEta1,lDummyPhi1,0.);

  LorentzVector lVis0; lVis0.SetCoordinates(lD0.Px(),lD0.Py(),lD0.Pz(),lD0.E());
  LorentzVector lVis1; lVis1.SetCoordinates(lD1.Px(),lD1.Py(),lD1.Pz(),lD1.E());
  lVisible.push_back(lVis0);
  lVisible.push_back(lVis1);
  //Calculate the MVA
  std::pair<LorentzVector,TMatrixD> lMVAMetInfo = fMVAMet->GetMet(lVisible,lJetInfo,lPFInfo,lVtxInfo,false);
  std::cout << "Met---> " << lMVAMetInfo.first.pt() << " -- " << lMVAMetInfo.first.phi() 
	    << " Cov matrix " << lMVAMetInfo.second(0,0)  << " -- " << lMVAMetInfo.second(0,1)  << " -- " << lMVAMetInfo.second(1,0)  << " -- " << lMVAMetInfo.second(1,1) << std::endl;

  //Add a PF Met object and put it in the event
  PFMET lDummy;
  PFMET lMVAMet(lDummy.getSpecific(),lPFMET.at(0).sumEt(),lMVAMetInfo.first,lPV->position()); //Use PFMET sum Et
  std::auto_ptr<reco::PFMETCollection> lPFMetColl;
  lPFMetColl.reset     ( new reco::PFMETCollection);
  lPFMetColl->push_back( lMVAMet);
  iEvent.put( lPFMetColl );
}
void MVAMetProducer::makeJets(std::vector<MetUtilities::JetInfo> &iJetInfo,PFJetCollection &iUCJets,PFJetCollection &iCJets,VertexCollection &iVertices,double iRho) { 
  for(int i0   = 0; i0 < (int) iUCJets.size(); i0++) {   // uncorrecte jets collection                                           
    const PFJet       *pUCJet = &(iUCJets.at(i0));
    for(int i1 = 0; i1 < (int) iCJets .size(); i1++) {   // corrected jets collection                                         
      const PFJet     *pCJet  = &(iCJets.at(i1));
      if(       pUCJet->jetArea() != pCJet->jetArea()                  ) continue;
      if( fabs(pUCJet->eta() - pCJet->eta())         > 0.01            ) continue;
      if( pUCJet->pt()                               < fJetPtMin       ) continue;
      //if( fabs(pCJet->eta())                         > 4.99            ) continue;
      if( !passPFLooseId(pUCJet)                                       ) continue;
      double lJec = pCJet ->pt()/pUCJet->pt();
      double lMVA = jetMVA(pUCJet,lJec,iVertices[0],iVertices,false);
      double lNeuFrac = (pUCJet->neutralEmEnergy()/pUCJet->energy() + pUCJet->neutralHadronEnergy()/pUCJet->energy());
      MetUtilities::JetInfo pJetObject; 
      pJetObject.p4       = pCJet->p4(); 
      pJetObject.mva      = lMVA;
      pJetObject.neutFrac = lNeuFrac;
      //Following discussion with Jet MET convenors < 10 GeV use rho correction only below 10
      if(pUCJet->pt() < 10) { 
	TLorentzVector lCorrPt ;  lCorrPt.SetPtEtaPhiM(max(pUCJet->pt()-pCJet->jetArea()*iRho,0.),pCJet->eta(),pCJet->phi(),pCJet->mass());
	pJetObject.p4.SetCoordinates(lCorrPt.Px(),lCorrPt.Py(),lCorrPt.Pz(),lCorrPt.E());
      }
      //pJetObject.p4       = pUCJet->p4(); 
      iJetInfo.push_back(pJetObject);
      break;
    }
  }
}
void MVAMetProducer::makeCandidates(std::vector<std::pair<LorentzVector,double> > &iPFInfo,PFCandidateCollection &iCands,Vertex *iPV) { 
  for(int i0 = 0; i0 < (int)iCands.size(); i0++) {
    const PFCandidate*  pflowCand = &(iCands.at(i0));
    double pDZ = -999;
    if(iPV != 0) pDZ  = pfCandDz(pflowCand,iPV); //If there is no track return negative number -999
    //LorentzVector pVec; pVec.SetCoordinates(pflowCand->pt(),pflowCand->eta(),pflowCand->phi(),pflowCand->mass());
    std::pair<LorentzVector,double> pPFObject(pflowCand->p4(),pDZ);
    iPFInfo.push_back(pPFObject);
  }
}
void MVAMetProducer::makeVertices(std::vector<Vector>        &iPVInfo,VertexCollection &iVertices) { 
  for(int i0    = 0; i0 < (int)iVertices.size(); i0++) {
    const Vertex       *pVertex = &(iVertices.at(i0));
    Vector pVec; pVec.SetCoordinates(pVertex->x(),pVertex->y(),pVertex->z());
    iPVInfo.push_back(pVec);
  }
}
bool MVAMetProducer::passPFLooseId(const PFJet *iJet) { 
  //std::cout << " ==> " << iJet->pt() << " - " << iJet->eta() << " : " << iJet->energy() << " - " << iJet->neutralEmEnergy() << " - " << iJet->nConstituents() << " - " << iJet->chargedHadronEnergy() << " _ " << iJet->chargedEmEnergy() << " - " << iJet->chargedMultiplicity() << std::endl;
  if(iJet->energy()== 0)                                  return false;
  if(iJet->neutralHadronEnergy()/iJet->energy() > 0.99)   return false;
  if(iJet->neutralEmEnergy()/iJet->energy()     > 0.99)   return false;
  if(iJet->nConstituents() <  2)                          return false;
  if(iJet->chargedHadronEnergy()/iJet->energy() <= 0 && fabs(iJet->eta()) < 2.4 ) return false;
  if(iJet->chargedEmEnergy()/iJet->energy() >  0.99  && fabs(iJet->eta()) < 2.4 ) return false;
  if(iJet->chargedMultiplicity()            < 1      && fabs(iJet->eta()) < 2.4 ) return false;
  return true;
}
double MVAMetProducer::pfCandDz(const PFCandidate* iPFCand, const Vertex *iPV) { 
  double lDz = -999;
  if(iPFCand->trackRef().isNonnull())    lDz = fabs(iPFCand->   trackRef()->dz(iPV->position()));
  if(iPFCand->gsfTrackRef().isNonnull()) lDz = fabs(iPFCand->gsfTrackRef()->dz(iPV->position()));
  return lDz;
}
double MVAMetProducer::jetMVA (const PFJet *iCorrJet,double iJec, const Vertex iPV, const reco::VertexCollection &iAllvtx,bool iPrintDebug) { 
  PileupJetIdentifier lPUJetId =  fPUJetIdAlgo->computeIdVariables(iCorrJet,iJec,&iPV,iAllvtx,true);
  if(iPrintDebug) { std::cout << "Debug Jet MVA: "
			      << lPUJetId.nvtx()      << " "
			      << iCorrJet->pt()       << " "
			      << lPUJetId.jetEta()    << " "
			      << lPUJetId.jetPhi()    << " "
			      << lPUJetId.d0()        << " "
			      << lPUJetId.dZ()        << " "
			      << lPUJetId.beta()      << " "
			      << lPUJetId.betaStar()  << " "
			      << lPUJetId.nCharged()  << " "
			      << lPUJetId.nNeutrals() << " "
			      << lPUJetId.dRMean()    << " "
			      << lPUJetId.frac01()    << " "
			      << lPUJetId.frac02()    << " "
			      << lPUJetId.frac03()    << " "
			      << lPUJetId.frac04()    << " "
			      << lPUJetId.frac05()
			      << " === : === "
			      << lPUJetId.mva() << " " << endl;
  }

  return lPUJetId.mva();
}
