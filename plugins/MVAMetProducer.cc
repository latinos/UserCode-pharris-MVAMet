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
  fCorrJetName    = iConfig.getParameter<edm::InputTag>("JetName");
  fUnCorrJetName  = iConfig.getParameter<edm::InputTag>("CorrJetName");
  fPFCandName     = iConfig.getParameter<edm::InputTag>("PFCandidateName");
  fVertexName     = iConfig.getParameter<edm::InputTag>("VertexName");
  fJetPtMin       = iConfig.getParameter<double>       ("JetPtMin");
  fDZMin          = iConfig.getParameter<double>       ("dZMin");
  fPUJetIdAlgo    = new PileupJetIdAlgo(iConfig);
  fMVAMet         = new MVAMet(fDZMin);
  fMVAMet         ->Initialize(iConfig,
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmet.root"))),
			       TString((getenv("CMSSW_BASE")+string("/src/pharris/MVAMet/data/gbrmetphi.root")))
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
  VertexCollection::const_iterator lIter = lHVertices->begin();
  Vertex lPV; if (lHVertices->size()>0) lPV = *lIter;
  VertexCollection lVertices = *lHVertices;

  //Make Generic Objects
  std::vector<LorentzVector >                   lVisible;
  std::vector<std::pair<LorentzVector,double> > lPFInfo;  makeCandidates(lPFInfo, lCands,&lPV);
  std::vector<std::pair<LorentzVector,double> > lJetInfo; makeJets      (lJetInfo,lUCJets,lCJets,lVertices);
  std::vector<Vector>                           lVtxInfo; makeVertices  (lVtxInfo,lVertices);
  
  //Dummy visible stuff
  double lDummyPt0  = 10; 
  double lDummyPhi0 = 0;
  double lDummyEta0 = 1;
  double lDummyPt1  = 10;
  double lDummyPhi1 = 2.; 
  double lDummyEta1 = -1;
  TLorentzVector lD0; lD0.SetPtEtaPhiM(lDummyPt0,lDummyPhi0,lDummyEta0,0.);
  TLorentzVector lD1; lD1.SetPtEtaPhiM(lDummyPt1,lDummyPhi1,lDummyEta1,0.);

  LorentzVector lVis0; lVis0.SetCoordinates(lD0.Px(),lD0.Py(),lD0.Pz(),lD0.E());
  LorentzVector lVis1; lVis1.SetCoordinates(lD1.Px(),lD1.Py(),lD1.Pz(),lD1.E());
  lVisible.push_back(lVis0);
  lVisible.push_back(lVis1);
  //Calculate the MVA
  std::pair<LorentzVector,double> lMVAMetInfo = fMVAMet->GetMet(lVisible,lPFInfo,lJetInfo,lVtxInfo);

  //Add a PF Met object and put it in the event
  PFMET lDummy;
  PFMET lMVAMet(lDummy.getSpecific(),lMVAMetInfo.second,lMVAMetInfo.first,lPV.position());
  std::auto_ptr<reco::PFMETCollection> pfmetcoll;
  pfmetcoll.reset     ( new reco::PFMETCollection);
  pfmetcoll->push_back( lMVAMet);
  iEvent.put( pfmetcoll );
}
void MVAMetProducer::makeJets(std::vector<std::pair<LorentzVector,double> > &iJetInfo,PFJetCollection &iUCJets,PFJetCollection &iCJets,VertexCollection &iVertices) { 
  for(int i0   = 0; i0 < (int) iUCJets.size(); i0++) {   // uncorrecte jets collection                                                                                                        
    const PFJet       *pUCJet = &(iUCJets.at(i0));
    for(int i1 = 0; i1 < (int) iCJets .size(); i1++) {   // corrected jets collection                                                                                                   
      const PFJet     *pCJet  = &(iCJets.at(i1));
      if(       pUCJet->jetArea() != pCJet->jetArea()                  ) continue;
      if( fabs(pUCJet->eta() - pCJet->eta())         < 0.01            ) continue;
      if( pCJet->pt()                                < fJetPtMin       ) continue;
      if( !passPFLooseId(pCJet)                                        ) continue;
      double lJec = pCJet ->pt()/pUCJet->pt();
      double lMVA = jetMVA(pUCJet,lJec,iVertices.at(0),iVertices);
      double lNeuFrac = (pCJet->neutralEmEnergy()/pCJet->energy() + pCJet->neutralHadronEnergy()/pCJet->energy());
      LorentzVector pVec; pVec.SetCoordinates(pCJet->pt()*lNeuFrac,pCJet->eta(),pCJet->phi(),pCJet->mass());
      std::pair<LorentzVector,double>  pJetObject(pVec,lMVA); 
      iJetInfo.push_back(pJetObject);
    }
  }
}
void MVAMetProducer::makeCandidates(std::vector<std::pair<LorentzVector,double> > &iPFInfo,PFCandidateCollection &iCands,Vertex *iPV) { 
  for(int i0 = 0; i0 < (int)iCands.size(); i0++) {
    const PFCandidate*  pflowCand = &(iCands.at(i0));
    double pDZ = -999;
    if(iPV == 0) pDZ  = pfCandDz(pflowCand,iPV); //If there is no track return negative number -999
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
  if(iJet->energy()== 0)                                  return false;
  if(iJet->neutralHadronEnergy()/iJet->energy() > 0.99)   return false;
  if(iJet->neutralEmEnergy()/iJet->energy()     > 0.99)   return false;
  if(iJet->nConstituents() <  2)                          return false;
  if(iJet->chargedHadronEnergy()/iJet->energy() <= 0 && fabs(iJet->eta()) < 2.4 ) return false;
  if(iJet->chargedEmEnergy()/iJet->energy() >  0.99  && fabs(iJet->eta()) < 2.4 ) return false;
  if(iJet->chargedMultiplicity()            < 1      && fabs(iJet->eta()) < 2.4 ) return false;
  return true;
}
float MVAMetProducer::pfCandDz(const PFCandidate* iPFCand, const Vertex *iPV) { 
  float lDz = -999;
  if(iPFCand->trackRef().isNonnull()) 
    lDz = fabs(iPFCand->trackRef()->dz(iPV->position()));
  else if(iPFCand->gsfTrackRef().isNonnull())
    lDz = fabs(iPFCand->gsfTrackRef()->dz(iPV->position()));
  //else if(iPFCand->muonRef().isNonnull() && iPFCand->muonRef()->innerTrack().isNonnull())
  //lDz = fabs(iPFCand->muonRef()->innerTrack()->dz(iPV->position()));
  return lDz;
}
double MVAMetProducer::jetMVA (const PFJet *iuncorrJet,double iJec, const Vertex iPV, const reco::VertexCollection &iAllvtx) { 
  PileupJetIdentifier lPUJetId =  fPUJetIdAlgo->computeIdVariables(iuncorrJet,iJec,&iPV,iAllvtx,true);
  return lPUJetId.mva();
}
