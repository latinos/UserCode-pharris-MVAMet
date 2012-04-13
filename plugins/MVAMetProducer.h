// user include files
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "CMGTools/External/interface/PileupJetIdAlgo.h"

#include "pharris/MVAMet/interface/MetUtilities.h"
#include "pharris/MVAMet/interface/MVAMet.h"

using namespace reco;

class MVAMetProducer : public edm::EDProducer {
 public:
  typedef math::XYZTLorentzVector LorentzVector;
  typedef math::XYZVector         Vector;
  explicit MVAMetProducer(const edm::ParameterSet&);
  ~MVAMetProducer();
  
 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag fCorrJetName;     
  edm::InputTag fUnCorrJetName; 
  edm::InputTag fPFCandName; 
  edm::InputTag fVertexName; 

  double                  fJetPtMin;
  double                  fDZMin;
  PileupJetIdAlgo        *fPUJetIdAlgo;
  MVAMet                 *fMVAMet;

  void makeJets      (std::vector<std::pair<std::pair<LorentzVector,double>,double> > &iJetInfo,PFJetCollection     &iUCJets,PFJetCollection &iCJets,VertexCollection &iVertices); 
  void makeCandidates(std::vector<          std::pair<LorentzVector,double> >         &iPFInfo,PFCandidateCollection &iCands,Vertex *iPV);
  void makeVertices  (std::vector<Vector>        &iPVInfo,VertexCollection &iVertices);

  bool   passPFLooseId(const PFJet *iJet);
  double pfCandDz(const PFCandidate* iPFCand, const Vertex *iPV) ;
  double jetMVA  (const PFJet *iuncorrJet,double iJec, const Vertex iPV, const reco::VertexCollection &iAllvtx);
};
