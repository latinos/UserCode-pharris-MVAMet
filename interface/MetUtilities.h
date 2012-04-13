#ifndef PH_METUTILITIES_H
#define PH_METUTILITIES_H

#include <utility>
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class MetUtilities {
 public:
  typedef math::XYZVector         Vector;
  typedef math::XYZTLorentzVector LorentzVector;

  MetUtilities();
  virtual ~MetUtilities();

  bool              passMVA  (std::pair<LorentzVector,double> iJet);
  LorentzVector    *leadPt   (std::vector<std::pair<LorentzVector,double> > &iJets,bool iFirst);
  int               NJets    (std::vector<std::pair<LorentzVector,double> > &iJets,double iPt);
  double            deltaR   (LorentzVector &iVec1,LorentzVector &iVec2);
  void              cleanJets(std::vector<LorentzVector> &iVis,std::vector<std::pair<LorentzVector,double> > &iJets);
  std::pair<LorentzVector,double> TKMet   (std::vector<std::pair<LorentzVector,double> > &iCands,int iDZ,bool iLowDz);
  std::pair<LorentzVector,double> JetMet  (std::vector<std::pair<LorentzVector,double> > &iJets ,bool iPassMVA);
  std::pair<LorentzVector,double> NoPUMet (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<std::pair<LorentzVector,double> > &iJets,double iDZ);
  std::pair<LorentzVector,double> PUMet   (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<std::pair<LorentzVector,double> > &iJets,double iDZ);
  std::pair<LorentzVector,double> PUCMet  (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<std::pair<LorentzVector,double> > &iJets,double iDZ);

  std::pair<LorentzVector,double> PFRecoil  (double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ);
  std::pair<LorentzVector,double> TKRecoil  (double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ);
  std::pair<LorentzVector,double> NoPURecoil(double iSumEt,LorentzVector iVis,
					     std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<std::pair<LorentzVector,double> > &iJets,double iDZ);
  std::pair<LorentzVector,double> PUCRecoil(double iSumEt,LorentzVector iVis,
					    std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<std::pair<LorentzVector,double> > &iJets,double iDZ);
 protected:
  // PU jet identifier 
  Float_t fMVACut[3][4][4];  //Jet Id MVA
};

#endif
