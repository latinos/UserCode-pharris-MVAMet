//--------------------------------------------------------------------------------------------------
// $Id $
//
// Met Regression
//
// Authors: P. Harris
//--------------------------------------------------------------------------------------------------

#ifndef PH_MVAMet_H
#define PH_MVAMet_H

#include <utility>
#include <vector>
#include <TString.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "pharris/MVAMet/interface/GBRForest.h"
#include "pharris/MVAMet/interface/MetUtilities.h"

class MVAMet {
 public:
  typedef math::XYZVector         Vector;
  typedef math::XYZTLorentzVector LorentzVector;

  MVAMet(double iDZCut=0.1);
  ~MVAMet();
  
  enum MVAType {
    kBaseline = 0
  };
  
  void    Initialize(const edm::ParameterSet &iConfig, 
		     TString iU1Weights    ="$CMSSW_BASE/src/MitPhysics/data/gbrmet.root",
		     TString iPhiWeights   ="$CMSSW_BASE/src/MitPhysics/data/gbrmetphi.root",
		     MVAMet::MVAType  iType=kBaseline);
  
  Bool_t   IsInitialized() const { return fIsInitialized; }
  Double_t evaluatePhi();
  Double_t evaluateU1();
  Double_t MVAValue(  bool iPhi,
		      Float_t iPFSumEt, 
		      Float_t iU      ,
		      Float_t iUPhi   ,
		      Float_t iTKSumEt,
		      Float_t iTKU    ,
		      Float_t iTKUPhi ,
		      Float_t iNPSumEt,
		      Float_t iNPU    ,
		      Float_t iNPUPhi ,
		      Float_t iPUSumEt,
		      Float_t iPUMet  ,
		      Float_t iPUMetPhi,
		      Float_t iPCSumEt,
		      Float_t iPCU    ,
		      Float_t iPCUPhi ,
		      Float_t iJSPt1  ,
		      Float_t iJSEta1 ,
		      Float_t iJSPhi1 ,
		      Float_t iJSPt2  ,
		      Float_t iJSEta2 ,
		      Float_t iJSPhi2 ,
		      Float_t iNJet   ,
		      Float_t iNAllJet,
		      Float_t iNPV    );
  
  std::pair<LorentzVector,double>  GetMet(std::vector<LorentzVector>                                       &iVis,
					  std::vector<std::pair<std::pair<LorentzVector,double>,double> >  &iJets,
					  std::vector<std::pair<LorentzVector,double> >                    &iCands,
					  std::vector<Vector>                                              &iVertices,
					  bool iPrintDebug=false);
  
  MetUtilities *fUtils;
    
  protected:
    TString      fPhiMethodName;
    TString      fU1MethodName;
    Bool_t       fIsInitialized;
    MVAType      fType;
    double  fDZCut  ;
    Float_t fU      ;
    Float_t fUPhi   ;
    Float_t fTKSumEt;
    Float_t fTKU    ;
    Float_t fTKUPhi ;
    Float_t fNPSumEt;
    Float_t fNPU    ;
    Float_t fNPUPhi ;
    Float_t fPUSumEt;
    Float_t fPUMet  ;
    Float_t fPUMetPhi ;
    Float_t fPCSumEt;
    Float_t fPCU    ;
    Float_t fPCUPhi ;
    Float_t fJSPt1  ;
    Float_t fJSEta1 ;
    Float_t fJSPhi1 ;
    Float_t fJSPt2  ;
    Float_t fJSEta2 ;
    Float_t fJSPhi2 ;
    Float_t fNJet   ;
    Float_t fNAllJet;
    Float_t fNPV    ;
    Float_t fUPhiMVA;
    
    Float_t* fPhiVals;
    Float_t* fU1Vals;
    
    
    GBRForest *fPhiReader;
    GBRForest *fU1Reader;
};
#endif
