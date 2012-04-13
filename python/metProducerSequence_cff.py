import FWCore.ParameterSet.Config as cms
from CMGTools.External.puJetIDAlgo_cff import *

mvaMet   = cms.EDProducer("MVAMetProducer",
                          CorrJetName     = cms.InputTag("ak5PFJetsL1L2L3"),
                          JetName         = cms.InputTag("ak5PFJets"),
                          PFCandidateName = cms.InputTag("particleFlow"),
                          VertexName      = cms.InputTag("offlinePrimaryVertices"),
                          JetPtMin        = cms.double(1.0),
                          dZMin           = cms.double(0.2),
                          impactParTkThreshold = cms.untracked.double(1.) ,
                          tmvaWeights = cms.untracked.string("CMGTools/External/data/mva_JetID_v1.weights.xml"),
                          tmvaMethod  = cms.untracked.string("JetID"),
                          version = cms.untracked.int32(-1),
                          tmvaVariables = cms.untracked.vstring(
    "nvtx",
    "jetPt",
    "jetEta",
    "jetPhi",
    "dZ",
    "d0",
    "beta",
    "betaStar",
    "nCharged",
    "nNeutrals",
    "dRMean",
    "frac01",
    "frac02",
    "frac03",
    "frac04",
    "frac05",
    ),
                          JetIdParams = JetIdParams
                          )

                          

metSequence  = cms.Sequence ( mvaMet )
