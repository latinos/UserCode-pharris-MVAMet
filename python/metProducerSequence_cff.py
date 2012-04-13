import FWCore.ParameterSet.Config as cms
from CMGTools.External.JetIdParams_cfi import *
from CMGTools.External.puJetIDAlgo_cff import PhilV1

mvaMet   = cms.EDProducer("MVAMetProducer",
                          CorrJetName     = cms.InputTag("ak5PFJetsL2L3"),
                          JetName         = cms.InputTag("ak5PFJets"),
                          PFCandidateName = cms.InputTag("particleFlow"),
                          VertexName      = cms.InputTag("offlinePrimaryVertices"),
                          JetPtMin        = cms.double(1.),
                          dZMin           = cms.double(0.1),
                          puJetIDAlgo     = PhilV1,
                          JetIdParams     = JetIdParams
                          )

                          

metSequence  = cms.Sequence ( mvaMet )
