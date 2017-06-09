import FWCore.ParameterSet.Config as cms

rpcintegration = cms.EDAnalyzer('RPCIntegration',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('histos.root'),
    verbosity = cms.untracked.int32(0),
)

ntuplemaker = cms.EDAnalyzer('NtupleMaker',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('ntuple.root'),
    verbosity = cms.untracked.int32(0),
)
