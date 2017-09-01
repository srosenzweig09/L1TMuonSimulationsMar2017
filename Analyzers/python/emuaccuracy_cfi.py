import FWCore.ParameterSet.Config as cms

emuaccuracy = cms.EDAnalyzer('EmuAccuracy',
    unpHitTag = cms.InputTag('emtfStage2Digis'),
    unpTrackTag = cms.InputTag('emtfStage2Digis'),
    #emuHitTag = cms.InputTag('simEmtfDigis', 'EMTF'),
    #emuTrackTag = cms.InputTag('simEmtfDigis', 'EMTF'),
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    emuHitTag2 = cms.InputTag('simEmtfDigis'),
    emuTrackTag2 = cms.InputTag('simEmtfDigis'),
    verbosity = cms.untracked.int32(0),
)
