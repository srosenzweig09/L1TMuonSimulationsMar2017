import FWCore.ParameterSet.Config as cms


rpcintegration = cms.EDAnalyzer('RPCIntegration',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('histos.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_rpcintegration(process):
    process.TFileService.fileName = process.rpcintegration.outFileName.value()

trackcounting = cms.EDAnalyzer('TrackCounting',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    gmtMuonTag = cms.InputTag('simGmtStage2Digis'),
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('histos.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_trackcounting(process):
    process.TFileService.fileName = process.trackcounting.outFileName.value()

ntupler = cms.EDAnalyzer('NtupleMaker',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    tkTrackTag = cms.InputTag('TTTracksFromTracklet', 'Level1TTTracks'),
    tkTrackAssocTag = cms.InputTag('TTTrackAssociatorFromPixelDigis', 'Level1TTTracks'),
    genPartTag = cms.InputTag('genParticles'),
    trkPartTag = cms.InputTag('mix', 'MergedTrackTruth'),
    pileupInfoTag = cms.InputTag('addPileupInfo'),

    cscSimHitsTag = cms.InputTag('g4SimHits','MuonCSCHits'),
    cscSimHitsXFTag = cms.InputTag('mix', 'g4SimHitsMuonCSCHits'),
    cscStripSimLinksTag = cms.InputTag('simMuonCSCDigis', 'MuonCSCStripDigiSimLinks'),
    cscWireSimLinksTag = cms.InputTag('simMuonCSCDigis', 'MuonCSCWireDigiSimLinks'),
    rpcSimHitsTag = cms.InputTag('g4SimHits','MuonRPCHits'),
    rpcSimHitsXFTag = cms.InputTag('mix', 'g4SimHitsMuonRPCHits'),
    rpcDigiSimLinksTag = cms.InputTag('simMuonRPCDigis', 'RPCDigiSimLink'),
    gemSimHitsTag = cms.InputTag('g4SimHits','MuonGEMHits'),
    gemSimHitsXFTag = cms.InputTag('mix', 'g4SimHitsMuonGEMHits'),
    gemDigiSimLinksTag = cms.InputTag('simMuonGEMDigis', 'GEM'),
    gemStripSimLinksTag = cms.InputTag('simMuonGEMDigis', 'GEM'),
    me0SimHitsTag = cms.InputTag('g4SimHits','MuonME0Hits'),
    me0SimHitsXFTag = cms.InputTag('mix', 'g4SimHitsMuonME0Hits'),
    me0DigiSimLinksTag = cms.InputTag('simMuonME0Digis', 'ME0'),
    me0StripSimLinksTag = cms.InputTag('simMuonME0Digis', 'ME0'),
    dtSimHitsTag = cms.InputTag('g4SimHits','MuonDTHits'),
    dtSimHitsXFTag = cms.InputTag('mix', 'g4SimHitsMuonDTHits'),
    dtDigiSimLinksTag = cms.InputTag('simMuonDTDigis'),
    crossingFrame = cms.bool(True),

    outFileName = cms.string('ntuple.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_ntupler(process):
    process.TFileService.fileName = process.ntupler.outFileName.value()

minbiasmuonanalyzer = cms.EDAnalyzer('MinBiasMuonAnalyzer',
    simTrackTag = cms.InputTag('g4SimHits'),
    trkPartTag = cms.InputTag('mix', 'MergedTrackTruth'),
    outFileName = cms.string('histos.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_minbiasmuonanalyzer(process):
    process.TFileService.fileName = process.minbiasmuonanalyzer.outFileName.value()

ttstubntupler = cms.EDAnalyzer('TTStubNtupleMaker',
    ttstubTag = cms.InputTag('TTStubsFromPhase2TrackerDigis', 'StubAccepted'),
    ttstubAssocTag = cms.InputTag('TTStubAssociatorFromPixelDigis', 'StubAccepted'),
    genPartTag = cms.InputTag('genParticles'),
    trkPartTag = cms.InputTag('mix', 'MergedTrackTruth'),
    outFileName = cms.string('ntuple.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_ttstubntupler(process):
    process.TFileService.fileName = process.ttstubntupler.outFileName.value()

