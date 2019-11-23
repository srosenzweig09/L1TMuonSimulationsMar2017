import FWCore.ParameterSet.Config as cms

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
    crossingFrame = cms.bool(False),

    outFileName = cms.string('ntuple.root'),
    verbosity = cms.untracked.int32(0),
)