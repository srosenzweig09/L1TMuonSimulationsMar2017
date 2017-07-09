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
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('histos.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_trackcounting(process):
  process.TFileService.fileName = process.trackcounting.outFileName.value()

ntupler = cms.EDAnalyzer('NtupleMaker',
    emuHitTag = cms.InputTag('simEmtfDigis'),
    emuTrackTag = cms.InputTag('simEmtfDigis'),
    genPartTag = cms.InputTag('genParticles'),
    outFileName = cms.string('ntuple.root'),
    verbosity = cms.untracked.int32(0),
)

def use_fs_ntupler(process):
  process.TFileService.fileName = process.ntupler.outFileName.value()

