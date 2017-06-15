import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:SingleMuon_PositiveEndCap_RPCInt.0.root'),
    #fileNames = cms.untracked.vstring('file:SingleMuon_PositiveEndCap_RPCInt.1.root'),
    fileNames = cms.untracked.vstring('file:SingleMuon_PositiveEndCap_RPCInt.2.root'),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.options = cms.untracked.PSet()


## Plugin: RPCIntegration
#process.load("L1TMuonSimulations.Analyzers.rpcintegration_cfi")
#process.rpcintegration.outFileName = "histos.root"
#process.rpcintegration.verbosity = 1

# Plugin: NtupleMaker
process.load("L1TMuonSimulations.Analyzers.rpcintegration_cfi")
process.ntupler.outFileName = "ntuple.root"
process.ntupler.verbosity = 1

# TFileService (needed for NtupleMaker)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(process.ntupler.outFileName._value)
)

# Paths
#process.p = cms.Path(process.rpcintegration)
process.p = cms.Path(process.ntupler)

