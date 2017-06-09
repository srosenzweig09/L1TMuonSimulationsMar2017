import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.root'),
)

process.options = cms.untracked.PSet()


# Plugin: RPCIntegration
process.load("L1TMuonSimulations.RPCIntegrationMar2017.rpcintegration_cfi")
process.rpcintegration.outFileName = "histos.root"
process.rpcintegration.verbosity = 1

# Plugin: NtupleMaker
#process.load("L1TMuonSimulations.RPCIntegrationMar2017.rpcintegration_cfi")
#process.ntuplemaker.outFileName = "ntuple.root"
#process.ntuplemaker.verbosity = 1

# TFileService (needed for NtupleMaker)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(process.ntuplemaker.outFileName._value)
)

# Paths
process.p = cms.Path(process.rpcintegration)
#process.p = cms.Path(process.ntuplemaker)

