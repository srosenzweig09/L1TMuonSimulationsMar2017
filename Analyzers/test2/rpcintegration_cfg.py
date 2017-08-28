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
    #fileNames = cms.untracked.vstring('file:SingleMuon_PositiveEndCap_RPCInt.2.root'), # default copad, 170614 samples
    #fileNames = cms.untracked.vstring('file:SingleMuon_PositiveEndCap_RPCInt.3.root'), # no copad, 170614 samples
    #fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.20170822.root'),
    fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.20170827.root'),
    #fileNames = cms.untracked.vstring('file:SingleNeutrino_PU50_RPCInt.0.root'),
    #fileNames = cms.untracked.vstring('file:SingleNeutrino_PU100_RPCInt.0.root'),
    #fileNames = cms.untracked.vstring('file:SingleNeutrino_PU140_RPCInt.0.root'),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.options = cms.untracked.PSet()

# TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("histos.root"),
    closeFileFast = cms.untracked.bool(True),
)

# Load the cfi
from L1TMuonSimulations.Analyzers.rpcintegration_cfi import *
process.load("L1TMuonSimulations.Analyzers.rpcintegration_cfi")

# Plugin: RPCIntegration
process.rpcintegration.outFileName = "histos.root"
process.rpcintegration.verbosity = 1
process.p = cms.Path(process.rpcintegration)
use_fs_rpcintegration(process)

## Plugin: TrackCounting
#process.trackcounting.outFileName = "histos.root"
#process.trackcounting.verbosity = 1
#process.p = cms.Path(process.trackcounting)
#use_fs_trackcounting(process)

## Plugin: NtupleMaker
#process.ntupler.outFileName = "ntuple.root"
#process.ntupler.verbosity = 1
#process.p = cms.Path(process.ntupler)
#use_fs_ntupler(process)

