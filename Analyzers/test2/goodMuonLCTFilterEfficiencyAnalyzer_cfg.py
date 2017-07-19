# From GeneratorInterface/Core/test/runGenFilterEfficiencyAnalyzer_cfg.py
import FWCore.ParameterSet.Config as cms

process = cms.Process("GenFilterEfficiency")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:MinBias_L1TMuon.0.root')
)

process.dummy = cms.EDAnalyzer("GenFilterEfficiencyAnalyzer",
                               genFilterInfoTag = cms.InputTag("goodMuonLCTFilterEfficiencyProducer")
)

process.p = cms.Path(process.dummy)
