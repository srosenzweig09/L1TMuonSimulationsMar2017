import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Muon")

process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_v14', '')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.a1 = cms.EDAnalyzer("AnalyzerCSCGeometry",
    csv = cms.string('cscgeometry.csv'),
    verbosity = cms.int32(1),
)

process.p1 = cms.Path(process.a1)
process.schedule = cms.Schedule(process.p1)

