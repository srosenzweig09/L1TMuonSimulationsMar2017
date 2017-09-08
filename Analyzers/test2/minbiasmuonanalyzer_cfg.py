import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('Demo',eras.Phase2C2_timing)

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:dummy.root'),
)

process.options = cms.untracked.PSet()

# TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("minbiasmuonanalyzer.root"),
    closeFileFast = cms.untracked.bool(True),
)

# Load the input
from L1TMuonSimulations.Configuration.tools import *
#txt = 'L1TMuonSimulations/Configuration/data/input_MinBias.txt'
txt = 'L1TMuonSimulations/Configuration/data/input_MinBias_PhaseIITDRSpring17GS.txt'
txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
fileNames_txt = loadFromFile(txt)
process.source.fileNames = fileNames_txt
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# Load the cfi
process.load('L1TMuonSimulations.Analyzers.rpcintegration_cfi')

# Tracking particles
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
#
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(0)
process.mix.maxBunch = cms.int32(0)
process.mix.digitizers = cms.PSet(
    mergedtruth = cms.PSet(process.trackingParticles),
)
process.mix.mixObjects = cms.PSet()
for a in process.aliases: delattr(process, a)

# My paths and schedule definitions
process.step1 = cms.Path(process.mix * process.minbiasmuonanalyzer)

# Configure framework report and summary
process.options.allowUnscheduled = cms.untracked.bool(True)
process.options.wantSummary = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

