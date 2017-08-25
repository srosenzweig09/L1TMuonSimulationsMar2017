# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: L1TMuonSimulations/Configuration/python/SingleMuonFlatOneOverPt2To7000_PositiveEndCap_cfi.py --step GEN,SIM,DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,RAW2DIGI --mc --eventcontent RAWSIM --datatier GEN-SIM-DIGI-RAW --processName RAWSIM --conditions auto:phase2_realistic --beamspot HLLHC14TeV --geometry Extended2023D17 --era Phase2C2_timing --pileup NoPileUp --python_filename pset_SingleMuon_PositiveEndCap.py --fileout file:SingleMuon_PositiveEndCap.root --no_exec --nThreads 4 -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAWSIM',eras.Phase2C2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('L1TMuonSimulations/Configuration/python/SingleMuonFlatOneOverPt2To7000_PositiveEndCap_cfi.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('file:SingleMuon_PositiveEndCap.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer2",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MaxPt = cms.double(7000.0),
        MinEta = cms.double(1.2),
        MinPhi = cms.double(-3.14159265359),
        MinPt = cms.double(2.0),
        PartID = cms.vint32(-13),
        PtSpectrum = cms.string('flatOneOverPt'),
        RandomCharge = cms.bool(True)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single muon+/- pt 2 to 7000 flat in 1/pt positive endcap')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.L1TrackTrigger_step,process.digi2raw_step,process.raw2digi_step,process.endjob_step,process.RAWSIMoutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# ______________________________________________________________________________
# Modify output
process.RAWSIMoutput.outputCommands += ['keep *_mix_MergedTrackTruth_*', 'keep *_mix_Tracker_*']
process.RAWSIMoutput.outputCommands += ['keep *_genParticles_*_*', 'keep *_simDtTriggerPrimitiveDigis_*_*', 'keep *_simCscTriggerPrimitiveDigis_*_*', 'keep *_simMuonDTDigis_*_*', 'keep *_simMuonCSCDigis_*_*', 'keep *_simMuonRPCDigis_*_*', 'keep *_simMuonGEMDigis_*_*', 'keep *_simMuonGEMPadDigis_*_*', 'keep *_simMuonME0Digis_*_*', 'keep *_simMuonME0ReDigis_*_*', 'keep *_simMuonME0PadDigis_*_*', 'keep *_simEmtfDigis*_*_*', 'keep *_simGmtStage2Digis_*_*', 'keep *_TTClustersFromPhase2TrackerDigis_*_*', 'keep *_TTStubsFromPhase2TrackerDigis_*_*']

# My paths and schedule definitions
print("[INFO] Using GlobalTag: %s" % process.GlobalTag.globaltag.value())
if False:
    from Configuration.StandardSequences.SimL1Emulator_cff import simCscTriggerPrimitiveDigis
    process.simCscTriggerPrimitiveDigis = simCscTriggerPrimitiveDigis
    process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")
    process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    from L1Trigger.L1TMuonEndCap.simEmtfDigis_cfi import simEmtfDigisMC
    process.simEmtfDigis = simEmtfDigisMC
    process.simEmtfDigis.verbosity = cms.untracked.int32(1)
if False:
    process.simEmtfDigis.spPCParams16.ZoneBoundaries = [0,36,54,96,127]
    process.simEmtfDigis.spPCParams16.UseNewZones    = True
    process.simEmtfDigis.RPCEnable                   = True
    process.simEmtfDigis.GEMEnable                   = True
    process.simEmtfDigis.Era                         = cms.string('Phase2C2')
    #process.simEmtfDigis.spPAParams16.PtLUTVersion   = cms.int32(5)
if True:
    # Make ME0 pads
    process.load('RecoLocalMuon.GEMRecHit.me0RecHits_cfi')
    process.load('RecoLocalMuon.GEMSegment.me0Segments_cfi')
    process.fakeSimMuonME0PadDigis = cms.EDProducer("FakeME0PadDigiProducer", InputCollection = cms.InputTag("me0Segments"))
    process.me0DigiRecoSequence = cms.Sequence(process.me0RecHits * process.me0Segments * process.fakeSimMuonME0PadDigis)
    process.muonDigi += process.me0DigiRecoSequence
    process.RAWSIMoutput.outputCommands += ['keep *_me0RecHits_*_*', 'keep *_me0Segments_*_*', 'keep *_fakeSimMuonME0PadDigis_*_*']
if True:
    from L1Trigger.L1TMuonEndCap.customise_Phase2C2 import customise as customise_Phase2C2
    process = customise_Phase2C2(process)
process.step1 = cms.Path((process.simCscTriggerPrimitiveDigis) + process.simEmtfDigis)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

