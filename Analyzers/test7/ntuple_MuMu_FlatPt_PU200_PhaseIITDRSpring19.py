# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --step L1 --mc --eventcontent FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --datatier GEN-SIM-DIGI-RAW --conditions auto:phase2_realistic --geometry Extended2023D41 --era Phase2C8_timing_layer_bar --filein file:step0.root --fileout file:step1.root --no_exec --nThreads 4 -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C8_timing_layer_bar_cff import Phase2C8_timing_layer_bar

process = cms.Process('L1',Phase2C8_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIITDRSpring19DR/Nu_E10-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v3/70000/C433BE7F-B206-A948-BAF8-05EEE3F2C4F7.root',
    ),
    secondaryFileNames = cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1NtupleMC nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('l1NtupleMC_L1.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

## Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

##Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# ______________________________________________________________________________
# Modify input files
if True:
  from L1TMuonSimulations.Configuration.tools import *
  txt = 'L1TMuonSimulations/Configuration/data/input_Mu_FlatPt2to100-pythia8-gun_PhaseIITDRSpring19_PU200.txt'
  txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
  fileNames_txt = loadFromFile(txt, fmt='')
  process.source.fileNames = fileNames_txt

# ______________________________________________________________________________
# Drop obsolete input branches
if True:
    process.source.inputCommands = cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016s_*_*_*',
        'drop l1tEMTFHit2016Extras_*_*_*',
        'drop l1tEMTFTrack2016s_*_*_*',
        'drop l1tEMTFTrack2016Extras_*_*_*',
        'drop *_simCscTriggerPrimitiveDigis_*_*',
        'drop *_simBmtfDigis_*_*',
        'drop *_simEmtfDigis_*_*',
        'drop *_simOmtfDigis_*_*',
        'drop *_simGmtStage2Digis_*_*',
        'drop *_simGtStage2Digis_*_*',
    )

# ______________________________________________________________________________
# Modify EMTF
if True:
    from L1Trigger.L1TMuonEndCap.customise_Phase2 import customise as customise_Phase2
    process = customise_Phase2(process)

# ______________________________________________________________________________
# Modify paths and schedule definitions
print("[INFO] Using GlobalTag: %s" % process.GlobalTag.globaltag.value())
if True:
    # Ntuplize
    process.load('L1TMuonSimulations.Analyzers.rpcintegration_cfi')
    process.ntupler.outFileName = 'ntuple_MuMu_FlatPt_PU200.root'
    process.ntupler.verbosity = 0
    process.TFileService = cms.Service('TFileService', fileName = cms.string(process.ntupler.outFileName.value()))
    # Modify sequences without any consequences
    #process.doAllDigiTask = cms.Task(process.generatorSmeared, process.muonDigiTask)
    process.SimL1TMuonTask = cms.Task(process.SimL1TMuonCommonTask, process.me0TriggerPseudoDigiTask, process.me0TriggerPseudoDigiTask105X, process.rpcRecHits, process.simBmtfDigis, process.simEmtfDigis, process.simOmtfDigis, process.simTwinMuxDigis)
    process.SimL1EmulatorCoreTask = cms.Task(process.SimL1TMuonTask)
    process.SimL1EmulatorTask = cms.Task(process.SimL1EmulatorCoreTask)
    process.ntuple_step = cms.Path(process.ntupler)
    process.schedule = cms.Schedule(process.L1simulation_step, process.ntuple_step)


# ______________________________________________________________________________
# Configure framework report and summary
process.options.wantSummary = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open('dump.py', 'w') as f:
    f.write(process.dumpPython())
