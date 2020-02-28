# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: L1TMuonSimulations/Configuration/python/SingleMuonFlatOneOverPt2To7000_PositiveEndCap_cfi.py --step GEN,SIM,DIGI:pdigi_valid,L1,L1TrackTrigger --conditions auto:phase2_realistic --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --era Phase2C8_timing_layer_bar --beamspot HLLHC14TeV --geometry Extended2023D41 --customise SimGeneral/MixingModule/customiseStoredTPConfig.higherPtTP,SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --python_filename pset_SingleMuon_Endcap_2GeV_PhaseIITDRSpring19.py --fileout file:SingleMuon_Endcap.root --mc --processName L1 --no_exec --nThreads 4 -n 100
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
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
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

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:SingleMuon_Slow.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

def get_fragment():
  # Based on L1T-PhaseIITDRSpring19GS-00167-fragment.py
  FLAVOR = 'stau'
  COM_ENERGY = 14000.
  MASS_POINT = 871   # GeV
  PROCESS_FILE = 'SimG4Core/CustomPhysics/data/RhadronProcessList.txt'
  PARTICLE_FILE = 'Configuration/Generator/data/particles_%s_%d_GeV.txt' % (FLAVOR, MASS_POINT)
  SLHA_FILE = 'Configuration/Generator/data/isa-slha%dGeV.out' % MASS_POINT
  PDT_FILE = 'Configuration/Generator/data/hscppythiapdt%s%d.tbl'  % (FLAVOR, MASS_POINT)
  USE_REGGE = False


  #import FWCore.ParameterSet.Config as cms

  #source = cms.Source("EmptySource")

  #from Configuration.Generator.Pythia8CommonSettings_cfi import *
  #from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

  generator = cms.EDFilter("Pythia8GeneratorFilter",
      pythiaPylistVerbosity = cms.untracked.int32(0),
      filterEfficiency = cms.untracked.double(1),
      pythiaHepMCVerbosity = cms.untracked.bool(False),
      comEnergy = cms.double(COM_ENERGY),
      crossSection = cms.untracked.double(-1),
      maxEventsToPrint = cms.untracked.int32(0),
      SLHAFileForPythia8 = cms.string('%s' % SLHA_FILE),
      PythiaParameters = cms.PSet(
          pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
          processParameters = cms.vstring(
              'SUSY:all = on',
              'SUSY:idA = 1000015',
              'SUSY:idB = 1000015',
              '1000015:mayDecay = off'
               ),
          parameterSets = cms.vstring(
              'pythia8CommonSettings',
              'pythia8CP5Settings',
              'processParameters'
          )
      )
   )

  generator.hscpFlavor = cms.untracked.string(FLAVOR)
  generator.massPoint = cms.untracked.int32(MASS_POINT)
  generator.slhaFile = cms.untracked.string(SLHA_FILE)
  generator.processFile = cms.untracked.string(PROCESS_FILE)
  generator.particleFile = cms.untracked.string(PARTICLE_FILE)
  generator.pdtFile = cms.FileInPath(PDT_FILE)
  generator.useregge = cms.bool(USE_REGGE)

  #ProductionFilterSequence = cms.Sequence(generator)
  return generator

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
process.generator = get_fragment()

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

## Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.L1TrackTrigger_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

##Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path).insert(0, process.generator)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.customiseStoredTPConfig
from SimGeneral.MixingModule.customiseStoredTPConfig import higherPtTP

#call to customisation function higherPtTP imported from SimGeneral.MixingModule.customiseStoredTPConfig
process = higherPtTP(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# ______________________________________________________________________________
# Modify EMTF
if True:
    from L1Trigger.L1TMuonEndCap.customise_Phase2 import customise as customise_Phase2
    process = customise_Phase2(process)

# ______________________________________________________________________________
# Modify paths and schedule definitions
print("[INFO] Using GlobalTag: %s" % process.GlobalTag.globaltag.value())
print("[INFO] Using random number seed: %d" % process.RandomNumberGeneratorService.generator.initialSeed.value())
if True:
    # Ntuplize
    process.load('L1TMuonSimulations.Analyzers.ntupler_cfi')
    process.TFileService = cms.Service('TFileService', fileName = process.ntupler.outFileName)
    # Modify sequences without any consequences
    process.doAllDigiTask = cms.Task(process.generatorSmeared, process.muonDigiTask)
    process.SimL1TMuonTask = cms.Task(process.SimL1TMuonCommonTask, process.me0TriggerPseudoDigiTask, process.me0TriggerPseudoDigiTask105X, process.rpcRecHits, process.simBmtfDigis, process.simEmtfDigis, process.simOmtfDigis, process.simTwinMuxDigis)
    process.SimL1EmulatorCoreTask = cms.Task(process.SimL1TMuonTask)
    process.SimL1EmulatorTask = cms.Task(process.SimL1EmulatorCoreTask)
    process.ntuple_step = cms.Path(process.ntupler)
    process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step, process.ntuple_step)


# ______________________________________________________________________________
# Configure framework report and summary
process.options.wantSummary = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open('dump.py', 'w') as f:
    f.write(process.dumpPython())
