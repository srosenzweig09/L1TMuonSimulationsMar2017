# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: L1TMuonSimulations/Configuration/python/PPD-PhaseIITDRSpring17GS-00003-fragment.py --step GEN,SIM --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions auto:phase2_realistic --beamspot HLLHC14TeV --geometry Extended2023D17 --era Phase2C2_timing --python_filename pset_MinBias_TuneCUETP8M1_14TeV-pythia8.py --fileout file:MinBias.root --no_exec --nThreads 4 -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C2_timing)

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
    annotation = cms.untracked.string('L1TMuonSimulations/Configuration/python/PPD-PhaseIITDRSpring17GS-00003-fragment.py nevts:100'),
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
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('file:MinBias.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CUEP8M1Settings', 
            'processParameters'),
        processParameters = cms.vstring('SoftQCD:nonDiffractive = on', 
            'SoftQCD:singleDiffractive = on', 
            'SoftQCD:doubleDiffractive = on'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on')
    ),
    comEnergy = cms.double(14000.0),
    crossSection = cms.untracked.double(71390000000.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


if True:
  ## CMSSW stuff
  process.load('Configuration.StandardSequences.Digi_cff')
  process.load('Configuration.StandardSequences.SimL1Emulator_cff')
  delattr(process, "emtfForestsDB")
  process.mix.digitizers = cms.PSet()
  process.mix.mixObjects = cms.PSet(
      mixSH = cms.PSet(
          crossingFrames = cms.untracked.vstring('MuonCSCHits', 'MuonDTHits', 'MuonRPCHits', 'MuonGEMHits', 'MuonME0Hits'),
          input = cms.VInputTag(cms.InputTag("g4SimHits","MuonCSCHits"), cms.InputTag("g4SimHits","MuonDTHits"), cms.InputTag("g4SimHits","MuonRPCHits"), cms.InputTag("g4SimHits","MuonGEMHits"), cms.InputTag("g4SimHits","MuonME0Hits")),
          subdets = cms.vstring('MuonCSCHits', 'MuonDTHits', 'MuonRPCHits', 'MuonGEMHits', 'MuonME0Hits'),
          type = cms.string('PSimHit'),
      )
  )
  for a in process.aliases: delattr(process, a)
  process.pdigi_valid = cms.Sequence(process.generatorSmeared+process.fixGenInfo+cms.SequencePlaceholder("randomEngineStateProducer")+cms.SequencePlaceholder("mix")+process.muonDigi)
  process.digitisation_step = cms.Path(process.pdigi_valid)
  process.L1simulation_step = cms.Path(process.SimL1Emulator)

  ## My stuff
  # If you see the following error, just check out DataFormats/CSCDigi and add the dictionary
  # > No data dictionary found for the following classes:
  # >
  # >   edm::Wrapper<std::vector<CSCCorrelatedLCTDigi> >
  process.goodMuonLCTs = cms.EDFilter("MuonLCTDigiSelector",
      src = cms.InputTag('simCscTriggerPrimitiveDigis','MPCSORTED'),
      cut = cms.string(''),
      filter = cms.bool(True),
  )

  process.goodMuonLCTFilterEfficiencyProducer = cms.EDProducer('GenFilterEfficiencyProducer',
      filterPath = cms.string('L1TMuon_step')
  )

  process.RAWSIMoutput.SelectEvents.SelectEvents = cms.vstring('L1TMuon_step')
  process.RAWSIMoutput.fileName = cms.untracked.string('file:MinBias_L1TMuon.root')
  process.RAWSIMoutput.outputCommands.extend(cms.untracked.vstring(
      'drop *_mix_*_*',
      'drop *_addPileupInfo_*_*',
      'drop *_*Digis_*_*',
  ))

  process.L1TMuon_step = cms.Path(process.simCscTriggerPrimitiveDigis*process.goodMuonLCTs)
  process.mixfiltersummary_step = cms.EndPath(process.goodMuonLCTFilterEfficiencyProducer)
  process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1TMuon_step,process.mixfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

