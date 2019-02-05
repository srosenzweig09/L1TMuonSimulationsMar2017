# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: L1TMuonSimulations/Configuration/python/SingleMuonFlatOneOverPt2To7000_PositiveEndCap_cfi.py --step GEN,SIM,DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,RAW2DIGI --mc --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --processName L1 --conditions auto:phase2_realistic --beamspot HLLHC14TeV --geometry Extended2023D17 --era Phase2_timing --pileup NoPileUp --customise SimGeneral/MixingModule/customiseStoredTPConfig.higherPtTP,SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --python_filename pset_SingleMuon_PositiveEndCap.py --fileout file:SingleMuon_PositiveEndCap.root --no_exec --nThreads 4 -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2_timing)

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

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:SingleMuon_PositiveEndCap.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
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
        MaxEta = cms.double(1.4),
        MaxPhi = cms.double(3.14159265359),
        MaxPt = cms.double(7000.0),
        MinEta = cms.double(0.8),
        MinPhi = cms.double(-3.14159265359),
        MinPt = cms.double(3.0),
        PartID = cms.vint32(-13),
        PtSpectrum = cms.string('flatOneOverPt'),
        RandomCharge = cms.bool(True)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single muon+/- pt 3 to 7000 flat in 1/pt positive overlap')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
#process.pL1TkPrimaryVertex = cms.Path(process.L1TkPrimaryVertex)
#process.pL1TkElectrons = cms.Path(process.L1TkElectrons)
#process.pL1TkJets = cms.Path(process.L1TkJets)
#process.pL1TkPhotons = cms.Path(process.L1TkPhotons)
#process.pL1TkMuon = cms.Path(process.L1TkMuons)
#process.pL1TrkMET = cms.Path(process.L1TkEtMiss)
#process.pL1TkTauFromCalo = cms.Path(process.L1TkTauFromCalo)
#process.pL1TkHTMissVtx = cms.Path(process.L1TkHTMissVtx)
#process.pL1TkIsoElectrons = cms.Path(process.L1TkIsoElectrons)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

## Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.L1TrackTrigger_step,process.digi2raw_step,process.raw2digi_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

##Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq

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
# Check LCT BX shift
import os
if 'CMSSW_VERSION' not in os.environ:
    raise RunTimeError('Could not determine CMSSW version.')
cmssw_version = os.environ['CMSSW_VERSION']
cmssw_version = cmssw_version[6:].split('_')[:3]
cmssw_version = tuple(int(x) for x in cmssw_version)
if cmssw_version < (10, 2, 0):
    process.simEmtfDigis.CSCInputBXShift = cms.int32(-6)

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
    process.load('L1TMuonSimulations.Analyzers.rpcintegration_cfi')
    process.ntupler.outFileName = 'ntuple_SingleMuon_Overlap.root'
    process.ntupler.verbosity = 0
    process.TFileService = cms.Service('TFileService', fileName = cms.string(process.ntupler.outFileName.value()))
    # Modify sequences without any consequences
    process.doAllDigi = cms.Sequence(process.generatorSmeared+process.muonDigi)
    process.SimL1TMuon = cms.Sequence(process.SimL1TMuonCommon+process.rpcRecHits+process.simTwinMuxDigis+process.me0TriggerPseudoDigiSequence+process.simEmtfDigis+process.simOmtfDigis)
    process.SimL1EmulatorCore = cms.Sequence(process.SimL1TMuon)
    process.SimL1Emulator = cms.Sequence(process.SimL1EmulatorCore)
    process.ntuple_step = cms.Path(process.ntupler)
    process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step, process.ntuple_step)


# ______________________________________________________________________________
# Configure framework report and summary
process.options.wantSummary = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open('dump.py', 'w') as f:
    f.write(process.dumpPython())
