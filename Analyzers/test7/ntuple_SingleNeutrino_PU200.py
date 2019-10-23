# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1NtupleMC --step L1 --mc --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --conditions auto:phase2_realistic --beamspot HLLHC14TeV --geometry Extended2023D17 --era Phase2_timing --filein file:dummy.root --no_exec --nThreads 4 -n 100
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
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/1808B0B0-6A5C-E811-9CFF-0025904C66E4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/2485A703-DE5B-E811-BF8B-0CC47AFB7FFC.root',
    ),
    secondaryFileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring(),
    #skipEvents = cms.untracked.uint32(1090),
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
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

## Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

##Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)


# ##############################################################################
# From Vladimir
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_10_1_5

if False:
    process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')

    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

    process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
    process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)

    process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
    process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

    #process.TTClusterAssociatorFromPixelDigis.digiSimLinks          = cms.InputTag( "simSiPixelDigis","Tracker" )
    process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)

    # Path and EndPath definitions
    process.L1simulation_step = cms.Path(process.SimL1Emulator)
    process.endjob_step = cms.EndPath(process.endOfProcess)
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

    # Schedule definition
    #process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
    process.schedule = cms.Schedule(process.EcalEBtp_step,process.hgcl1tpg_step,process.L1simulation_step,process.L1TrackTrigger_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
# ##############################################################################

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# ______________________________________________________________________________
# Modify input files
if True:
  from L1TMuonSimulations.Configuration.tools import *
  txt = 'L1TMuonSimulations/Configuration/data/input_SingleNeutrino_PhaseIIFall17D-L1TPU200.txt'
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
        'drop *_simEmtfDigis_*_*',
        'drop *_simOmtfDigis_*_*',
        'drop *_simGmtStage2Digis_*_*',
        'drop *_simGtStage2Digis_*_*',
    )
    #process.step1 = cms.Path(process.simCscTriggerPrimitiveDigis + process.simEmtfDigis)
    #process.schedule = cms.Schedule(process.step1)

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
if True:
    # Ntuplize
    process.load('L1TMuonSimulations.Analyzers.ntupler_cfi')
    process.TFileService = cms.Service('TFileService', fileName = process.ntupler.outFileName)
    # Modify sequences without any consequences
    #process.SimL1TMuon = cms.Sequence(process.simCscTriggerPrimitiveDigis + process.simEmtfDigis)
    process.SimL1TMuon = cms.Sequence(process.SimL1TMuonCommon + process.rpcRecHits + process.simTwinMuxDigis + process.me0TriggerPseudoDigiSequence + process.simEmtfDigis)
    process.SimL1EmulatorCore = cms.Sequence(process.SimL1TMuon)
    process.SimL1Emulator = cms.Sequence(process.SimL1EmulatorCore)
    process.ntuple_step = cms.Path(process.ntupler)
    process.schedule = cms.Schedule(process.L1simulation_step, process.ntuple_step)


# ______________________________________________________________________________
# Configure framework report and summary
process.options.wantSummary = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open('dump.py', 'w') as f:
    f.write(process.dumpPython())

