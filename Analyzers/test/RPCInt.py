# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1NtupleMC --step RAW2DIGI --mc --eventcontent RAWSIM --nThreads 4 --conditions 90X_upgrade2023_realistic_v9 --beamspot HLLHC14TeV --geometry Extended2023D4 --era Phase2C2_timing --filein file:dummy.root -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Phase2C2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:dummy.root'),
    secondaryFileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring('keep *', 'drop *_simEmtfDigis_*_*', 'drop *_simGmtStage2Digis_*_*', 'drop *_simGtStage2Digis_*_*'),
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

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('l1NtupleMC_RAW2DIGI.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v9', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.RAWSIMoutput_step)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



# ______________________________________________________________________________
# Modify output
process.RAWSIMoutput.outputCommands = ['drop *', 'keep *_genParticles_*_*', 'keep *_simCscTriggerPrimitiveDigis_*_*', 'keep *_simMuonRPCDigis_*_*', 'keep *_simMuonGEMDigis_*_*', 'keep *_simEmtfDigis_*_*']
#process.RAWSIMoutput.outputCommands.append('keep *_mix_MergedTrackTruth_*')
#for x in ['keep *_simMuonGEMDigis_*_*', 'keep *_simMuonGEMPadDigis_*_*', 'keep *_simMuonGEMPadDigiClusters_*_*']:
#    process.RAWSIMoutput.outputCommands.append(x)

# Modify source
fileNames = [
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/042D2403-4326-E711-94EC-5065F37D1132.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/063545D6-5326-E711-B3D7-5065F37D8152.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/063B7709-4326-E711-A2AF-24BE05C68681.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/06FEAF34-E025-E711-BF53-0CC47A0AD476.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/0833B49E-4926-E711-BEA0-5065F37D90C2.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/0AC074AE-4926-E711-8712-4C79BA3203F5.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/0E743386-E725-E711-A710-0025907859B8.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/10C4732E-0226-E711-9926-002590D9D8BA.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/1493CC7C-E725-E711-A323-003048CB7B30.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/14F582CF-5326-E711-836D-5065F37D21E2.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/16619DE4-4226-E711-B865-E0071B740D80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/1A98CBAF-FA25-E711-BDBA-0025907D2212.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/1CC4AEF9-5326-E711-95B8-A0000420FE80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/1E927703-4326-E711-9EC9-A0000420FE80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/20814EC1-1C26-E711-B6BA-0CC47A0AD6C4.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/2832FF6C-4D26-E711-A9D6-A0000420FE80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/289F5A13-4326-E711-9FF4-5065F3811272.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/2E8EADB5-4926-E711-BC98-A0000620FE80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/2EC94971-4D26-E711-B155-A0000420FE80.root",
    "/store/mc/PhaseIISpring17D/SingleMu_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/340AE58E-4D26-E711-A096-A0000420FE80.root",
]
fileNames2 = [
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_101.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_102.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_103.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_104.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_105.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_106.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_107.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_108.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_109.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_110.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_111.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_112.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_113.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_114.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_115.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_116.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_117.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_118.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_119.root",
    "/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/SingleMuon_PositiveEndCap_120.root",
]
fileNames3 = []
if False:
  import glob
  dirname = '/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/'
  fileNames3 = glob.glob('/cms/data'+dirname+'*.root')
  fileNames3 = [x[9:] for x in fileNames3]
process.source.fileNames = cms.untracked.vstring(fileNames2)


# My paths and schedule definitions
if False:
    process.simMuonRPCDigis.doBkgNoise               = False
    process.simMuonGEMDigis.doBkgNoise               = False
if True:
    from Configuration.StandardSequences.SimL1Emulator_cff import simCscTriggerPrimitiveDigis
    process.simCscTriggerPrimitiveDigis = simCscTriggerPrimitiveDigis
    process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")
    process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    from L1Trigger.L1TMuonEndCap.simEmtfDigis_cfi import simEmtfDigisMC
    process.simEmtfDigis = simEmtfDigisMC
    #process.simEmtfDigis.verbosity = cms.untracked.int32(1)
if True:
    process.load('L1Trigger.L1TMuonEndCap.fakeEmtfParams_cff')
    process.simEmtfDigis.spPCParams16.ZoneBoundaries = [0,36,54,96,127]
    process.simEmtfDigis.spPCParams16.UseNewZones    = True
    process.simEmtfDigis.spPCParams16.FixME11Edges   = True
    #process.simEmtfDigis.spPCParams16.CoordLUTDir    = 'ph_lut_v2'
    process.simEmtfDigis.GEMEnable                   = True
process.step1 = cms.Path((process.simCscTriggerPrimitiveDigis) + process.simEmtfDigis)
process.schedule = cms.Schedule(process.step1, process.RAWSIMoutput_step)

# Run triple EMTF emulators
process.simEmtfDigisCSC = process.simEmtfDigis.clone(RPCEnable = False, GEMEnable = False)
process.simEmtfDigisRPC = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = False)
process.step1.replace(process.simEmtfDigis, process.simEmtfDigisCSC*process.simEmtfDigisRPC*process.simEmtfDigis)
process.RAWSIMoutput.outputCommands.append('keep *_simEmtfDigis*_*_*')


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

