import FWCore.ParameterSet.Config as cms
import pset_SingleMuon_PositiveEndCap as pset
process = pset.process
process.setName_('RAWSIM2')


# ______________________________________________________________________________
# Modify source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:dummy.root'),
    secondaryFileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring('keep *', 'drop *_simEmtfDigis_*_*', 'drop *_simGmtStage2Digis_*_*', 'drop *_simGtStage2Digis_*_*'),
)
process.maxEvents.input = cms.untracked.int32(-1)

fileNames = [
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_1.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_2.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_3.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_4.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_5.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_6.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_7.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_8.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_9.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/SingleMuon_PositiveEndCap_10.root",
]
process.source.fileNames = cms.untracked.vstring(fileNames)
if False:
  import glob
  dirname = '/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170610_075611/0000/'
  fileNames_glob = glob.glob('/cms/data'+dirname+'*.root')
  fileNames_glob = [x[9:] for x in fileNames_glob]
  process.source.fileNames = fileNames_glob
if True:
  process.RAWSIMoutput.fileName = cms.untracked.string('file:SingleMuon_PositiveEndCap_RPCInt.root')
  process.RAWSIMoutput.outputCommands  = ['drop *']
  process.RAWSIMoutput.outputCommands += ['keep *_genParticles_*_*', 'keep *_simCscTriggerPrimitiveDigis_*_*', 'keep *_simMuonRPCDigis_*_*', 'keep *_simMuonGEMDigis_*_*', 'keep *_simMuonGEMPadDigis_*_*', 'keep *_simMuonGEMPadDigiClusters_*_*', 'keep *_simEmtfDigis_*_*']

# My paths and schedule definitions
process.step1 = cms.Path((process.simCscTriggerPrimitiveDigis) + process.simEmtfDigis)
process.RAWSIMoutput.SelectEvents.SelectEvents = cms.vstring('step1')
process.schedule = cms.Schedule(process.step1, process.RAWSIMoutput_step)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

