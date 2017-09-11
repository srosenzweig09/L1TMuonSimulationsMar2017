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
#process.maxEvents.input = cms.untracked.int32(-1)
process.maxEvents.input = cms.untracked.int32(200000)

fileNames = [
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_71.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_72.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_73.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_74.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_75.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_76.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_77.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_78.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_79.root",
    "/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/SingleMuon_PositiveEndCap_80.root",
]
process.source.fileNames = cms.untracked.vstring(fileNames)

if False:
  import glob
  dirname = '/store/user/jiafulow/L1MuonTrigger/P2_9_2_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170613_001230/0000/'
  fileNames_glob = glob.glob('/cms/data'+dirname+'*.root')
  fileNames_glob = [x[9:] for x in fileNames_glob]
  process.source.fileNames = fileNames_glob

if True:
  # Input
  from L1TMuonSimulations.Configuration.tools import *
  txt = 'L1TMuonSimulations/Configuration/data/input_SingleMuon_PositiveEndCap.txt'
  #txt = 'L1TMuonSimulations/Configuration/data/input_SingleNeutrino_PU50.txt'
  #txt = 'L1TMuonSimulations/Configuration/data/input_SingleNeutrino_PU100.txt'
  #txt = 'L1TMuonSimulations/Configuration/data/input_SingleNeutrino_PU140.txt'
  txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
  fileNames_txt = loadFromFile(txt)
  process.source.fileNames = fileNames_txt
  process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
  # Output
  process.RAWSIMoutput.fileName = cms.untracked.string('file:SingleMuon_PositiveEndCap_RPCInt.root')
  #process.RAWSIMoutput.fileName = cms.untracked.string('file:SingleNeutrino_PU50_RPCInt.root')
  #process.RAWSIMoutput.fileName = cms.untracked.string('file:SingleNeutrino_PU100_RPCInt.root')
  #process.RAWSIMoutput.fileName = cms.untracked.string('file:SingleNeutrino_PU140_RPCInt.root')
  process.RAWSIMoutput.outputCommands  = ['drop *']
  process.RAWSIMoutput.outputCommands += ['keep *_genParticles_*_*', 'keep *_mix_MergedTrackTruth_*', 'keep *_simCscTriggerPrimitiveDigis_*_*', 'keep *_simMuonRPCDigis_*_*', 'keep *_simMuonGEMDigis_*_*', 'keep *_simMuonGEMPadDigis_*_*', 'keep *_simMuonGEMPadDigiClusters_*_*', 'keep *_simEmtfDigis*_*_*']

if True:
  process.simEmtfDigisCSC  = process.simEmtfDigis.clone(RPCEnable = False, GEMEnable = False, IRPCEnable = False, ME0Enable = False, TTEnable = False)
  process.simEmtfDigisRPC  = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = False, IRPCEnable = False, ME0Enable = False, TTEnable = False)
  process.simEmtfDigisGEM  = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = True , IRPCEnable = False, ME0Enable = False, TTEnable = False)
  process.simEmtfDigisIRPC = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = True , IRPCEnable = True , ME0Enable = False, TTEnable = False)

if True:
  # From https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTrackingTruth
  #process.load("SimGeneral.MixingModule.digitizers_cfi")
  process.mix.playback = cms.untracked.bool(True)
  process.trackingParticles.simHitCollections = cms.PSet(
    pixel = cms.VInputTag(process.trackingParticles.simHitCollections.pixel)
  )
  process.mix.digitizers = cms.PSet(
    pixel = cms.PSet(process.pixelDigitizer),
    mergedtruth = cms.PSet(process.trackingParticles),
  )
  for a in process.aliases: delattr(process, a)

# My paths and schedule definitions
#process.step1 = cms.Path((process.simCscTriggerPrimitiveDigis) + process.simEmtfDigis)
#process.step1 = cms.Path(process.simEmtfDigis)
process.step1 = cms.Path(process.simEmtfDigisCSC*process.simEmtfDigisRPC*process.simEmtfDigisGEM*process.simEmtfDigisIRPC*process.simEmtfDigis)
process.RAWSIMoutput.SelectEvents.SelectEvents = cms.vstring('step1')
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.digitisation_step = cms.Path(cms.SequencePlaceholder("mix"))  # only needed for SingleMuon 170614 sample
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)        # only needed for SingleMuon 170614 sample
process.schedule = cms.Schedule(process.digitisation_step, process.L1TrackTrigger_step, process.step1, process.RAWSIMoutput_step)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Run in unscheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

