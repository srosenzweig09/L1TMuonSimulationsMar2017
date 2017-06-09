import FWCore.ParameterSet.Config as cms
import RPCInt as pset; process = pset.process

if True:
  import glob
  dirname = '/store/user/jiafulow/L1MuonTrigger/9_0_0/SingleMuon_PositiveEndCap/ParticleGuns/CRAB3/170510_154046/0000/'
  fileNames3 = glob.glob('/cms/data'+dirname+'*.root')
  fileNames3 = [x[9:] for x in fileNames3]
process.source.fileNames = cms.untracked.vstring(fileNames3)

process.RAWSIMoutput.fileName = cms.untracked.string('l1NtupleMC_RAW2DIGI.full.root')

process.maxEvents.input = cms.untracked.int32(-1)
