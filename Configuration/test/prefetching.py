import FWCore.ParameterSet.Config as cms

process = cms.Process("Whatever")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# ______________________________________________________________________________
# Modify input files
if True:
  from L1TMuonSimulations.Configuration.tools import *
  txt = 'L1TMuonSimulations/Configuration/data/DisplacedSUSY/input_DisplacedSUSY_SmuonToMuNeutralino_M-500_CTau-1000.txt'
  txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
  fileNames_txt = loadFromFile(txt, fmt='')
  process.source.fileNames = fileNames_txt

# ______________________________________________________________________________
import ROOT
good_list = []
for fname in fileNames_txt:
  redirector = 'root://cmsxrootd.fnal.gov/'
  fname = redirector + fname
  root_file = ROOT.TFile.Open(fname, 'READ')
  if not root_file:
    print('[ERROR] Failed to open the file: %s' % fname)
  else:
    good_list.append(fname)

print('[INFO] The following are good kids:')
for fname in good_list:
  print(fname)
