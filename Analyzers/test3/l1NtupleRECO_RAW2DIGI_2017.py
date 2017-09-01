# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1NtupleRECO -s RAW2DIGI --era=Run2_2017 --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --conditions=auto:run2_data -n 1000 --data --no_exec --no_output --filein=/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/678/00000/6E1B76C0-A466-E711-A521-02163E0128F4.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/678/00000/6E1B76C0-A466-E711-A521-02163E0128F4.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/679/00000/9AD7D2B4-A766-E711-B44E-02163E011E6F.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/681/00000/B2370662-AE66-E711-9BC6-02163E01A4E3.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/009F39A0-8A68-E711-806D-02163E01A1BE.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/04870657-8968-E711-96B0-02163E019CB3.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1C9DD9FE-9E68-E711-BA29-02163E01A3B2.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1E2469ED-9768-E711-AE47-02163E01A5AC.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1E44504F-8968-E711-9E26-02163E01A737.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/2EC43DA1-8C68-E711-BDC0-02163E01A6F7.root',
        '/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/30A9DAE4-9268-E711-A52C-02163E0124B2.root',
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1NtupleRECO nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODRAWEMU 

#call to customisation function L1NtupleAODRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleAODRAWEMU(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


if True:
    #process.load('L1Trigger.L1TMuonEndCap.fakeEmtfParams_2017_data_cff')
    #delattr(process, 'emtfForestsSource')
    #delattr(process, 'emtfForestsDB')
    delattr(process, 'emtfParamsSource')
    delattr(process, 'emtfParams')
if True:
    process.TFileService.fileName = cms.string('L1Ntuple_2017.root')
    process.L1TReEmul.remove(process.simEcalTriggerPrimitiveDigis)
    process.L1TReEmul.remove(process.simHcalTriggerPrimitiveDigis)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

