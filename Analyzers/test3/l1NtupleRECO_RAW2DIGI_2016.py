# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1NtupleRECO -s RAW2DIGI --era=Run2_2016 --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --conditions=auto:run2_data -n 1000 --data --no_exec --no_output --filein=/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/9ED0BEE9-2E3F-E711-B04D-A0369F6369D2.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Run2_2016)

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
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/0001016E-AD3F-E711-BB8F-6CC2173BB830.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/000E8ACA-103F-E711-9EC1-001E67E6F503.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/0018E116-F03E-E711-AA00-A0369F836280.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/007D19E2-8B3F-E711-A9AD-001E677926FC.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/00975D8D-C83F-E711-85DD-008CFAF7245E.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/009D720B-DB3E-E711-829B-001E677926B4.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/00CFD3D4-FB3E-E711-9568-F04DA2752F68.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/00DC7AD3-4641-E711-A688-0090FAA57B20.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/00DFD501-943F-E711-8D06-7CD30AD08EBA.root',
        '/store/data/Run2016G/SingleMuon/RAW-RECO/ZMu-18Apr2017-v1/00000/00E361DA-D13E-E711-9FAB-3417EBE644B3.root',
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
    process.load('L1Trigger.L1TMuonEndCap.fakeEmtfParams_2016_data_cff')
    #delattr(process, 'emtfForestsSource')
    #delattr(process, 'emtfForestsDB')
    delattr(process, 'emtfParamsSource')
    delattr(process, 'emtfParams')
if True:
    process.emtfForestsDB = cms.ESSource(
        "EmptyESSource",
        recordName = cms.string('L1TMuonEndCapForestRcd'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
    )
    process.emtfForests = cms.ESProducer(
        "L1TMuonEndCapForestESProducer",
        PtAssignVersion = cms.int32(5),
        bdtXMLDir = cms.string("v_16_02_21")  # corresponding to pT LUT v5
    )
if True:
    process.TFileService.fileName = cms.string('L1Ntuple_2016.root')
    process.L1TReEmul.remove(process.simEcalTriggerPrimitiveDigis)
    process.L1TReEmul.remove(process.simHcalTriggerPrimitiveDigis)
if True:
    process.load('L1TMuonSimulations.Analyzers.emuaccuracy_cfi')
    process.emuaccuracy.verbosity = 1
    process.p = cms.Path(process.emuaccuracy)
    process.schedule.append(process.p)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())

