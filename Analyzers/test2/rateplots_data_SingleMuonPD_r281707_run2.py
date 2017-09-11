# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple --step RAW2DIGI --data --eventcontent RAW --era Run2_2016 --conditions auto:run2_data --customise L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --filein /store/data/Run2016H/SingleMuon/RAW/v1/000/281/707/00000/14FEA8EB-B984-E611-9269-02163E012B55.root --no_exec -n 100
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
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2016H/SingleMuon/RAW/v1/000/281/707/00000/14FEA8EB-B984-E611-9269-02163E012B55.root'),
    secondaryFileNames = cms.untracked.vstring(),
    lumisToProcess = cms.untracked.VLuminosityBlockRange(),
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('l1Ntuple_RAW2DIGI.root'),
    outputCommands = process.RAWEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWoutput_step = cms.EndPath(process.RAWoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.RAWoutput_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# customisation of the process.

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
    from L1TMuonSimulations.Configuration.tools import *
    txt = 'L1TMuonSimulations/Configuration/data/input_SingleMuon_Run2016H_r281707.txt'
    txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
    fileNames_txt = loadFromFile(txt)
    process.source.fileNames = fileNames_txt
    process.source.lumisToProcess = ['281707:99-281707:982', '281707:1000-281707:1065']
if True:
    from L1Trigger.L1TMuonEndCap.customise_Phase2C2 import customise as customise_Phase2C2
    process = customise_Phase2C2(process)
if True:
    #process.schedule.remove(process.L1TrackTrigger_step)
    process.simEmtfDigis.GEMEnable                   = False
    process.simEmtfDigis.IRPCEnable                  = False
    process.simEmtfDigis.TTEnable                    = False
if True:
    process.TFileService = cms.Service("TFileService",
        fileName = cms.string("histos.root"),
        closeFileFast = cms.untracked.bool(True),
    )
    from L1TMuonSimulations.Analyzers.rpcintegration_cfi import *
    process.load("L1TMuonSimulations.Analyzers.rpcintegration_cfi")
    process.trackcounting.outFileName = "rateplots_data.root"
    process.trackcounting.verbosity = 1
    process.p = cms.Path(process.simEmtfDigis * process.trackcounting)
    use_fs_trackcounting(process)

process.schedule = cms.Schedule(process.raw2digi_step, process.L1TReEmulPath, process.p)
#process.maxEvents.input = -1


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Run in unscheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())
