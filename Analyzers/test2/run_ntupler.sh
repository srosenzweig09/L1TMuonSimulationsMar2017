#!/bin/bash

PIP=RPCInt.py
PIPPO=pippo_RPCInt.py


echo "running SingleNeutrino_PU140"
cp $PIP $PIPPO
cat <<EOF >> $PIPPO
if True:
  txt = 'L1TMuonSimulations/Configuration/data/input_SingleNeutrino_PU140.txt'
  txt = os.path.join(os.environ['CMSSW_BASE'], 'src', txt)
  fileNames_txt = loadFromFile(txt, fmt='root://cmsxrootd.fnal.gov/%s')
  process.source.fileNames = fileNames_txt
  #
  process.load('L1TMuonSimulations.Analyzers.rpcintegration_cfi')
  process.ntupler.outFileName = 'ntuple_SingleNeutrino_PU140.root'
  process.ntupler.docString = 'SingleNeutrino_PU140'
  process.ntupler.verbosity = 0
  process.TFileService = cms.Service('TFileService', fileName = cms.string(process.ntupler.outFileName.value()))
  #
  process.ntuple_step = cms.Path(process.ntupler)
  process.step1 = cms.Path(process.simEmtfDigis)
  process.schedule = cms.Schedule(process.step1, process.ntuple_step)
  #
  process.maxEvents.input = cms.untracked.int32(-1)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000
EOF
cmsRun $PIPPO >pippo_RPCInt.log 2>&1

