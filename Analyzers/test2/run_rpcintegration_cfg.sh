#!/bin/bash

PIP=rpcintegration_cfg.py
PIPPO=pippo_rpcintegration_cfg.py
#INPUT="process.source.fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.root')"
INPUT=


echo "running csc"
cp $PIP $PIPPO
echo $INPUT >> $PIPPO
echo "process.rpcintegration.verbosity = cms.untracked.int32(0)" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigisCSC')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigisCSC')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.csc.root')" >> $PIPPO
echo "use_fs_rpcintegration(process)" >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

echo "running rpc"
cp $PIP $PIPPO
echo $INPUT >> $PIPPO
echo "process.rpcintegration.verbosity = cms.untracked.int32(0)" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigisRPC')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigisRPC')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.rpc.root')" >> $PIPPO
echo "use_fs_rpcintegration(process)" >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

echo "running gem"
cp $PIP $PIPPO
echo $INPUT >> $PIPPO
echo "process.rpcintegration.verbosity = cms.untracked.int32(0)" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigisGEM')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigisGEM')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.gem.root')" >> $PIPPO
echo "use_fs_rpcintegration(process)" >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

echo "running irpc"
cp $PIP $PIPPO
echo $INPUT >> $PIPPO
echo "process.rpcintegration.verbosity = cms.untracked.int32(0)" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigis')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigis')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.irpc.root')" >> $PIPPO
echo "use_fs_rpcintegration(process)" >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

