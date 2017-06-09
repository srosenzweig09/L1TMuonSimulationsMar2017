#!/bin/bash

PIPPO=pippo_rpcintegration_cfg.py


echo "running csc"
cp rpcintegration_cfg.py $PIPPO
echo "process.source.fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.root')" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigisCSC')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigisCSC')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.csc.root')"   >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

echo "running rpc"
cp rpcintegration_cfg.py $PIPPO
echo "process.source.fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.root')" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigisRPC')" >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigisRPC')" >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.rpc.root')"   >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

echo "running gem"
cp rpcintegration_cfg.py $PIPPO
echo "process.source.fileNames = cms.untracked.vstring('file:l1NtupleMC_RAW2DIGI.full.root')" >> $PIPPO
echo "process.rpcintegration.emuHitTag   = cms.InputTag('simEmtfDigis')"    >> $PIPPO
echo "process.rpcintegration.emuTrackTag = cms.InputTag('simEmtfDigis')"    >> $PIPPO
echo "process.rpcintegration.outFileName = cms.string('histos.gem.root')"   >> $PIPPO
cmsRun $PIPPO >/dev/null 2>&1

