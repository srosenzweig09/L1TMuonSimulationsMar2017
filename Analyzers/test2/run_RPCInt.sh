#!/bin/bash

PIP=RPCInt.py
PIPPO=pippo_RPCInt.py


echo "running simEmtfDigisCSC*simEmtfDigisRPC*simEmtfDigisGEM*simEmtfDigis"
cp $PIP $PIPPO
cat <<EOF >>$PIPPO
# Run quadruple EMTF emulators
process.simEmtfDigisCSC = process.simEmtfDigis.clone(RPCEnable = False, GEMEnable = False, IRPCEnable = False, TTEnable = False)
process.simEmtfDigisRPC = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = False, IRPCEnable = False, TTEnable = False)
process.simEmtfDigisGEM = process.simEmtfDigis.clone(RPCEnable = True , GEMEnable = True , IRPCEnable = False, TTEnable = False)
process.step1.replace(process.simEmtfDigis, process.simEmtfDigisCSC*process.simEmtfDigisRPC*process.simEmtfDigisGEM*process.simEmtfDigis)
process.RAWSIMoutput.outputCommands.append('keep *_simEmtfDigis*_*_*')

process.RAWSIMoutput.fileName = cms.untracked.string(process.RAWSIMoutput.fileName.value().replace('.root','.full.root'))

process.maxEvents.input = cms.untracked.int32(-1)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
EOF
cmsRun $PIPPO >pippo_RPCInt.log 2>&1

