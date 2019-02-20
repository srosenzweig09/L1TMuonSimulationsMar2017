# L1TMuonSimulationsMar2017

Package to do the Phase 2 studies for Endcap Muon Track Finder (EMTF).

[![Build Status](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017.svg)](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017)
[![CMSSW version](https://img.shields.io/badge/cmssw-CMSSW__10__4__0-002963.svg)](https://github.com/cms-sw/cmssw)
[![Latest tag](https://img.shields.io/github/tag/jiafulow/L1TMuonSimulationsMar2017.svg)](https://github.com/jiafulow/L1TMuonSimulationsMar2017)

## Build

```shell
export SCRAM_ARCH=slc6_amd64_gcc700
scram p -n P2_CMSSW_10_4_0 CMSSW_10_4_0
cd P2_CMSSW_10_4_0/src
cmsenv

# Checkout the emulator in the 'phase2-develop' branch
git clone -b phase2-develop git@github.com:jiafulow/DataFormatsSep2016.git DataFormats
git clone -b phase2-develop git@github.com:jiafulow/L1TriggerSep2016.git L1Trigger
# Checkout this repository
git clone git@github.com:jiafulow/L1TMuonSimulationsMar2017 L1TMuonSimulations
# Compile
scram b -j 10
```

## Version

- v1.2.0 (2019-02-19): Move to CMSSW_10_4_0. Extend to overlap region.

- v1.1.0 (2018-11-06): First results. Using CMSSW_10_1_5. Include GEM, ME0, iRPC. Add CSC bend info.
