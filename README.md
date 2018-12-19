# L1TMuonSimulationsMar2017

Package to do the Phase 2 studies for Endcap Muon Track Finder (EMTF).

[![Build Status](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017.svg)](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017)
[![CMSSW version](https://img.shields.io/badge/cmssw-CMSSW__10__1__5-002963.svg)](https://github.com/cms-sw/cmssw)
[![Latest tag](https://img.shields.io/github/tag/jiafulow/L1TMuonSimulationsMar2017.svg)](https://github.com/jiafulow/L1TMuonSimulationsMar2017)

## Build

```shell
export SCRAM_ARCH=slc6_amd64_gcc630
scram p -n P2_CMSSW_10_1_5 CMSSW_10_1_5
cd P2_CMSSW_10_1_5/src
cmsenv

# Checkout the emulator in the 'phase2-develop' branch
git clone -b phase2-develop git@github.com:jiafulow/DataFormatsSep2016.git DataFormats
git clone -b phase2-develop git@github.com:jiafulow/L1TriggerSep2016.git L1Trigger
# Checkout this repository
git clone git@github.com:jiafulow/L1TMuonSimulationsMar2017 L1TMuonSimulations
# Compile
scram b -j 8
```

## Develop

Please do not work on the 'master' branch directly. Create a new branch for new features.

## Version

- v1.1.0 (2018-11-06): First results.
