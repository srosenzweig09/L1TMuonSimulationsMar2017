# L1TMuonSimulationsMar2017

Package to do Phase 2 studies for L1 Endcap Muon Track Finder (EMTF).

[![Build Status](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017.svg)](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017)
[![CMSSW version](https://img.shields.io/badge/cmssw-CMSSW__9__2__X-002963.svg)](https://github.com/cms-sw/cmssw)
[![Latest tag](https://img.shields.io/github/tag/jiafulow/L1TMuonSimulationsMar2017.svg)](https://github.com/jiafulow/L1TMuonSimulationsMar2017)

## Build

```shell
export SCRAM_ARCH=slc6_amd64_gcc530
scram p -n P2_CMSSW_9_2_3_patch1 CMSSW_9_2_3_patch1
cd P2_CMSSW_9_2_3_patch1/src
cmsenv

# L1TriggerSep2016 using branch 'phase2-develop'
git clone -b phase2-develop git@github.com:jiafulow/L1TriggerSep2016.git L1Trigger
# DataFormatsSep2016 using branch 'l1t-integration-CMSSW_9_2_0'
git clone -b l1t-integration-CMSSW_9_2_0 git@github.com:jiafulow/DataFormatsSep2016.git DataFormats
# This repository
git clone git@github.com:jiafulow/L1TMuonSimulationsMar2017 L1TMuonSimulations
# Compile
scram b -j 8
```

## Develop

Please do not work on the 'master' branch directly. Create a new branch for new features.

## Versions

- v0.0.X: Initial development
