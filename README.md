# L1TMuonSimulationsMar2017

Software package to do the Phase 2 studies for the Level-1 Endcap Muon Track Finder (EMTF).

[![Build Status](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017.svg)](https://travis-ci.org/jiafulow/L1TMuonSimulationsMar2017)
[![CMSSW version](https://img.shields.io/badge/cmssw-CMSSW__10__6__3-002963.svg)](https://github.com/cms-sw/cmssw)
[![Latest tag](https://img.shields.io/github/tag/jiafulow/L1TMuonSimulationsMar2017.svg)](https://github.com/jiafulow/L1TMuonSimulationsMar2017)

## Build

### SLC7

``` shell
export SCRAM_ARCH=slc7_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram p -n P2_CMSSW_10_6_3 CMSSW_10_6_3
cd P2_CMSSW_10_6_3/src
cmsenv

# Checkout the emulator in the 'phase2-develop' branch
git clone -b phase2-develop git@github.com:jiafulow/DataFormatsSep2016.git DataFormats
git clone -b phase2-develop git@github.com:jiafulow/L1TriggerSep2016.git L1Trigger
# Checkout this repository
git clone git@github.com:jiafulow/L1TMuonSimulationsMar2017 L1TMuonSimulations
# Compile
scram b -j 10
```

### SLC6

```shell
export SCRAM_ARCH=slc6_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
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

Please do not develop on the 'master' branch. The 'master' branch is frequently changed/rebased. Create a new branch and do any development there.

## Version

- v3.1.0 (2019-09-13): Move to CMSSW_10_6_3 with the change to use ME0TriggerDigi and GEMPadDigiCluster. Switch to use PhaseIITDRSpring19 samples. No change in NN.

- v3.0.0 (2019-09-04): Remove PU discr. Retrain NN to do 3 parameters at the same time: pT, displaced pT, and d0.

- v2.0.0 (2019-06-27): Codes are re-organized/re-factored so that they can be ported to CMSSW.

- v1.5.0 (2019-05-06): Add displaced muon patterns. NN development incomplete.

- v1.4.0 (2019-03-20): Add feasiblity study with displaced muons.

- v1.3.0 (2019-03-12): Results using CMSSW_10_4_0 including overlap region.

- v1.2.0 (2019-02-19): Move to CMSSW_10_4_0. Extend to overlap region. The NN has been modified to use 36 features and 30/25/20 nodes in the 3 hidden layers.

- v1.1.0 (2018-11-06): First results. Using CMSSW_10_1_5. Include GEM, ME0, iRPC. Add CSC bend info.
