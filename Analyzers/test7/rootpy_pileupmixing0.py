import numpy as np
np.random.seed(2026)

import os, sys, datetime
from six.moves import range, zip, map, filter

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency
from rootpy.tree import Tree, TreeChain
from rootpy.io import root_open
from ROOT import gROOT
gROOT.SetBatch(True)


# ______________________________________________________________________________
# Utilities

# Enums
kDT, kCSC, kRPC, kGEM, kME0 = 0, 1, 2, 3, 4


# ______________________________________________________________________________
# Analysis: picking pileup events
#
# The goal is to pick a pileup event and split it into 12 new events - there
# are 12 sectors, each sector can be used as a new event. Collect them into
# a "pool". This pool will be used as the noise to mix with the signal events
# (i.e. muon + 0PU events) to make them look like signal+noise events
# (i.e. muon + 200PU events).

class PickingAnalysis(object):
  def run(self, pileup=200):
    # Load tree
    tree = load_minbias_batch(jobid, pileup=pileup)

    # Event range
    #maxEvents = -1
    maxEvents = 1000

    # __________________________________________________________________________
    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxEvents != -1 and ievt == maxEvents:
        break

      print("Processing event: {0}".format(ievt))

      # Hits
      for ihit, hit in enumerate(evt.hits):
        print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.fr, hit.sim_phi, hit.sim_theta, hit.time, hit.sim_tp1, hit.sim_tp2))
      # Tracks
      for itrk, trk in enumerate(evt.tracks):
        print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7}".format(itrk, trk.sector, trk.mode, trk.pt, trk.phi, trk.eta, trk.theta, trk.q))
      # Particles
      for ipart, part in enumerate(evt.particles):
        print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))

    # End loop over events
    unload_tree()


# ______________________________________________________________________________
# Analysis: mixing pileup events
#
# The goal is to mix the signal (muon + 0PU) events. For each muon, find the
# sector where most of the hits from the muon belong to. Pull a random event
# from the noise pool (see PickingAnalysis). Add all the hits from signal and
# noise together. Be careful about positive or negative endcap - if signal is
# in the positive endcap, use only noise in the positive endcap, because some
# hit variables are not invariant under swapping of +/- endcap.
#
# Once the mixing is done, run the pattern recognition to find the roads (or
# proto-tracks). Save them into an output file to be used for training NN.

class MixingAnalysis(object):
  def run(self):
    # Load tree
    tree = load_pgun_batch(jobid)

    # Event range
    #maxEvents = -1
    maxEvents = 1000

    # __________________________________________________________________________
    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxEvents != -1 and ievt == maxEvents:
        break

      print("Processing event: {0}".format(ievt))

      # Hits
      for ihit, hit in enumerate(evt.hits):
        print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.fr, hit.sim_phi, hit.sim_theta, hit.time, hit.sim_tp1, hit.sim_tp2))
      # Tracks
      for itrk, trk in enumerate(evt.tracks):
        print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7}".format(itrk, trk.sector, trk.mode, trk.pt, trk.phi, trk.eta, trk.theta, trk.q))
      # Particles
      for ipart, part in enumerate(evt.particles):
        print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))

    # End loop over events
    unload_tree()


# ______________________________________________________________________________
# Settings

# Condor or not
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  os.environ['ROOTPY_GRIDMODE'] = 'true'

# Algorithm (pick one)
algo = 'default'  # phase 2
#algo = 'run3'
#algo = 'omtf'
if use_condor:
  algo = sys.argv[1]

# Analysis mode (pick one)
#analysis = 'picking'
analysis = 'mixing'
if use_condor:
  analysis = sys.argv[2]

# Job id
jobid = 0
if use_condor:
  jobid = int(sys.argv[3])


# ______________________________________________________________________________
# Input files

infile_r = None  # input file handle

def define_collections(tree):
  tree.define_collection(name='hits', prefix='vh_', size='vh_size')
  tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
  tree.define_collection(name='particles', prefix='vp_', size='vp_size')
  tree.define_collection(name='evt_info', prefix='ve_', size='ve_size')
  return

def load_tree_single(infile):
  print('[INFO] Opening file: %s' % infile)
  global infile_r
  infile_r = root_open(infile)
  tree = infile_r.ntupler.tree
  define_collections(tree)
  return tree

def load_tree_multiple(infiles):
  print('[INFO] Opening file: %s' % ' '.join(infiles))
  tree = TreeChain('ntupler/tree', infiles)
  define_collections(tree)
  return tree

def unload_tree():
  global infile_r
  try:
    infile_r.close()
  except:
    pass

def load_pgun_batch(k):
  jj = np.split(np.arange(100), 100)[k]
  infiles = []
  for j in jj:
    infiles.append('root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleMuon_Endcap_2GeV/ParticleGuns/CRAB3/190416_194707/%04i/ntuple_SingleMuon_Endcap_%i.root' % ((j+1)/1000, (j+1)))
    infiles.append('root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleMuon_Endcap2_2GeV/ParticleGuns/CRAB3/190416_194826/%04i/ntuple_SingleMuon_Endcap2_%i.root' % ((j+1)/1000, (j+1)))
  return load_tree_multiple(infiles)

def load_minbias_batch(k, pileup=200):
  if pileup == 140:
    pufiles = ['root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleNeutrino_PU140/SingleNeutrino/CRAB3/190416_160046/0000/ntuple_SingleNeutrino_PU140_%i.root' % (i+1) for i in xrange(56)]
  elif pileup == 200:
    pufiles = ['root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleNeutrino_PU200/SingleNeutrino/CRAB3/190416_160207/0000/ntuple_SingleNeutrino_PU200_%i.root' % (i+1) for i in xrange(63)]
  elif pileup == 250:
    pufiles = ['root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleNeutrino_PU250/SingleNeutrino/CRAB3/190416_160323/0000/ntuple_SingleNeutrino_PU250_%i.root' % (i+1) for i in xrange(50)]
  elif pileup == 300:
    pufiles = ['root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/ntuple_SingleNeutrino_PU300/SingleNeutrino/CRAB3/190416_160441/0000/ntuple_SingleNeutrino_PU300_%i.root' % (i+1) for i in xrange(53)]
  else:
    raise RunTimeError('Cannot recognize pileup: {0}'.format(pileup))
  #
  infile = pufiles[k]
  return load_tree_single(infile)


# ______________________________________________________________________________
# Main

if __name__ == "__main__":
  start_time = datetime.datetime.now()
  print('[INFO] Current time    : {0}'.format(start_time))
  print('[INFO] Using cmssw     : {0}'.format(os.environ['CMSSW_VERSION']))
  print('[INFO] Using condor    : {0}'.format(use_condor))
  print('[INFO] Using algo      : {0}'.format(algo))
  print('[INFO] Using analysis  : {0}'.format(analysis))
  print('[INFO] Using job id    : {0}'.format(jobid))

  if analysis == 'picking':
    analysis = PickingAnalysis()
    analysis.run()
  elif analysis == 'mixing':
    analysis = MixingAnalysis()
    analysis.run()
  else:
    raise RunTimeError('Cannot recognize analysis: {0}'.format(analysis))

  stop_time = datetime.datetime.now()
  print('[INFO] Elapsed time    : {0}'.format(stop_time - start_time))
  # DONE!
