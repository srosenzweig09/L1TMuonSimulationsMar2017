"""Data analysis (quick)."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtf_algos import *


# ______________________________________________________________________________
# Analyses

class _BaseAnalysis(object):
  """Abstract base class"""
  pass

# ______________________________________________________________________________
class DummyAnalysis(_BaseAnalysis):
  """Dummy analysis.

  Description.
  """

  def run(self, algo):
    # Load tree
    tree = load_pgun_test()
    verbosity = 1

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      if verbosity >= kINFO:
        print('evt {0} has {1} particles and {2} hits'.format(ievt, len(evt.particles), len(evt.hits)))

      # Particles
      part = evt.particles[0]  # particle gun
      if verbosity >= kINFO:
        print('.. part {0} {1} {2} {3} {4} {5}'.format(0, part.pt, part.eta, part.phi, part.invpt, part.d0))

      # Sim hits
      if verbosity >= kINFO:
        for isimhit, simhit in enumerate(evt.simhits):
          simhit_id = (simhit.type, simhit.layer, simhit.chamber)
          print('.. simhit {0} {1} {2} {3}'.format(isimhit, simhit_id, simhit.phi, simhit.theta))

      # Trigger primitives
      if verbosity >= kINFO:
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.type, hit.station, hit.ring, find_endsec(hit.endcap, hit.sector), hit.fr, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.type == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6}'.format(ihit, hit_id, hit.emtf_phi, hit.emtf_theta, hit.bend, hit.quality, hit_sim_tp))

    # End loop over events
    return

# ______________________________________________________________________________
class ZoneAnalysis(_BaseAnalysis):
  """Find zone boundaries.

  Description.
  """

  def run(self, algo, pileup=200):
    # Overwrite eta bins
    eta_bins = (0.8, 1.2, 1.55, 1.98, 2.5)
    eta_bins = eta_bins[::-1]

    def find_eta_bin(eta):
      ieta = np.digitize((np.abs(eta),), eta_bins[1:])[0]  # skip lowest edge
      ieta = np.clip(ieta, 0, len(eta_bins)-2)
      return ieta

    nzones = len(eta_bins) - 1
    out = {}  # dict of dict
    for zone in range(nzones):
      out[zone] = {}

    # __________________________________________________________________________
    # Load tree
    tree = load_pgun_batch(jobid)
    verbosity = 1

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Particles (pT > 4 GeV)
      part = evt.particles[0]  # particle gun
      if not (part.pt > 4):
        continue

      zone = find_eta_bin(part.eta)

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        lay = find_emtf_layer(hit)
        d = out[zone]
        if lay not in d:
          d[lay] = []
        d[lay].append(hit.emtf_theta)

    # End loop over events

    # __________________________________________________________________________
    # Print results
    for zone in range(nzones):
      d = out[zone]
      keys = sorted(d.keys())
      for k in keys:
        lay = k
        alist = d[lay]
        n = len(alist)
        if n > 100:
          p = np.percentile(alist, [1,2,2.5,3], overwrite_input=True)
          q = np.percentile(alist, [97,97.5,98,99], overwrite_input=True)
          print(zone, '%03i' % lay, '%5i' % n, p, q)
    return


# ______________________________________________________________________________
# Main

import os, sys, datetime

# Algorithm (pick one)
algo = 'default'  # phase 2
#algo = 'run3'

# Analysis mode (pick one)
#analysis = 'dummy'
analysis = 'zone'

# Job id (pick an integer)
jobid = 0

# Max num of events (-1 means all events)
maxevents = -1

# Condor or not
# if 'CONDOR_EXEC' is defined, overwrite the 3 arguments (algo, analysis, jobid)
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  os.environ['ROOTPY_GRIDMODE'] = 'true'
  algo = sys.argv[1]
  analysis = sys.argv[2]
  jobid = int(sys.argv[3])
  maxevents = -1

# Main function
def main():
  start_time = datetime.datetime.now()
  print('[INFO] Current time    : {0}'.format(start_time))
  print('[INFO] Using cmssw     : {0}'.format(os.environ['CMSSW_VERSION']))
  print('[INFO] Using condor    : {0}'.format(use_condor))
  print('[INFO] Using algo      : {0}'.format(algo))
  print('[INFO] Using analysis  : {0}'.format(analysis))
  print('[INFO] Using jobid     : {0}'.format(jobid))
  print('[INFO] Using maxevents : {0}'.format(maxevents))

  if analysis == 'dummy':
    anna = DummyAnalysis()
    anna.run(algo=algo)

  elif analysis == 'zone':
    anna = ZoneAnalysis()
    anna.run(algo=algo)

  else:
    raise RuntimeError('Cannot recognize analysis: {0}'.format(analysis))

  # DONE!
  stop_time = datetime.datetime.now()
  print('[INFO] Elapsed time    : {0}'.format(stop_time - start_time))

if __name__ == '__main__':
  main()
