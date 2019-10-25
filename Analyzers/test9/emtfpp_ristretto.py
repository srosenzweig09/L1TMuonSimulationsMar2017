"""Data preparation."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtfpp_algos import *


# ______________________________________________________________________________
# Classes

class _BaseAnalysis(object):
  """Abstract base class"""
  pass

# ______________________________________________________________________________
class SignalAnalysis(_BaseAnalysis):
  """Prepare signal data used for training.

  Description.
  """

  def run(self, algo):
    return

# ______________________________________________________________________________
class NoiseAnalysis(_BaseAnalysis):
  """Prepare noise data used for training.

  Description.
  """

  def run(self, algo):
    return

# ______________________________________________________________________________
class EffieAnalysis(_BaseAnalysis):
  """Prepare muon+200PU data used for evaluating efficiency.

  Description.
  """

  def run(self, algo):
    return

# ______________________________________________________________________________
class RatesAnalysis(_BaseAnalysis):
  """Prepare neutrino+200PU data used for evaluating trigger rates.

  Description.
  """

  def run(self, algo):
    return


# ______________________________________________________________________________
# Main

import os, sys, datetime

# Algorithm (pick one)
algo = 'default'  # phase 2
#algo = 'run3'

# Analysis mode (pick one)
analysis = 'signal'
#analysis = 'noise'
#analysis = 'effie'
#analysis = 'rates'

# Job id (pick an integer)
jobid = 0

# Max num of events (-1 means all events)
maxevents = 100

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

  if analysis == 'signal':
    anna = SignalAnalysis()
    anna.run(algo=algo)

  elif analysis == 'noise':
    anna = NoiseAnalysis()
    anna.run(algo=algo)

  elif analysis == 'effie':
    anna = EffieAnalysis()
    anna.run(algo=algo)

  elif analysis == 'rates':
    anna = RatesAnalysis()
    anna.run(algo=algo)

  else:
    raise RuntimeError('Cannot recognize analysis: {0}'.format(analysis))

  # DONE!
  stop_time = datetime.datetime.now()
  print('[INFO] Elapsed time    : {0}'.format(stop_time - start_time))

if __name__ == '__main__':
  main()
