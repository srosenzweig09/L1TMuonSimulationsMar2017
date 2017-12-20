"""
# Setup working directory (first time)
mkdir muon_work
cd muon_work
export SCRAM_ARCH=slc6_amd64_gcc630
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_9_4_1
cd CMSSW_9_4_1/src
cmsenv
cd ../..
VIRTUALENV_PATH=venv
virtualenv $VIRTUALENV_PATH
source $VIRTUALENV_PATH/bin/activate
pip install --upgrade pip
pip install uproot
"""

"""
# Setup working directory (next)
cd muon_work
cd CMSSW_9_4_1/src
cmsenv
cd ../..
VIRTUALENV_PATH=venv
source $VIRTUALENV_PATH/bin/activate
"""

import numpy as np
np.random.seed(2023)

import uproot
from itertools import izip
from rootpy.plotting import Hist, Hist2D
from rootpy.tree import Tree, TreeModel, TreeChain, FloatCol, IntCol, ShortCol
from rootpy.io import root_open

# ______________________________________________________________________________
# Setup
class Event(object):
  pass

class TreeCollectionObject(object):
  pass

class TreeCollection(object):
  def __init__(self, keys, values, prefix=""):
    assert(len(values) > 0)
    self.prefix = prefix
    self.size = len(values[0])
    self.collection = [TreeCollectionObject() for i in xrange(self.size)]
    for k, v in izip(keys, values):  # for each variable
      assert(len(v) == self.size)
      k1 = self.remove_prefix(k, self.prefix)
      for obj, v1 in izip(self.collection, v):  # for each object
        setattr(obj, k1, v1)
  def remove_prefix(self, text, prefix):
    return text[text.startswith(prefix) and len(prefix):]
  def __len__(self):
    return self.size
  def __iter__(self):
    for index in range(len(self)):
      yield self.__getitem__(index)
  def __getitem__(self, index):
    return self.collection[index]

# ______________________________________________________________________________
# Analyzer

# Enums
kDT, kCSC, kRPC, kGEM, kTT = 0, 1, 2, 3, 20

# Lambdas
deg_to_rad = lambda x: x * np.pi/180.

rad_to_deg = lambda x: x * 180./np.pi

# Functions
def delta_phi(lhs, rhs):  # in radians
  rad = lhs - rhs
  while rad <  -np.pi:  rad += np.pi*2
  while rad >= +np.pi:  rad -= np.pi*2
  return rad

def delta_theta(lhs, rhs):  # in radians
  rad = lhs - rhs
  return rad

def select_by_eta(eta):
  return 1.24 <= abs(eta) < 2.4

def select_by_bx(bx):
  return bx == 0

def select_by_vertex(vx, vy, vz):
  return np.sqrt(vx*vx + vy*vy) < 15. and abs(vz) < 50.

# Book histograms
histograms = {}
histogram2Ds = {}

hname, htitle = "muon_eta_vs_pt", "; 1/p_{T} [1/GeV]; |#eta|"
histogram2Ds[hname] = Hist2D(100, -0.2, 0.2, 65, 1.2, 2.5, name=hname, title=htitle, type='F')


# ______________________________________________________________________________
# Open file with uproot
#print 'Using numpy ver {0}'.format(np.__version__)  # 1.12.1
#print 'Using uproot ver {0}'.format(uproot.__version__)  # 2.5.15

f = uproot.open('rateplots_mc_r305310_run2_all.1.root')
tree = f["ntupler/tree"]
#print tree.keys()
print '[INFO] Opening file: {0}'.format(f.name)

# Define collections
vh_names = tree.allkeys(filtername=lambda x: x.startswith("vh_") and x != "vh_size")
vt_names = tree.allkeys(filtername=lambda x: x.startswith("vt_") and x != "vt_size")
vp_names = tree.allkeys(filtername=lambda x: x.startswith("vp_") and x != "vp_size")

#vh_arrays = tree.arrays(vh_names, outputtype=tuple)
#vt_arrays = tree.arrays(vt_names, outputtype=tuple)
#vp_arrays = tree.arrays(vp_names, outputtype=tuple)

vh_arrays = tree.lazyarrays(vh_names, outputtype=tuple)
vt_arrays = tree.lazyarrays(vt_names, outputtype=tuple)
vp_arrays = tree.lazyarrays(vp_names, outputtype=tuple)

vh_iter = iter(izip(*vh_arrays))
vt_iter = iter(izip(*vt_arrays))
vp_iter = iter(izip(*vp_arrays))

# Get number of events
#maxEvents = -1
maxEvents = 100000
ievt = 0
nevts = len(vh_arrays[0])
print '[INFO] Getting num of events: {0}'.format(nevts)


# ______________________________________________________________________________
# Loop over events
for ievt in xrange(nevts):
  if maxEvents != -1 and ievt == maxEvents:
    break

  # ____________________________________________________________________________
  # Setup

  #evt = Event()
  #evt.hits = TreeCollection(vh_names, next(vh_iter), prefix="vh_")
  #evt.tracks = TreeCollection(vt_names, next(vt_iter), prefix="vt_")
  #evt.particles = TreeCollection(vp_names, next(vp_iter), prefix="vp_")

  # ____________________________________________________________________________
  # Verbose

  verbose = False

  if verbose:
    if (ievt % 1 == 0):  print("Processing event: {0}".format(ievt))
    # Hits
    for ihit, hit in enumerate(evt.hits):
      print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.fr, hit.sim_phi, hit.sim_theta, hit.sim_tp1, hit.sim_tp2))
    # Tracks
    for itrk, trk in enumerate(evt.tracks):
      print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7}".format(itrk, trk.sector, trk.mode, trk.pt, trk.phi, trk.eta, trk.theta, trk.q))
    # Gen particles
    for ipart, part in enumerate(evt.particles):
      print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))
  else:
    if (ievt % 1000 == 0):  print("Processing event: {0}".format(ievt))

  # ____________________________________________________________________________
  # Analysis

  #h = histogram2Ds["muon_eta_vs_pt"]
  #for ipart, part in enumerate(evt.particles):
  #  if part.pt > 5.:
  #    if select_by_eta(part.eta):
  #      if select_by_bx(part.bx):
  #        if select_by_vertex(part.vx, part.vy, part.vz):
  #          h.fill(float(part.q)/part.pt, abs(part.eta))

  h = histogram2Ds["muon_eta_vs_pt"]
  vp_data = next(vp_iter)
  nvariables = len(vp_data)
  nparticles = len(vp_data[0])
  def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

  for ipart in xrange(nparticles):
    part = TreeCollectionObject()
    for k, v in izip(vp_names, vp_data):
      k1 = remove_prefix(k, "vp_")
      v1 = v[ipart]
      setattr(part, k1, v1)

    if part.pt > 5.:
      if select_by_eta(part.eta):
        if select_by_bx(part.bx):
          if select_by_vertex(part.vx, part.vy, part.vz):
            h.fill(float(part.q)/part.pt, abs(part.eta))




# ______________________________________________________________________________
# End job

with root_open("histos_uproot.root", "recreate") as f:
  for k, v in histograms.iteritems():
    v.Write()
  for k, v in histogram2Ds.iteritems():
    v.Write()

