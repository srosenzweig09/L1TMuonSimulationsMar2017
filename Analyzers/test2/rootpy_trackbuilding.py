import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open


# ______________________________________________________________________________
# Tree models
#   see: http://www.rootpy.org/auto_examples/tree/model.html

class Hit(TreeModel):
  pass

class Track(TreeModel):
  pass

class Particle(TreeModel):
  pass


# ______________________________________________________________________________
# Analyzer

# Open file
infile = root_open('ntuple.4.root')
tree = infile.ntupler.tree
#maxEvents = -1
maxEvents = 20

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM, kTT = 0, 1, 2, 3, 20

# Functions
def delta_phi(lhs, rhs):  # in degrees
  deg = lhs - rhs
  while deg <  -180.:  deg += 360.
  while deg >= +180.:  deg -= 360.
  return deg

def get_pt_bin(pt):
  ipt = -1
  if ((1.0/20 - 0.005) < 1.0/pt <= (1.0/20)):
    ipt = 4
  return ipt

def get_eta_bin(eta):
  ieta = -1
  if (1.64 < abs(eta) <= 1.74):
    ieta = 0
  elif (1.74 < abs(eta) <= 1.84):
    ieta = 1
  elif (1.84 < abs(eta) <= 1.94):
    ieta = 2
  elif (1.94 < abs(eta) <= 2.04):
    ieta = 3
  elif (2.04 < abs(eta) <= 2.14):
    ieta = 4
  return ieta


# Book histograms
histograms = {}
histogram2Ds = {}

# CSC
for i in xrange(4):
  for j in xrange(5):
    hname = "dphi_csc_st%i_eta%i" % (i,j)
    histograms[hname] = Hist(100, -360., 360., name=hname, title="; #Delta#phi [deg]", type='F')

# TT
for i in xrange(5):
  for j in xrange(5):
    hname = "dphi_tt_st%i_eta%i" % (i,j)
    histograms[hname] = Hist(100, -360., 360., name=hname, title="; #Delta#phi [deg]", type='F')


# ______________________________________________________________________________
# Loop over events
for ievt, evt in enumerate(tree):
  if maxEvents != -1 and ievt == maxEvents:
    break

  # ____________________________________________________________________________
  # Verbose

  verbose = True

  if verbose:
    if (ievt % 1 == 0):  print("Processing event: {0}".format(ievt))

    # Hits
    for ihit, hit in enumerate(evt.hits):
      print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sim_phi, hit.sim_theta, hit.fr))
    # Tracks
    for itrk, trk in enumerate(evt.tracks):
      print(".. trk  {0} {1} {2} {3} {4} {5} {6}".format(itrk, trk.pt, trk.phi, trk.eta, trk.theta, trk.q, trk.mode))
    # Gen particles
    for ipart, part in enumerate(evt.genparticles):
      print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))
  else:
    if (ievt % 1000 == 0):  print("Processing event: {0}".format(ievt))

  # ____________________________________________________________________________
  # Make plots

  no_genparticles = (len(evt.genparticles) == 0)

  if no_genparticles: continue
  assert len(evt.genparticles) == 1

  mypart = evt.genparticles[0]
  ipt = get_pt_bin(mypart.pt)
  ieta = get_eta_bin(mypart.eta)

  if ipt == -1 or ieta == -1: continue



  continue  # end loop over event

