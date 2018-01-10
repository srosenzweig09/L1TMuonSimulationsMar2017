import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from ROOT import gROOT
gROOT.SetBatch(True)


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

def range_phi_deg(deg):
  while deg <  -180.:
    deg += 360.
  while deg >= +180.:
    deg -= 360.
  return deg

def calc_phi_loc_deg_from_glob(glob, sector):
  # glob in deg, sector [1-6]
  glob = range_phi_deg(glob)
  loc = glob - 15. - (60. * (sector-1))
  return loc

def calc_phi_loc_int(glob, sector):
  # glob in deg, sector [1-6]
  loc = calc_phi_loc_deg_from_glob(glob, sector)
  if (loc + 22.) < 0.:
    loc += 360.
  loc = (loc + 22.) * 60.
  phi_int = int(round(loc))
  return phi_int

def calc_theta_int(theta, endcap):
  # theta in deg, endcap [-1,+1]
  if endcap == -1:
    theta = 180. - theta
  theta = (theta - 8.5) * 128./(45.0-8.5)
  theta_int = int(round(theta))
  return theta_int

def calc_theta_rad_from_eta(eta):
  theta = np.arctan2(1.0, np.sinh(eta))
  return theta

def calc_theta_deg_from_eta(eta):
  return np.rad2deg(calc_theta_rad_from_eta(eta))

def extrapolate_to_emtf(phi, invpt, eta):  # phi in radians
  # 1.204 is the magic constant at eta of 1.9
  eta_sf = np.sinh(1.9) / np.sinh(eta)
  return phi - 1.204 * invpt * eta_sf

def find_sector(phi):  # phi in radians
  dphi = delta_phi(phi, np.pi/12)  # sector 1 starts at 15 deg
  dphi = int(np.floor(dphi/(np.pi/3)))  # divide by 60 deg
  if dphi < 0:
    sector = 7 + dphi
  else:
    sector = 1 + dphi
  return sector


# Globals
eta_bins = [1.2, 1.40943, 1.58427, 1.76857, 1.94529, 2.15444, 2.5]
pt_bins = [-0.2, -0.191121, -0.181153, -0.1684, -0.156863, -0.144086, -0.125538, -0.0946667, 0.0784762, 0.116727, 0.138507, 0.152444, 0.1666, 0.177728, 0.190803, 0.2]
assert(len(eta_bins) == 6+1)
assert(len(pt_bins) == 15+1)

#nlayers = (9+10+3)*2  # (CSC+RPC+GEM)x(F/R)
nlayers = 25


# More functions

# From https://stackoverflow.com/a/31539746
def weighted_percentile(data, percents, weights=None):
  ''' percents in units of 1%
  weights specifies the frequency (count) of data.
  '''
  if weights is None:
    return np.percentile(data, percents)
  ind=np.argsort(data)
  d=data[ind]
  w=weights[ind]
  p=1.*w.cumsum()/w.sum()*100
  y=np.interp(percents, p, d)
  return y


def wrapper_emtf_layer(f):
  # [hit_type][station][ring][fr]
  # hit_type = DT (0), CSC (1), RPC (2), GEM (3)
  # station = 1, 2, 3, 4
  # ring = 1, 2, 3, 4
  # fr = 0, 1
  cache = {}
  def wrapper(x):
    if "lut" not in cache:
      lut = np.zeros((4,5,5,2), dtype=np.int32) - 99
      indices = [
        # CSC
        (1,1,1,0),  # ME1/1r
        (1,1,1,1),  # ME1/1f
        (1,1,2,0),  # ME1/2r
        (1,1,2,1),  # ME1/2f
        (1,1,3,0),  # ME1/3
        (1,2,1,0),  # ME2/1r
        (1,2,1,1),  # ME2/1f
        (1,2,2,0),  # ME2/2
        (1,3,1,0),  # ME3/1
        (1,3,2,0),  # ME3/2
        (1,4,1,0),  # ME4/1
        (1,4,2,0),  # ME4/2
        # RPC
        (2,1,2,0),  # RE1/2
        (2,1,3,0),  # RE1/3
        (2,2,2,0),  # RE2/2
        (2,2,3,0),  # RE2/3
        (2,3,1,0),  # RE3/1
        (2,3,2,0),  # RE3/2
        (2,3,3,0),  # RE3/3
        (2,4,1,0),  # RE4/1
        (2,4,2,0),  # RE4/2
        (2,4,3,0),  # RE4/3
        # GEM and ME0
        (3,1,1,0),  # GE1/1
        (3,2,1,0),  # GE2/1
        (3,1,4,0),  # ME0
      ]
      for i, ind in enumerate(indices):
        lut[ind] = i
      assert(np.max(lut) == nlayers-1)
      cache["lut"] = lut
      print "LUT: ", lut  #FIXME
    lut = cache["lut"]  # get the LUT
    ind = f(x)  # get the index of the LUT
    return lut[ind]
  return wrapper

@wrapper_emtf_layer
def emtf_layer(hit):
  # Builds the index of the LUT.
  # The LUT is cached in wrapper_emtf_layer()
  if hit.type == kCSC and hit.station == 1 and hit.ring == 4:  # special case: ME1/1a
    hit_ring = 1
  else:
    hit_ring = hit.ring
  if hit.type == kCSC and hit.station == 1:  # special case: ME1/*
    hit_fr = int(hit.fr)
  elif hit.type == kCSC and hit.station == 2 and hit.ring == 1:  # special case: ME2/1
    hit_fr = int(hit.fr)
  else:
    hit_fr = 0
  ind = (hit.type, hit.station, hit_ring, hit_fr)
  return ind


def wrapper_emtf_pgun_weight(f):
  # [ipt][ieta]
  # ipt = ...
  # ieta = ...
  cache = {}
  def wrapper(x):
    if "lut" not in cache:
      a = pt_bins
      b = eta_bins
      na = len(a)
      nb = len(b)
      c = np.zeros((na+1,nb+1), dtype=np.float32)

      for i, _ in np.ndenumerate(c):
        ia = i[0]
        ib = i[1]
        if ia == 0 or ib == 0 or ia == na or ib == nb:  continue
        xa = (a[ia] - a[ia-1]) / (a[-1] - a[0])
        xb = (b[ib] - b[ib-1]) / (b[-1] - b[0])
        c[i] = (xa * xb)  # weight
        c[i] = 1.0/c[i]  # 1/weight

      cache["lut"] = c
      print "LUT: ", c  #FIXME
    lut = cache["lut"]
    ind = f(x)
    return lut[ind]
  return wrapper

@wrapper_emtf_pgun_weight
def emtf_pgun_weight(part):
  # Builds the index of the LUT
  # The LUT is cached in wrapper_emtf_pgun_weight()
  ipt = np.digitize([part.invpt], pt_bins)[0]
  ieta = np.digitize([abs(part.eta)], eta_bins)[0]
  ind = (ipt, ieta)
  return ind

# Book histograms
histograms = {}
histogram2Ds = {}


# ______________________________________________________________________________
# Open file
infile = root_open('ntuple_SingleMuon_Toy_5GeV_add.2.root')
tree = infile.ntupler.tree
print "[INFO] Opening file: %s" % infile

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='particles', prefix='vp_', size='vp_size')

# Get number of events
#maxEvents = -1
maxEvents = 400000
#maxEvents = 10000
print "[INFO] Using max events: %i" % maxEvents

# ______________________________________________________________________________
# Loop over events

patterns_phi = []
patterns_theta = []

for i in xrange(len(pt_bins)):
  patterns_phi.append([])
  patterns_theta.append([])
  for j in xrange(len(eta_bins)):
    patterns_phi[-1].append([])
    patterns_theta[-1].append([])
    for lay in xrange(nlayers+1):
      patterns_phi[-1][-1].append([])
      patterns_theta[-1][-1].append([])

for ievt, evt in enumerate(tree):
  if maxEvents != -1 and ievt == maxEvents:
    break

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

  part = evt.particles[0]  # particle gun

  part.invpt = np.true_divide(part.q, part.pt)
  part.exphi = extrapolate_to_emtf(part.phi, part.invpt, part.eta)
  part.sector = find_sector(part.exphi)
  part.endcap = 1 if part.eta >= 0. else -1
  part.emtf_phi = calc_phi_loc_int(np.rad2deg(part.exphi), part.sector)
  part.emtf_theta = calc_theta_int(calc_theta_deg_from_eta(part.eta), part.endcap)

  smear = True
  if smear:
    # one sector is 60 deg + 20 deg from neighbor
    # one emtf_phi unit is 1/60 deg
    # so one sector covers 80 * 60 = 4800 units
    quadstrip = 4 * 8 / np.sqrt(12)
    doublestrip = 2 * 8 / np.sqrt(12)
    smear = doublestrip * np.random.normal()
    part.emtf_phi_nosmear = part.emtf_phi
    part.emtf_phi += smear

  if ievt < 20:
    print("evt {0} has {1} particles and {2} hits".format(ievt, len(evt.particles), len(evt.hits)))
    print(".. part invpt: {0} eta: {1} phi: {2} exphi: {3} se: {4} ph: {5} th: {6}".format(part.invpt, part.eta, part.phi, part.exphi, part.sector, part.emtf_phi, part.emtf_theta))

  ipt = np.digitize([part.invpt], pt_bins)[0]
  ieta = np.digitize([abs(part.eta)], eta_bins)[0]
  the_patterns_phi = patterns_phi[ipt][ieta]
  the_patterns_theta = patterns_theta[ipt][ieta]

  if ievt < 20:
    print(".. .. ipt: {0} ieta: {1}".format(ipt, ieta))

  #pgun_weight = emtf_pgun_weight(part)

  # Loop over hits
  for ihit, hit in enumerate(evt.hits):
    lay = emtf_layer(hit)
    assert(lay != -99)
    if ievt < 20:
      print(".. hit {0} type: {1} st: {2} ri: {3} fr: {4} lay: {5} se: {6} ph: {7} th: {8} tp: {9}".format(ihit, hit.type, hit.station, hit.ring, hit.fr, lay, hit.sector, hit.emtf_phi, hit.emtf_theta, hit.sim_tp1))

    if hit.sector == part.sector and hit.sim_tp1 == 0:
      the_patterns_phi[lay].append(hit.emtf_phi - part.emtf_phi)
      the_patterns_theta[lay].append(hit.emtf_theta - part.emtf_theta)

      if ievt < 20:
        print(".. .. dphi: {0} dtheta: {1}".format(hit.emtf_phi - part.emtf_phi, hit.emtf_theta - part.emtf_theta))

  the_patterns_phi[nlayers].append(part.emtf_phi)
  the_patterns_theta[nlayers].append(part.emtf_theta)


# ______________________________________________________________________________
# End job

#check = np.zeros((len(pt_bins), len(eta_bins), nlayers+1), dtype=np.int32)
#for i in xrange(len(pt_bins)):
#  for j in xrange(len(eta_bins)):
#    for k in xrange(nlayers+1):
#      check[(i, j, k)] = len(patterns_phi[i][j][k])


# Plot histograms
with root_open("histos_tb.root", "recreate") as f:
  for i in xrange(len(pt_bins)):
    for j in xrange(len(eta_bins)):
      for k in xrange(nlayers+1):
        hname = "patterns_phi_%i_%i_%i" % (i,j,k)
        h1a = Hist(201, -402, 402, name=hname, title=hname, type='F')
        for x in patterns_phi[i][j][k]:  h1a.fill(x)
        h1a.Write()

        hname = "patterns_theta_%i_%i_%i" % (i,j,k)
        h1b = Hist(81, -40.5, 40.5, name=hname, title=hname, type='F')
        for x in patterns_theta[i][j][k]:  h1b.fill(x)
        h1b.Write()

# Save objects
import pickle
with open('histos_tb.pkl', 'wb') as f:
  pickle.dump([patterns_phi, patterns_theta], f)
