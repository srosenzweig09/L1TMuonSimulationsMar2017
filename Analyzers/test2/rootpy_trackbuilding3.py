import numpy as np
np.random.seed(2023)

from itertools import izip
import time
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
eta_bins = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
pt_bins = [-0.2, -0.190937, -0.180533, -0.169696, -0.158343, -0.143231, -0.123067, -0.0936418, 0.0895398, 0.123191, 0.142493, 0.157556, 0.169953, 0.180755, 0.190829, 0.2]
assert(len(eta_bins) == 6+1)
assert(len(pt_bins) == 15+1)

nlayers = 25  # 13 (CSC) + 9 (RPC) + 3 (GEM)


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
        (1,1,1,1),  # ME1/1f
        (1,1,1,0),  # ME1/1r
        (1,1,2,1),  # ME1/2f
        (1,1,2,0),  # ME1/2r
        (1,1,3,0),  # ME1/3
        (1,2,1,1),  # ME2/1f
        (1,2,1,0),  # ME2/1r
        (1,2,2,1),  # ME2/2f
        (1,2,2,0),  # ME2/2r
        (1,3,1,0),  # ME3/1
        (1,3,2,0),  # ME3/2
        (1,4,1,0),  # ME4/1
        (1,4,2,0),  # ME4/2
        # RPC
        (2,1,2,1),  # RE1/2f
        (2,1,2,0),  # RE1/2r
        (2,2,2,0),  # RE2/2
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
      for i, index in enumerate(indices):
        lut[index] = i
      assert(np.max(lut) == nlayers-1)
      cache["lut"] = lut
      #print "LUT: ", lut
    lut = cache["lut"]  # get the LUT
    index = f(x)  # get the index of the LUT
    return lut[index]
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
  elif hit.type == kCSC and hit.station == 2:  # special case: ME2/*
    hit_fr = int(hit.fr)
  elif hit.type == kRPC and hit.station == 1:  # special case: RE1/*
    hit_fr = int(hit.fr)
  else:
    hit_fr = 0
  index = (hit.type, hit.station, hit_ring, hit_fr)
  return index


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
      #print "LUT: ", c
    lut = cache["lut"]
    index = f(x)
    return lut[index]
  return wrapper

@wrapper_emtf_pgun_weight
def emtf_pgun_weight(part):
  # Builds the index of the LUT
  # The LUT is cached in wrapper_emtf_pgun_weight()
  ipt = np.digitize([np.true_divide(part.q, part.pt)], pt_bins[1:])[0]  # skip lowest edge
  ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge
  index = (ipt, ieta)
  return index

# ______________________________________________________________________________
# Classes

class Particle(object):
  def __init__(self, pt, eta, phi, q):
    self.pt = pt
    self.eta = eta
    self.phi = phi
    self.q = q

class Pattern(object):
  def __init__(self, xmin, xmed, xmax, ymin, ymed, ymax):
    self.xmin = xmin
    self.xmed = xmed
    self.xmax = xmax
    self.ymin = ymin
    self.ymed = ymed
    self.ymax = ymax

class PatternBank(object):
  def __init__(self, patterns_phi, patterns_theta):
    self.x_array = np.array(patterns_phi, dtype=np.int32) # make a copy
    self.y_array = np.array(patterns_theta, dtype=np.int32)
    assert(self.x_array.shape == (len(pt_bins)-1, len(eta_bins)-1, nlayers, 3))
    assert(self.y_array.shape == (len(pt_bins)-1, len(eta_bins)-1, nlayers, 3))

class Hit(object):
  def __init__(self, _id, emtf_layer, emtf_phi, emtf_theta, emtf_bend):
    self.id = _id
    self.emtf_layer = emtf_layer
    self.emtf_phi = emtf_phi
    self.emtf_theta = emtf_theta
    self.emtf_bend = emtf_bend

class Road(object):
  def __init__(self, _id, hits, mode, quality):
    self.id = _id
    self.hits = hits
    self.mode = mode
    self.quality = quality

class PatternRecognition(object):
  def __init__(self, bank):
    self.bank = bank

  def _apply_patterns(self, endcap, sector, ipt_range, ieta_range, iphi_range, sector_hits):
    amap = {}

    # Retrieve patterns with (ipt, ieta, lay, pattern)
    ipt_slice = slice(ipt_range[0], ipt_range[-1]+1, None)
    ieta_slice = slice(ieta_range[0], ieta_range[-1]+1, None)
    pattern_x = self.bank.x_array[ipt_slice, ieta_slice, np.newaxis, :, :]
    pattern_y = self.bank.y_array[ipt_slice, ieta_slice, np.newaxis, :, :]
    pattern_iphi = np.arange(iphi_range[0], iphi_range[-1]+1, dtype=np.int32)

    # Loop over hits
    for ihit, hit in enumerate(sector_hits):
      # Make hit coordinates
      hit_x = hit.emtf_phi - pattern_iphi * 16  # multiply by 'doublestrip' unit (2 * 8)
      hit_y = hit.emtf_theta
      hit_lay = hit.lay

      # Match patterns
      mask = (pattern_x[...,hit_lay,0] <= hit_x) & (hit_x < pattern_x[...,hit_lay,2]) & (pattern_y[...,hit_lay,0] - 2 <= hit_y) & (hit_y < pattern_y[...,hit_lay,2] + 2)

      # Create a hit (for output)
      hit_id = (hit.type, hit.station, hit.ring, hit.fr)
      myhit = Hit(hit_id, hit_lay, hit.emtf_phi, hit.emtf_theta, hit.pattern)

      # Get results
      for index, condition in np.ndenumerate(mask):
        if condition:  # good hit
          ipt, ieta, iphi = index
          road_id = (endcap, sector, ipt, ieta, iphi)
          amap.setdefault(road_id, []).append(myhit)  # append hit to road

    # Create a road
    roads = []
    for road_id, road_hits in amap.iteritems():
      road_mode = 0
      for hit in road_hits:
        _type, station, ring, fr = hit.id
        road_mode |= (1 << (4 - station))
      #
      endcap, sector, ipt, ieta, iphi = road_id
      road_quality = 15//2 - abs(ipt - 15//2)
      #
      if road_mode in (11, 13, 14, 15):  # single-muon modes
        myroad = Road(road_id, road_hits, road_mode, road_quality)
        roads.append(myroad)
    return roads

  def run(self, hits, part=None):
    roads = []

    fake_modes = np.zeros(12, dtype=np.int32)
    for ihit, hit in enumerate(hits):
      hit.endsec = (hit.sector - 1) if hit.endcap == 1 else (hit.sector - 1 + 6)
      hit.lay = emtf_layer(hit)
      assert(hit.lay != -99)
      if hit.type == kCSC:  # at least 2 CSC hits
        fake_modes[hit.endsec] |= (1 << (4 - hit.station))

    # Loop over sector processors
    for endcap in [-1, +1]:
      for sector in [1, 2, 3, 4, 5, 6]:
        endsec = (sector - 1) if endcap == 1 else (sector - 1 + 6)
        fake_mode = fake_modes[endsec]
        early_exit = sum([(fake_mode>>3) & 1, (fake_mode>>2) & 1, (fake_mode>>1) & 1, (fake_mode>>0) & 1]) < 2  # at least 2 CSC hits
        if early_exit:  continue

        # Patterns to run
        ipt_range = xrange(len(pt_bins))
        ieta_range = xrange(len(eta_bins))
        iphi_range = xrange(4800/16)  # divide by 'doublestrip' unit (2 * 8)

        # Hits
        sector_hits = [hit for hit in hits if hit.endsec == endsec]

        # Cheat using gen particle info
        if part is not None:
          ipt_range = [part.ipt]
          ieta_range = [part.ieta]
          ipt_range = [x for x in xrange(part.ipt-1, part.ipt+1+1) if 0 <= x < len(pt_bins)-1]
          ieta_range = [x for x in xrange(part.ieta-1, part.ieta+1+1) if 0 <= x < len(eta_bins)-1]
          tmp_phis = [hit.emtf_phi for hit in sector_hits if hit.type == kCSC]
          tmp_phi = np.mean(tmp_phis)
          iphi = int(tmp_phi/16)
          iphi_range = xrange(iphi-4, iphi+4+1)

        roads_tmp = self._apply_patterns(endcap, sector, ipt_range, ieta_range, iphi_range, sector_hits)
        roads += roads_tmp
    return roads

class RoadCleaning(object):
  def __init__(self):
    pass

  def _groupby(self, data):
    def is_adjacent(prev, curr, length):
      return prev[:-1] == curr[:-1] and (prev[-1] + length) == curr[-1]

    if data:
      data.sort()
      myiter = iter(data)
      prev = curr = next(myiter)
      # Iterate over data
      while True:
        group = []
        stop = False
        # Iterate until the next value is different
        while is_adjacent(prev, curr, len(group)):
          try:
            group.append(curr)
            curr = next(myiter)
          except StopIteration:
            stop = True
            break
        # Output group
        yield group
        prev = curr
        if stop:
          return

  def run(self, roads):
    amap = {road.id : road for road in roads}

    clean_roads = []
    for group in self._groupby(amap.keys()):
      assert(len(group))
      pivot = (len(group) - 1) // 2
      road_id = group[pivot]
      road = amap[road_id]
      clean_roads.append(road)
    return clean_roads


# ______________________________________________________________________________
# Book histograms
histograms = {}
histogram2Ds = {}

# Efficiency
eff_pt_bins = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50., 60., 70., 100.]
for k in ["denom", "numer"]:
  hname = "eff_vs_genpt_%s" % k
  histograms[hname] = Hist(eff_pt_bins, name=hname, title="; gen p_{T} [GeV]", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_geneta_%s" % k
  histograms[hname] = Hist(26, 1.2, 2.5, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_genphi_%s" % k
  histograms[hname] = Hist(32, -3.2, 3.2, name=hname, title="; gen |#phi|", type='F')
  histograms[hname].Sumw2()


# ______________________________________________________________________________
# Open file
infile = root_open('ntuple_SingleMuon_Toy_5GeV_add.3.root')
tree = infile.ntupler.tree
print('[INFO] Opening file: %s' % infile)

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='particles', prefix='vp_', size='vp_size')

# Get number of events
#maxEvents = -1
#maxEvents = 1000000
maxEvents = 10000
print('[INFO] Using max events: %i' % maxEvents)

# Analysis mode
#analysis = "training"
analysis = "application"
print('[INFO] Using analysis mode: %s' % analysis)

# Other stuff
bankfile = 'histos_tb.3.pkl'


# ______________________________________________________________________________
# Begin job

# Analysis: training
if analysis == "training":

  # 3-dimensional arrays of lists
  # [ipt][ieta][lay]
  patterns_phi = [[[[] for k in xrange(nlayers)] for j in xrange(len(eta_bins)-1)] for i in xrange(len(pt_bins)-1)]
  patterns_theta = [[[[] for k in xrange(nlayers)] for j in xrange(len(eta_bins)-1)] for i in xrange(len(pt_bins)-1)]

# Analysis: application
elif analysis == "application":

  import pickle
  with open(bankfile, 'rb') as f:
    data = pickle.load(f)
    patterns_phi, patterns_theta = data

  bank = PatternBank(patterns_phi, patterns_theta)
  recog = PatternRecognition(bank)
  clean = RoadCleaning()
  out_particles = []
  out_roads = []
  npassed, ntotal = 0, 0


# ______________________________________________________________________________
# Loop over events

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
  # Analysis: training

  if analysis == "training":

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
      quadstrip = (4 * 8) / np.sqrt(12)
      doublestrip = (2 * 8) / np.sqrt(12)
      smear = doublestrip * np.random.normal()
      part.emtf_phi_nosmear = part.emtf_phi
      part.emtf_phi += smear

    if ievt < 20:
      print("evt {0} has {1} particles and {2} hits".format(ievt, len(evt.particles), len(evt.hits)))
      print(".. part invpt: {0} eta: {1} phi: {2} exphi: {3} se: {4} ph: {5} th: {6}".format(part.invpt, part.eta, part.phi, part.exphi, part.sector, part.emtf_phi, part.emtf_theta))

    part.ipt = np.digitize([np.true_divide(part.q, part.pt)], pt_bins[1:])[0]  # skip lowest edge
    part.ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge
    the_patterns_phi = patterns_phi[part.ipt][part.ieta]
    the_patterns_theta = patterns_theta[part.ipt][part.ieta]

    #pgun_weight = emtf_pgun_weight(part)

    # Loop over hits
    for ihit, hit in enumerate(evt.hits):
      lay = emtf_layer(hit)
      assert(lay != -99)
      if ievt < 20:
        print(".. hit {0} type: {1} st: {2} ri: {3} fr: {4} lay: {5} se: {6} ph: {7} th: {8} tp: {9}".format(ihit, hit.type, hit.station, hit.ring, hit.fr, lay, hit.sector, hit.emtf_phi, hit.emtf_theta, hit.sim_tp1))

      if hit.sector == part.sector and hit.bx == 0 and hit.sim_tp1 == 0 and hit.sim_tp2 == 0:
        the_patterns_phi[lay].append(hit.emtf_phi - part.emtf_phi)
        #the_patterns_theta[lay].append(hit.emtf_theta - part.emtf_theta)
        the_patterns_theta[lay].append(hit.emtf_theta)

    #the_patterns_phi[nlayers].append(part.emtf_phi)
    #the_patterns_theta[nlayers].append(part.emtf_theta)


  # ____________________________________________________________________________
  # Analysis: application

  elif analysis == "application":

    #roads = recog.run(evt.hits)

    # Cheat using gen particle info
    part = evt.particles[0]  # particle gun
    part.invpt = np.true_divide(part.q, part.pt)
    part.ipt = np.digitize([np.true_divide(part.q, part.pt)], pt_bins[1:])[0]  # skip lowest edge
    part.ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge

    roads = recog.run(evt.hits, part)

    clean_roads = clean.run(roads)
    clean_roads.sort(key=lambda x: (x.mode, x.quality), reverse=True)

    if len(clean_roads) > 0:
      mypart = Particle(
        pt = part.pt,
        eta = part.eta,
        phi = part.phi,
        q = part.q,
      )
      out_particles.append(mypart)
      out_roads.append(clean_roads[0])

    if ievt < 20:
      print("evt {0} has {1} roads and {2} clean roads".format(ievt, len(roads), len(clean_roads)))
      print(".. part invpt: {0} eta: {1} phi: {2} ipt: {3} ieta: {4}".format(part.invpt, part.eta, part.phi, part.ipt, part.ieta))
      for iroad, myroad in enumerate(roads):
        print(".. road {0} {1} {2} {3} {4}".format(iroad, myroad.id, len(myroad.hits), myroad.mode, myroad.quality))

      for iroad, myroad in enumerate(clean_roads):
        print(".. croad {0} {1} {2} {3} {4}".format(iroad, myroad.id, len(myroad.hits), myroad.mode, myroad.quality))
        for ihit, myhit in enumerate(myroad.hits):
          print(".. .. hit  {0} {1} {2} {3} {4}".format(ihit, myhit.id, myhit.emtf_phi, myhit.emtf_theta, myhit.emtf_bend))

    # Quick efficiency
    trigger = len(clean_roads) > 0

    hname = "eff_vs_genpt_denom"
    histograms[hname].fill(part.pt)
    if trigger:
      hname = "eff_vs_genpt_numer"
      histograms[hname].fill(part.pt)

    if part.pt > 20.:
      hname = "eff_vs_geneta_denom"
      histograms[hname].fill(abs(part.eta))
      if trigger:
        hname = "eff_vs_geneta_numer"
        histograms[hname].fill(abs(part.eta))

      hname = "eff_vs_genphi_denom"
      histograms[hname].fill(part.phi)
      if trigger:
        hname = "eff_vs_genphi_numer"
        histograms[hname].fill(part.phi)

    # Keep statistics
    ntotal += 1
    if trigger:
      npassed += 1

  continue  # end loop over events

maxEvents = ievt


# ______________________________________________________________________________
# End job

# Analysis: training
if analysis == "training":

  # Plot histograms
  print('[INFO] Creating file: histos_tb.root')
  with root_open('histos_tb.root', 'recreate') as f:
    for i in xrange(len(pt_bins)-1):
      for j in xrange(len(eta_bins)-1):
        for k in xrange(nlayers):
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
  print('[INFO] Creating file: histos_tb.pkl')
  with open('histos_tb.pkl', 'wb') as f:
    for i in xrange(len(pt_bins)-1):
      for j in xrange(len(eta_bins)-1):
        for k in xrange(nlayers):
          if len(patterns_phi[i][j][k]) > (0.001 * maxEvents):
            patterns_phi[i][j][k].sort()
            x = np.percentile(patterns_phi[i][j][k], [5, 50, 95])
            x = [int(round(xx)) for xx in x]
            patterns_phi[i][j][k] = x
          else:
            patterns_phi[i][j][k] = [0, 0, 0]

          if len(patterns_theta[i][j][k]) > (0.001 * maxEvents):
            patterns_theta[i][j][k].sort()
            x = np.percentile(patterns_theta[i][j][k], [2, 50, 98])
            x = [int(round(xx)) for xx in x]
            patterns_theta[i][j][k] = x
          else:
            patterns_theta[i][j][k] = [0, 0, 0]
    patterns_phi = np.array(patterns_phi, dtype=np.int32)
    patterns_theta = np.array(patterns_theta, dtype=np.int32)
    pickle.dump([patterns_phi, patterns_theta], f)


# Analysis: application
elif analysis == "application":

  # Plot histograms
  print('[INFO] Creating file: histos_tba.root')
  with root_open('histos_tba.root', 'recreate') as f:
    for hname in ["eff_vs_genpt", "eff_vs_geneta", "eff_vs_genphi"]:
      hname_numer = "%s_%s" % (hname, "numer")
      hname_denom = "%s_%s" % (hname, "denom")
      eff = Efficiency(histograms[hname_numer], histograms[hname_denom], name=hname)
      eff.Write()
    print('[INFO] npassed/ntotal: %i/%i = %f' % (npassed, ntotal, np.true_divide(npassed, ntotal)))

  # Save objects
  import pickle
  print('[INFO] Creating file: histos_tba.pkl')
  with open('histos_tba.pkl', 'wb') as f:
    assert(len(out_particles) == npassed)
    assert(len(out_roads) == npassed)
    pickle.dump([out_particles, out_roads], f)
