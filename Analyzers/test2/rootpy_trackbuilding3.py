import numpy as np
np.random.seed(2023)

from itertools import izip
import sys
import time
#import concurrent.futures
from rootpy.plotting import Hist, Hist2D, Graph, Efficiency
from rootpy.tree import Tree, TreeChain, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from ROOT import gROOT
gROOT.SetBatch(True)


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
  eta_sf = np.sinh(1.9) / np.sinh(abs(eta))
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
#pt_bins = [-0.2, -0.190937, -0.180533, -0.169696, -0.158343, -0.143231, -0.123067, -0.0936418, 0.0895398, 0.123191, 0.142493, 0.157556, 0.169953, 0.180755, 0.190829, 0.2]
pt_bins = [-0.2, -0.18, -0.16, -0.133333, -0.10, -0.05, 0.05, 0.10, 0.133333, 0.16, 0.18, 0.2]
assert(len(eta_bins) == 6+1)
assert(len(pt_bins) == 11+1)

nlayers = 25  # 13 (CSC) + 9 (RPC) + 3 (GEM)
are_csc_layers = np.zeros(nlayers, dtype=np.bool)
are_csc_layers[:13] = True

safe_ipt = lambda x: np.clip(x, 0, len(pt_bins)-2)


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

# Decide EMTF hit layer number
class EMTFLayer(object):
  def __init__(self):
    self.lut = np.zeros((4,5,5,2), dtype=np.int32) - 99
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
    assert(len(indices) == nlayers)
    for i, index in enumerate(indices):
      self.lut[index] = i

  def get(self, hit):
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
    return self.lut[index]

anemtflayer = EMTFLayer()
def emtf_layer(hit):
  return anemtflayer.get(hit)

# Decide EMTF hit bend
class EMTFBend(object):
  def __init__(self):
    self.lut = (5, -5, 4, -4, 3, -3, 2, -2, 1, -1, 0)

  def get(self, hit):
    if hit.type == kCSC:
      clct = int(hit.pattern)
    else:
      clct = 10
    return self.lut[clct]

anemtfbend = EMTFBend()
def emtf_bend(hit):
  return anemtfbend.get(hit)

# Decide EMTF road quality (by pT)
class EMTFRoadQuality(object):
  def __init__(self):
    self.best_ipt = np.digitize([0.], pt_bins[1:])[0]  # skip lowest edge

  def get(self, ipt):
    return self.best_ipt - abs(ipt - self.best_ipt)

anemtfroadquality = EMTFRoadQuality()
def emtf_road_quality(ipt):
  return anemtfroadquality.get(ipt)

# Decide EMTF road sort code
class EMTFRoadSortCode(object):
  def __init__(self):
    pass

  def get(self, road_mode, road_quality, road_hits):
    #def madorsky_code(mode, qual):
    #  code = 0
    #  code |= ((mode >> 3) & 1) << 6
    #  code |= ((qual >> 2) & 1) << 5
    #  code |= ((mode >> 2) & 1) << 4
    #  code |= ((qual >> 1) & 1) << 3
    #  code |= ((mode >> 1) & 1) << 2
    #  code |= ((qual >> 0) & 1) << 1
    #  code |= ((mode >> 0) & 1) << 0
    #return code

    def mlayer_code(hits, qual):
      # 11     10     9      8      7    6      5    4    3      2..0
      # GE1/1, ME1/1, ME1/2, GE2/1, ME2, RE1&2, ME3, ME4, RE3&4, qual
      hits_to_mlayer = (10,10,9,9,9,7,7,7,7,5,5,4,4,6,6,6,3,3,3,3,3,3,11,8,11)
      code = 0
      for hit in hits:
        hit_lay = hit.emtf_layer
        mlayer = hits_to_mlayer[hit_lay]
        code |= (1 << mlayer)
      code |= qual
      return code
    return mlayer_code(road_hits, road_quality)

anemtfroadsortcode = EMTFRoadSortCode()
def emtf_road_sort_code(road_mode, road_quality, road_hits):
  return anemtfroadsortcode.get(road_mode, road_quality, road_hits)

# Decide EMTF road mode
class EMTFRoadMode(object):
  def __init__(self):
    self.singlemu = (11,13,14,15)
    self.doublemu = (7,10,12) + self.singlemu
    self.muopen = (3,5,6,9) + self.doublemu

anemtfroadmode = EMTFRoadMode()
def emtf_is_singlemu(mode):
  return mode in anemtfroadmode.singlemu
def emtf_is_doublemu(mode):
  return mode in anemtfroadmode.doublemu
def emtf_is_muopen(mode):
  return mode in anemtfroadmode.muopen

# Decide EMTF pgun event weight
class EMTFPGunWeight(object):
  def __init__(self):
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
    self.lut = c

  def get(self, part):
    ipt = np.digitize([part.invpt], pt_bins[1:])[0]  # skip lowest edge
    ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge
    ipt = safe_ipt(ipt)
    index = (ipt, ieta)
    return self.lut[index]

anemtfpgunweight = EMTFPGunWeight()
def emtf_pgun_weight(part):
  return anemtfpgunweight.get(part)


class EMTFExtrapolation(object):
  def __init__(self):
    self.theta_bins = (14, 0.5, 1.9)
    self.pt_bins = (100, -0.2, 0.2)
    self.loaded = False

  def _find_bin(self, x, bins):
    x = np.clip(x, bins[1], bins[2])
    binx = (x - bins[1]) / (bins[2] - bins[1]) * bins[0]
    return int(binx)

  def find_theta_bin(self, part):
    x = np.sinh(1.8) / np.sinh(abs(part.eta))
    return self._find_bin(x, self.theta_bins)

  def find_pt_bin(self, part):
    x = part.invpt
    return self._find_bin(x, self.pt_bins)

  def get(self, part):
    if not self.loaded:
      with np.load(bankfile) as data:
        self.lut = data['patterns_exphi']
        self.loaded = True
    index = (self.find_theta_bin(part), self.find_pt_bin(part))
    c = self.lut[index]
    exphi = part.phi + c * (part.invpt * np.sinh(1.8) / np.sinh(abs(part.eta)))
    return exphi

anemtfextrapolation = EMTFExtrapolation()
def emtf_extrapolation(part):
  return anemtfextrapolation.get(part)


# ______________________________________________________________________________
# Classes

class Particle(object):
  def __init__(self, pt, eta, phi, q):
    self.pt = pt
    self.eta = eta
    self.phi = phi
    self.q = q

  def to_parameters(self):
    parameters = np.array((np.true_divide(part.q, part.pt), part.phi, part.eta), dtype=np.float32)
    return parameters

class Pattern(object):
  def __init__(self, xmin, xmed, xmax, ymin, ymed, ymax):
    self.xmin = xmin
    self.xmed = xmed
    self.xmax = xmax
    self.ymin = ymin
    self.ymed = ymed
    self.ymax = ymax

class PatternBank(object):
  def __init__(self, bankfile):
    with np.load(bankfile) as data:
      patterns_phi = data['patterns_phi']
      patterns_theta = data['patterns_theta']
    self.x_array = patterns_phi
    self.y_array = patterns_theta
    assert(self.x_array.dtype == np.int32)
    assert(self.y_array.dtype == np.int32)
    assert(self.x_array.shape == (len(pt_bins)-1, len(eta_bins)-1, nlayers, 3))
    assert(self.y_array.shape == (len(pt_bins)-1, len(eta_bins)-1, nlayers, 3))

class Hit(object):
  def __init__(self, _id, bx, emtf_layer, emtf_phi, emtf_theta, emtf_bend, sim_tp):
    self.id = _id  # (_type, station, ring, fr)
    self.bx = bx
    self.emtf_layer = emtf_layer
    self.emtf_phi = emtf_phi
    self.emtf_theta = emtf_theta
    self.emtf_bend = emtf_bend
    self.sim_tp = sim_tp

class Road(object):
  def __init__(self, _id, hits, mode, quality, sort_code):
    self.id = _id  # (endcap, sector, ipt, ieta, iphi)
    self.hits = hits
    self.mode = mode
    self.quality = quality
    self.sort_code = sort_code

  def to_variables(self, use_sim_tp=False):
    amap = {}
    #np.random.shuffle(self.hits)  # pick a random hit for now
    for hit in self.hits:
      hit_lay = hit.emtf_layer
      if (not use_sim_tp) or (use_sim_tp and hit.sim_tp == True):
        if hit_lay not in amap:
          amap[hit_lay] = hit
    #
    if not use_sim_tp:
      tmp_thetas =[v.emtf_theta for k, v in amap.iteritems()]
      tmp_theta =  np.median(tmp_thetas, overwrite_input=True)
      for hit in self.hits:
        hit_lay = hit.emtf_layer
        hit_check = amap[hit_lay]
        if abs(hit.emtf_theta - tmp_theta) < abs(hit_check.emtf_theta - tmp_theta):
          amap[hit_lay] = hit
    #
    hits_phi = np.zeros(nlayers, dtype=np.float32) + np.nan
    hits_theta = np.zeros(nlayers, dtype=np.float32) + np.nan
    hits_bend = np.zeros(nlayers, dtype=np.float32) + np.nan
    hits_mask = np.zeros(nlayers, dtype=np.float32) + 1.0
    for lay, hit in amap.iteritems():
      hits_phi[lay] = hit.emtf_phi
      hits_theta[lay] = hit.emtf_theta
      hits_bend[lay] = hit.emtf_bend
      hits_mask[lay] = False
    #
    (endcap, sector, ipt, ieta, iphi) = self.id
    road_info = [ipt, ieta, iphi]
    variables = np.hstack((hits_phi, hits_theta, hits_bend, hits_mask, road_info))
    return variables

class Track(object):
  def __init__(self, _id, hits, mode, pt, q, emtf_phi, emtf_theta, chi2, ndof):
    self.id = _id  # (endcap, sector)
    self.hits = hits
    self.mode = mode
    self.xml_pt = pt
    self.pt = pt * (1.0 + 0.22 * 1.28155)  # erf(1.28155/sqrt(2)) = 0.8 [90% upper limit from -1 to -1]
    self.q = q
    self.emtf_phi = emtf_phi
    self.emtf_theta = emtf_theta
    self.chi2 = chi2
    self.ndof = ndof

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
      #hit_x = hit.emtf_phi - (pattern_iphi * 32 - 16)  # multiply by 'quadstrip' unit (4 * 8)
      hit_x = hit.emtf_phi - (pattern_iphi * 16 - 8)  # multiply by 'doublestrip' unit (2 * 8)
      hit_y = hit.emtf_theta
      hit_lay = hit.lay

      # Match patterns
      mask = (pattern_x[...,hit_lay,0] <= hit_x) & (hit_x <= pattern_x[...,hit_lay,2]) & (pattern_y[...,hit_lay,0] - 2 <= hit_y) & (hit_y <= pattern_y[...,hit_lay,2] + 2)

      # Create a hit (for output)
      hit_id = (hit.type, hit.station, hit.ring, hit.fr)
      hit_sim_tp = (hit.sim_tp1 == 0 and hit.sim_tp2 == 0)
      myhit = Hit(hit_id, hit.bx, hit_lay, hit.emtf_phi, hit.emtf_theta, emtf_bend(hit), hit_sim_tp)

      # Get results
      for index, condition in np.ndenumerate(mask):
        if condition:  # good hit
          ipt = ipt_range[index[0]]
          ieta = ieta_range[index[1]]
          iphi = iphi_range[index[2]]
          road_id = (endcap, sector, ipt, ieta, iphi)
          amap.setdefault(road_id, []).append(myhit)  # append hit to road

    # Create a road
    roads = []
    for road_id, road_hits in amap.iteritems():
      # Try BX window (-1,0)
      road_mode = 0
      tmp_road_hits = []
      for hit in road_hits:
        if hit.bx in (-1,0):
          (_type, station, ring, fr) = hit.id
          road_mode |= (1 << (4 - station))
          tmp_road_hits.append(hit)

      if not emtf_is_singlemu(road_mode):
        # Try BX window (0,+1)
        road_mode = 0
        tmp_road_hits = []
        for hit in road_hits:
          if hit.bx in (0,+1):
            (_type, station, ring, fr) = hit.id
            road_mode |= (1 << (4 - station))
            tmp_road_hits.append(hit)

      if emtf_is_singlemu(road_mode):
        (endcap, sector, ipt, ieta, iphi) = road_id
        road_quality = emtf_road_quality(ipt)
        road_sort_code = emtf_road_sort_code(road_mode, road_quality, tmp_road_hits)
        myroad = Road(road_id, tmp_road_hits, road_mode, road_quality, road_sort_code)
        roads.append(myroad)
    return roads

  def run(self, hits, part=None):
    roads = []

    fake_modes = np.zeros(12, dtype=np.int32)
    for ihit, hit in enumerate(hits):
      if hit.bx in (-1,0,+1):
        hit.endsec = (hit.sector - 1) if hit.endcap == 1 else (hit.sector - 1 + 6)
        hit.lay = emtf_layer(hit)
        assert(hit.lay != -99)
        if hit.type == kCSC:  # at least 2 CSC hits
          fake_modes[hit.endsec] |= (1 << (4 - hit.station))

    # Loop over sector processors
    for endcap in (-1, +1):
      for sector in (1, 2, 3, 4, 5, 6):
        endsec = (sector - 1) if endcap == 1 else (sector - 1 + 6)
        fake_mode = fake_modes[endsec]
        early_exit = sum([(fake_mode>>3) & 1, (fake_mode>>2) & 1, (fake_mode>>1) & 1, (fake_mode>>0) & 1]) < 2  # at least 2 CSC hits
        if early_exit:  continue

        # Patterns to run
        ipt_range = xrange(len(pt_bins))
        ieta_range = xrange(len(eta_bins))
        #iphi_range = xrange(4928/32)  # divide by 'quadstrip' unit (4 * 8)
        iphi_range = xrange(4928/16)  # divide by 'doublestrip' unit (2 * 8)

        # Hits
        sector_hits = [hit for hit in hits if hit.bx in (-1,0,+1) and hit.endsec == endsec]

        # Cheat using gen particle info
        if part is not None:
          #ipt_range = [x for x in xrange(part.ipt-1, part.ipt+1+1) if 0 <= x < len(pt_bins)-1]
          ipt_range = xrange(0,(len(pt_bins)-1)//2+1) if part.q < 0 else xrange((len(pt_bins)-1)//2, len(pt_bins)-1)
          ieta_range = [x for x in xrange(part.ieta-1, part.ieta+1+1) if 0 <= x < len(eta_bins)-1]
          tmp_phis = [hit.emtf_phi for hit in sector_hits if hit.type == kCSC and hit.station >= 2]
          #tmp_phi = np.mean(tmp_phis)
          tmp_phi = np.median(tmp_phis, overwrite_input=True)
          #iphi = int(tmp_phi/32)  # divide by 'quadstrip' unit (4 * 8)
          iphi = int(tmp_phi/16)  # divide by 'doublestrip' unit (2 * 8)
          iphi_range = xrange(iphi-6, iphi+6+1)

        roads_tmp = self._apply_patterns(endcap, sector, ipt_range, ieta_range, iphi_range, sector_hits)
        roads += roads_tmp
    return roads

class RoadCleaning(object):
  def __init__(self):
    pass

  # https://stackoverflow.com/a/30396816
  def _iter_from_middle(self, lst):
    try:
      middle = len(lst)//2
      yield lst[middle]
      for shift in xrange(1, middle+1):
        # order is important!
        yield lst[middle - shift]
        yield lst[middle + shift]
    except IndexError: # occures on lst[len(lst)] or for empty list
      raise StopIteration

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

  def _sortby(self, clean_roads, groupinfo):
    if clean_roads:
      clean_roads.sort(key=lambda road: road.sort_code, reverse=True)

      # Iterate over clean_roads
      for i, road in enumerate(clean_roads):
        keep = True
        gi = groupinfo[road.id]

        # No intersect between two ranges (x1, x2), (y1, y2): (x2 < y1) || (x1 > y2)
        # Intersect: !((x2 < y1) || (x1 > y2)) = (x2 >= y1) and (x1 <= y2)
        for j, road_to_check in enumerate(clean_roads[:i]):
          gj = groupinfo[road_to_check.id]
          get_endsec = lambda x: x[:2]
          # Allow +/-2 due to extrapolation-to-ME3 error
          if (get_endsec(road.id) == get_endsec(road_to_check.id)) and (gi[1] >= gj[0]-2) and (gi[0] <= gj[1]+2):
            keep = False
            break

        if keep:
          yield road
      return

  def run(self, roads):
    # road_id = (endcap, sector, ipt, ieta, iphi)
    amap = {road.id : road for road in roads}

    # pick median in each iphi group
    clean_roads = []
    groupinfo = {}
    for group in self._groupby(amap.keys()):
      for index in self._iter_from_middle(range(len(group))):
        road_id = group[index]
        road = amap[road_id]
        keep = True
        if (0 <= index-1) and road.sort_code < amap[group[index-1]].sort_code:
          keep = False
        if (index+1 < len(group)) and road.sort_code < amap[group[index+1]].sort_code:
          keep = False
        if keep:
          break
      clean_roads.append(road)
      #
      get_iphi = lambda x: x[4]
      g = (get_iphi(group[0]), get_iphi(group[-1]))  # first and last road_id's in the iphi group
      groupinfo[road_id] = g

    # sort the roads by (mode, quality)
    sorted_clean_roads = []
    for road in self._sortby(clean_roads, groupinfo):
      sorted_clean_roads.append(road)
    return sorted_clean_roads

class PtAssignment(object):
  def __init__(self, kerasfile):
    (encoder, model, model_weights, bpca) = kerasfile
    with np.load(encoder) as loaded:
      self.x_mean = loaded['x_mean']
      self.x_std  = loaded['x_std']

    with np.load(bpca) as loaded:
      self.Vt_all = loaded['Vt_all']
      self.v_k_all = loaded['v_k_all']
      self.v_n_all = loaded['v_n_all']
      self.v_mean_all = loaded['v_mean_all']
      self.v_std_all = loaded['v_std_all']

    from keras.models import load_model
    import keras.backend as K
    def huber_loss(y_true, y_pred, delta=1.345):
      x = K.abs(y_true - y_pred)
      squared_loss = 0.5*K.square(x)
      absolute_loss = delta * (x - 0.5*delta)
      xx = K.switch(x < delta, squared_loss, absolute_loss)
      #return K.sum(xx, axis=-1)
      return K.mean(xx, axis=-1)

    self.loaded_model = load_model(model, custom_objects={'huber_loss': huber_loss})
    self.loaded_model.load_weights(model_weights)

  def run(self, x):
    x_new = np.array([], dtype=np.float32)
    y = np.array([], dtype=np.float32)
    z = np.array([], dtype=np.float32)
    if len(x) == 0:
      return (x_new, y, z)

    assert(len(x.shape) == 2)
    assert(x.shape[1] == (nlayers * 4) + 3)

    self.nentries = x.shape[0]
    self.x_copy = x.copy()

    # Get views
    self.x_phi   = self.x_copy[:, nlayers*0:nlayers*1]
    self.x_theta = self.x_copy[:, nlayers*1:nlayers*2]
    self.x_bend  = self.x_copy[:, nlayers*2:nlayers*3]
    self.x_mask  = self.x_copy[:, nlayers*3:nlayers*4].astype(np.bool)  # this makes a copy
    self.x_road  = self.x_copy[:, nlayers*4:nlayers*5]  # ipt, ieta, iphi

    # Subtract median phi from hit phis
    #self.x_phi_median    = self.x_road[:, 2] * 32 - 16  # multiply by 'quadstrip' unit (4 * 8)
    self.x_phi_median    = self.x_road[:, 2] * 16 - 8  # multiply by 'doublestrip' unit (2 * 8)
    self.x_phi_median    = self.x_phi_median[:, np.newaxis]
    self.x_phi          -= self.x_phi_median

    # Subtract median theta from hit thetas
    self.x_theta_median  = np.nanmedian(self.x_theta, axis=1)
    self.x_theta_median  = self.x_theta_median[:, np.newaxis]
    self.x_theta        -= self.x_theta_median

    # Zones
    self.x_ieta  = self.x_road[:, 1].astype(np.int32)

    # Standard scales
    self.x_copy -= self.x_mean
    self.x_copy /= self.x_std

    ## Remove outlier hits by checking hit thetas
    #x_theta_tmp = np.abs(self.x_theta) > 4.0
    #self.x_phi  [x_theta_tmp] = np.nan
    #self.x_theta[x_theta_tmp] = np.nan
    #self.x_bend [x_theta_tmp] = np.nan
    #self.x_mask [x_theta_tmp] = 1.0

    ## Something wrong with GE2/1?
    #bad_ge21 = 23
    #self.x_phi  [:, bad_ge21] = np.nan
    #self.x_theta[:, bad_ge21] = np.nan
    #self.x_bend [:, bad_ge21] = np.nan
    #self.x_mask [:, bad_ge21] = 1.0

    # Add variables: theta_median and mode variables
    self.x_theta_median -= 3  # scaled to [0,1]
    self.x_theta_median /= 83
    hits_to_station = np.array((5,5,1,1,1,2,2,2,2,3,3,4,4,1,1,2,3,3,3,4,4,4,5,2,1), dtype=np.int32)  # '5' denotes ME1/1
    assert(len(hits_to_station) == nlayers)
    self.x_mode_vars = np.zeros((self.nentries, 5), dtype=np.float32)
    self.x_mode_vars[:,0] = np.any(self.x_mask[:,hits_to_station == 5] == 0, axis=1)
    self.x_mode_vars[:,1] = np.any(self.x_mask[:,hits_to_station == 1] == 0, axis=1)
    self.x_mode_vars[:,2] = np.any(self.x_mask[:,hits_to_station == 2] == 0, axis=1)
    self.x_mode_vars[:,3] = np.any(self.x_mask[:,hits_to_station == 3] == 0, axis=1)
    self.x_mode_vars[:,4] = np.any(self.x_mask[:,hits_to_station == 4] == 0, axis=1)

    # Remove NaN
    #np.nan_to_num(self.x_copy, copy=False)
    self.x_copy[np.isnan(self.x_copy)] = 0.0

    # Get x
    #x_new = self.x_phi
    x_new = np.hstack((self.x_phi, self.x_theta, self.x_bend, self.x_theta_median, self.x_mode_vars))

    # Predict y
    y = self.loaded_model.predict(x_new)

    # PCA stuff
    z = self._bpca()
    return (x_new, y, z)

  def _bpca(self):
    x_copy = np.hstack((self.x_phi, self.x_theta, self.x_bend))
    x_mask_copy = self.x_mask.copy()

    result = np.zeros((x_copy.shape[0], 3), dtype=np.float32)

    for i in xrange(x_copy.shape[0]):
      zone = self.x_ieta[i]
      Vt = self.Vt_all[zone]
      v_k = self.v_k_all[zone]
      v_n = self.v_n_all[zone]
      v_mean = self.v_mean_all[zone]
      v_std = self.v_std_all[zone]

      x_transformed = np.dot(x_copy[i], Vt.T)
      x_transformed -= v_mean
      x_transformed /= v_std
      x_transformed **= 2
      x_transformed = x_transformed[v_k:v_n]

      chi2 = np.sum(x_transformed)
      ndof = np.sum(x_copy[i] != 0.0) - v_k
      ndof_mask = np.sum(~x_mask_copy[i])
      result[i] = (chi2, ndof, ndof_mask)
    return result


class TrackProducer(object):
  def __init__(self):
    pass

  def run(self, clean_roads, variables, predictions, chi2_vars):
    assert(len(clean_roads) == len(variables))
    assert(len(clean_roads) == len(predictions))
    assert(len(clean_roads) == len(chi2_vars))

    tracks = []

    for myroad, myvars, mypreds, mychi2 in izip(clean_roads, variables, predictions, chi2_vars):
      # Unpack variables
      assert(len(myvars.shape) == 1)
      assert(myvars.shape[0] == (nlayers * 3) + 6)
      x_phi          = myvars[nlayers*0:nlayers*1]
      x_theta        = myvars[nlayers*1:nlayers*2]
      x_bend         = myvars[nlayers*2:nlayers*3]
      x_theta_median = myvars[nlayers*3]
      x_mode_vars    = myvars[nlayers*3+1:nlayers*3+6]

      trk_mode = 0
      for i, x in enumerate(x_mode_vars):
        if i == 0:
          station = 1
        else:
          station = i
        if x:
          trk_mode |= (1 << (4 - station))

      trk_ipt = np.digitize([mypreds[0]], pt_bins[1:])[0]  # skip lowest edge
      trk_ipt = safe_ipt(trk_ipt)
      quality1 = myroad.quality
      quality2 = emtf_road_quality(trk_ipt)

      bx_counter1 = 0  # count hits with BX <= -1
      bx_counter2 = 0  # count hits with BX <= 0
      for hit in myroad.hits:
        if hit.bx <= -1:
          bx_counter1 += 1
        if hit.bx <= 0:
          bx_counter2 += 1
      trk_bx_zero = (bx_counter1 < 2 and bx_counter2 >= 2)

      if emtf_is_singlemu(trk_mode) and quality2 <= (quality1+1) and trk_bx_zero:
        (endcap, sector, ipt, ieta, iphi) = myroad.id
        trk_id = (endcap, sector)
        trk_pt = np.abs(1.0/mypreds[0])
        trk_q  = np.sign(mypreds[0])
        trk = Track(trk_id, myroad.hits, trk_mode, trk_pt, trk_q, iphi, x_theta_median, mychi2[0], mychi2[1])
        tracks.append(trk)
    return tracks

# ______________________________________________________________________________
# Book histograms
histograms = {}
histogram2Ds = {}

# Efficiency
eff_pt_bins = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50., 60., 80., 120.]
for k in ["denom", "numer"]:
  hname = "eff_vs_genpt_%s" % k
  histograms[hname] = Hist(eff_pt_bins, name=hname, title="; gen p_{T} [GeV]", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_geneta_%s" % k
  histograms[hname] = Hist(26, 1.2, 2.5, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_genphi_%s" % k
  histograms[hname] = Hist(32, -3.2, 3.2, name=hname, title="; gen #phi", type='F')
  histograms[hname].Sumw2()


# ______________________________________________________________________________
# Open file
infile = 'ntuple_SingleMuon_Toy_5GeV_add.3.root'
#infile_r = root_open(infile)
#tree = infile_r.ntupler.tree
tree = TreeChain('ntupler/tree', [infile])
print('[INFO] Opening file: %s' % infile)

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='particles', prefix='vp_', size='vp_size')

# Get number of events
#maxEvents = -1
#maxEvents = 2000000
maxEvents = 10000
print('[INFO] Using max events: %i' % maxEvents)

# Analysis mode
#analysis = "training"
#analysis = "application"
#analysis = "rates"
analysis = "effie"
print('[INFO] Using analysis mode: %s' % analysis)

# Other stuff
bankfile = 'histos_tb.8.npz'

kerasfile = ['encoder.8.npz', 'model.8.h5', 'model_weights.8.h5', 'bpca.8.npz']

#pufiles = ['root://cmsxrootd.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_9_2_3_patch1/ntuple_SingleNeutrino_PU140/ParticleGuns/CRAB3/180116_214607/0000/ntuple_SingleNeutrino_PU140_%i.root' % (i+1) for i in xrange(100)]
#pufiles = ['root://cmsxrootd.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_9_2_3_patch1/ntuple_SingleNeutrino_PU200/ParticleGuns/CRAB3/180116_214738/0000/ntuple_SingleNeutrino_PU200_%i.root' % (i+1) for i in xrange(100)]
pufiles = ['root://cmsio2.rc.ufl.edu//store/user/jiafulow/L1MuonTrigger/P2_9_2_3_patch1/ntuple_SingleNeutrino_PU200/ParticleGuns/CRAB3/180116_214738/0000/ntuple_SingleNeutrino_PU200_%i.root' % (i+1) for i in xrange(100)]


# ______________________________________________________________________________
# Begin job

# Analysis: training
if analysis == "training":

  # 3-dimensional arrays of lists
  # [ipt][ieta][lay]
  patterns_phi = [[[[] for k in xrange(nlayers)] for j in xrange(len(eta_bins)-1)] for i in xrange(len(pt_bins)-1)]
  patterns_theta = [[[[] for k in xrange(nlayers)] for j in xrange(len(eta_bins)-1)] for i in xrange(len(pt_bins)-1)]
  e = EMTFExtrapolation()
  patterns_exphi = [[[] for j in xrange(e.pt_bins[0])] for i in xrange(e.theta_bins[0])]

# Analysis: application
elif analysis == "application":

  # Workers
  bank = PatternBank(bankfile)
  recog = PatternRecognition(bank)
  clean = RoadCleaning()
  out_particles = []
  out_roads = []
  npassed, ntotal = 0, 0

# Analysis: rates, effie
elif analysis == "rates" or analysis == "effie":
  pass


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
    #part.exphi = extrapolate_to_emtf(part.phi, part.invpt, part.eta)
    part.exphi = emtf_extrapolation(part)
    part.sector = find_sector(part.exphi)
    part.endcap = 1 if part.eta >= 0. else -1
    part.emtf_phi = calc_phi_loc_int(np.rad2deg(part.exphi), part.sector)
    part.emtf_theta = calc_theta_int(calc_theta_deg_from_eta(part.eta), part.endcap)

    smear = True
    if smear:
      # CLCT spatial resolution (halfstrip) = (w/2)/sqrt(12)
      pitch = 2.3271e-3  # in radians
      sigma = (pitch/2)/np.sqrt(12)
      sigma *= 2  # this is an arbitrary scale factor
      smear = sigma * np.random.normal()
      part.emtf_phi_nosmear = part.emtf_phi
      part.emtf_phi = calc_phi_loc_int(np.rad2deg(part.exphi + smear), part.sector)

    if ievt < 20:
      print("evt {0} has {1} particles and {2} hits".format(ievt, len(evt.particles), len(evt.hits)))
      print(".. part invpt: {0} pt: {1} eta: {2} phi: {3} exphi: {4} se: {5} ph: {6} th: {7}".format(part.invpt, part.pt, part.eta, part.phi, part.exphi, part.sector, part.emtf_phi, part.emtf_theta))

    part.ipt = np.digitize([part.invpt], pt_bins[1:])[0]  # skip lowest edge
    part.ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge
    part.ipt = safe_ipt(part.ipt)

    the_patterns_phi = patterns_phi[part.ipt][part.ieta]
    the_patterns_theta = patterns_theta[0][part.ieta]  # no binning in pt

    e = EMTFExtrapolation()
    the_patterns_exphi = patterns_exphi[e.find_theta_bin(part)][e.find_pt_bin(part)]

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

        if hit.type == kCSC and hit.station == 3:  # extrapolate to ME3
          dphi = delta_phi(np.deg2rad(hit.sim_phi), part.phi)
          dphi /= (part.invpt * np.sinh(1.8) / np.sinh(abs(part.eta)))
          the_patterns_exphi.append(dphi)


  # ____________________________________________________________________________
  # Analysis: application

  elif analysis == "application":

    #roads = recog.run(evt.hits)

    # Cheat using gen particle info
    part = evt.particles[0]  # particle gun
    part.invpt = np.true_divide(part.q, part.pt)
    part.ipt = np.digitize([part.invpt], pt_bins[1:])[0]  # skip lowest edge
    part.ieta = np.digitize([abs(part.eta)], eta_bins[1:])[0]  # skip lowest edge
    part.ipt = safe_ipt(part.ipt)

    roads = recog.run(evt.hits, part)
    clean_roads = clean.run(roads)

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
      print(".. part invpt: {0} pt: {1} eta: {2} phi: {3} ipt: {4} ieta: {5}".format(part.invpt, part.pt, part.eta, part.phi, part.ipt, part.ieta))
      for iroad, myroad in enumerate(roads):
        print(".. road {0} {1} {2} {3} {4} {5}".format(iroad, myroad.id, len(myroad.hits), myroad.mode, myroad.quality, myroad.sort_code))
      for iroad, myroad in enumerate(clean_roads):
        print(".. croad {0} {1} {2} {3} {4} {5}".format(iroad, myroad.id, len(myroad.hits), myroad.mode, myroad.quality, myroad.sort_code))
        for ihit, myhit in enumerate(myroad.hits):
          print(".. .. hit  {0} {1} {2} {3} {4} {5}".format(ihit, myhit.id, myhit.emtf_layer, myhit.emtf_phi, myhit.emtf_theta, myhit.emtf_bend))

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


  # ____________________________________________________________________________
  # Analysis: rates, effie

  elif analysis == "rates" or analysis == "effie":
    break


  continue  # end loop over events

#infile_r.close()


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
  print('[INFO] Creating file: histos_tb.npz')
  if True:
    patterns_phi_tmp = patterns_phi
    patterns_theta_tmp = patterns_theta
    patterns_phi = np.zeros((len(pt_bins)-1, len(eta_bins)-1, nlayers, 3), dtype=np.int32)
    patterns_theta = np.zeros((len(pt_bins)-1, len(eta_bins)-1, nlayers, 3), dtype=np.int32)
    #
    for i in xrange(len(pt_bins)-1):
      for j in xrange(len(eta_bins)-1):
        for k in xrange(nlayers):
          patterns_phi_tmp_ijk = patterns_phi_tmp[i][j][k]
          if len(patterns_phi_tmp_ijk) > 1000:
            patterns_phi_tmp_ijk.sort()
            x = np.percentile(patterns_phi_tmp_ijk, [5, 50, 95])
            if k == 22 or k == 23 or k == 24:  # keep more GEMs
              x = np.percentile(patterns_phi_tmp_ijk, [3.5, 50, 96.5])
            x = [int(round(xx)) for xx in x]
            while (x[2] - x[0]) < 32:  # make sure the range is larger than twice the 'doublestrip' unit
              x[0] -= 1
              x[2] += 1
            patterns_phi[i,j,k] = x

          #patterns_theta_tmp_ijk = patterns_theta_tmp[i][j][k]
          patterns_theta_tmp_ijk = patterns_theta_tmp[0][j][k]  # no binning in pt
          if len(patterns_theta_tmp_ijk) > 1000:
            patterns_theta_tmp_ijk.sort()
            x = np.percentile(patterns_theta_tmp_ijk, [2, 50, 98])
            x = [int(round(xx)) for xx in x]
            patterns_theta[i,j,k] = x

    # Mask layers by (ieta, lay)
    valid_layers = [
      (0,2), (0,3), (0,4), (0,7), (0,8), (0,10), (0,12), (0,13), (0,14), (0,15), (0,18), (0,21),
      (1,2), (1,3), (1,7), (1,8), (1,10), (1,12), (1,13), (1,14), (1,15), (1,17), (1,18), (1,20), (1,21), (1,22),
      (2,0), (2,1), (2,5), (2,6), (2,9), (2,10), (2,11), (2,12), (2,17), (2,20), (2,22), (2,23),
      (3,0), (3,1), (3,5), (3,6), (3,9), (3,11), (3,16), (3,19), (3,20), (3,22), (3,23),
      (4,0), (4,1), (4,5), (4,6), (4,9), (4,11), (4,16), (4,19), (4,22), (4,23),
      (5,0), (5,1), (5,5), (5,6), (5,9), (5,11), (5,16), (5,19), (5,23),
    ]
    mask = np.ones_like(patterns_phi, dtype=np.bool)
    for valid_layer in valid_layers:
      mask[:,valid_layer[0],valid_layer[1],:] = False
    patterns_phi[mask] = 0
    patterns_theta[mask] = 0

    overwrite_extrapolation = False
    smooth_extrapolation = True
    if overwrite_extrapolation:
      e = EMTFExtrapolation()
      patterns_exphi_tmp = patterns_exphi
      patterns_exphi_tmp = np.asarray(patterns_exphi_tmp)
      patterns_exphi = np.zeros(patterns_exphi_tmp.shape, dtype=np.float32)
      for index, x in np.ndenumerate(patterns_exphi_tmp):
        if x:
          patterns_exphi[index] = np.median(x, overwrite_input=True)
      if smooth_extrapolation:
        from scipy.interpolate import Rbf
        patterns_exphi_tmp = patterns_exphi
        patterns_exphi = np.zeros(patterns_exphi_tmp.shape, dtype=np.float32)
        x = [e.pt_bins[1] + (i+0.5)/e.pt_bins[0]*(e.pt_bins[2] - e.pt_bins[1]) for i in xrange(e.pt_bins[0])]
        for index in np.ndindex(patterns_exphi_tmp.shape[:1]):
          rbf = Rbf(x, patterns_exphi_tmp[index], smooth = 0.3, function='multiquadric')
          patterns_exphi[index] = rbf(x)
    else:
      with np.load(bankfile) as data:
        patterns_exphi = data['patterns_exphi']

    outfile = 'histos_tb.npz'
    np.savez_compressed(outfile, patterns_phi=patterns_phi, patterns_theta=patterns_theta, patterns_exphi=patterns_exphi)


# Analysis: application
elif analysis == "application":

  # Plot histograms
  print('[INFO] Creating file: histos_tba.root')
  with root_open('histos_tba.root', 'recreate') as f:
    for hname in ["eff_vs_genpt", "eff_vs_geneta", "eff_vs_genphi"]:
      denom = histograms[hname + "_denom"]
      numer = histograms[hname + "_numer"]
      eff = Efficiency(numer, denom, name=hname)
      eff.SetStatisticOption(0)  # kFCP
      eff.SetConfidenceLevel(0.682689492137)  # one sigma
      eff.Write()
    print('[INFO] npassed/ntotal: %i/%i = %f' % (npassed, ntotal, float(npassed)/ntotal))

  # Save objects
  print('[INFO] Creating file: histos_tba.npz')
  if True:
    assert(len(out_particles) == npassed)
    assert(len(out_roads) == npassed)
    parameters = np.zeros((npassed, 3), dtype=np.float32)
    variables = np.zeros((npassed, (nlayers * 4) + (3)), dtype=np.float32)
    for i, (part, road) in enumerate(izip(out_particles, out_roads)):
      parameters[i] = part.to_parameters()
      #variables[i] = road.to_variables()
      variables[i] = road.to_variables(use_sim_tp=True)
    remove_zeros = ~variables[:,nlayers*3:nlayers*4].all(axis=1)  # when hits_mask is all 1
    outfile = 'histos_tba.npz'
    np.savez_compressed(outfile, parameters=parameters[remove_zeros], variables=variables[remove_zeros])


# ______________________________________________________________________________
# Analysis: rates

elif analysis == "rates":

  # Workers
  bank = PatternBank(bankfile)
  recog = PatternRecognition(bank)
  clean = RoadCleaning()
  ptassign = PtAssignment(kerasfile)
  trkprod = TrackProducer()

  def make_rates(infile_j):
    infile_r = root_open(infile_j)
    tree = infile_r.ntupler.tree

    # Define collection
    tree.define_collection(name='hits', prefix='vh_', size='vh_size')
    tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
    tree.define_collection(name='particles', prefix='vp_', size='vp_size')

    # Outputs
    out = {}
    out["nevents"] = []
    for m in ["emtf", "emtf2023"]:
      out["highest_%s_absEtaMin0_absEtaMax2.5_qmin12_pt" % m] = []

    # Loop over events
    for ievt, evt in enumerate(tree):

      roads = recog.run(evt.hits)
      clean_roads = clean.run(roads)
      variables = np.array([road.to_variables() for road in clean_roads], dtype=np.float32)
      variables_1, predictions, chi2_vars = ptassign.run(variables)
      emtf2023_tracks = trkprod.run(clean_roads, variables_1, predictions, chi2_vars)

      if ievt < 20 and False:
        print("evt {0} has {1} roads, {2} clean roads, {3} tracks, {4} old tracks".format(ievt, len(roads), len(clean_roads), len(emtf2023_tracks), len(evt.tracks)))
        for ipart, part in enumerate(evt.particles):
          if part.pt > 5.:
            part.invpt = np.true_divide(part.q, part.pt)
            print(".. part {0} {1} {2} {3} {4}".format(ipart, part.invpt, part.pt, part.eta, part.phi))
        for iroad, myroad in enumerate(clean_roads):
          print(".. croad {0} {1} {2} {3} {4} {5}".format(iroad, myroad.id, len(myroad.hits), myroad.mode, myroad.quality, myroad.sort_code))
          #for ihit, myhit in enumerate(myroad.hits):
          #  print(".. .. hit  {0} {1} {2} {3} {4} {5}".format(ihit, myhit.id, myhit.emtf_layer, myhit.emtf_phi, myhit.emtf_theta, myhit.emtf_bend))
        for itrk, mytrk in enumerate(emtf2023_tracks):
          print(".. trk {0} {1} {2} {3} {4} {5} {6}".format(itrk, mytrk.id, len(mytrk.hits), mytrk.mode, mytrk.pt, mytrk.chi2, mytrk.ndof))
          for ihit, myhit in enumerate(mytrk.hits):
            print(".. .. hit  {0} {1} {2} {3} {4} {5}".format(ihit, myhit.id, myhit.emtf_layer, myhit.emtf_phi, myhit.emtf_theta, myhit.emtf_bend))
        for itrk, mytrk in enumerate(evt.tracks):
          print(".. otrk {0} {1} {2}".format(itrk, mytrk.mode, mytrk.pt))


      # Fill histograms
      out["nevents"].append(1.0)

      def fill_highest_pt():
        highest_pt = -999999.
        for itrk, trk in enumerate(tracks):
          if select(trk):
            if highest_pt < trk.pt:  # using scaled pT
              highest_pt = trk.pt
        if highest_pt > 0.:
          highest_pt = min(100.-1e-3, highest_pt)
          out[hname].append(highest_pt)

      select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in (11,13,14,15))
      tracks = evt.tracks
      hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"
      fill_highest_pt()

      select = lambda trk: trk and (trk.chi2/trk.ndof) <= 8.5
      tracks = emtf2023_tracks
      hname = "highest_emtf2023_absEtaMin0_absEtaMax2.5_qmin12_pt"
      fill_highest_pt()

      continue  # end loop over events

    #infile_r.close()
    return out

  def make_rates_endjob(outputs):
    outputs_copy = []
    for j, out in enumerate(outputs):
      print("Retrieving job: {0}".format(j))
      outputs_copy.append(out)

    print('[INFO] Creating file: histos_tbb.root')
    with root_open('histos_tbb.root', 'recreate') as f:
      histograms_1 = {}
      hname = "nevents"
      histograms_1[hname] = Hist(5, 0, 5, name=hname, title="; count", type='F')
      hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"
      histograms_1[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
      hname = "highest_emtf2023_absEtaMin0_absEtaMax2.5_qmin12_pt"
      histograms_1[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

      for k, v in histograms_1.iteritems():
        for j, out in enumerate(outputs_copy):
          for x in out[k]:
            v.fill(x)
        v.Write()
    return

  # ____________________________________________________________________________
  # Parallel processing

  #with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
  #  outputs = executor.map(make_rates, pufiles[:4], chunksize=4)
  #  make_rates_endjob(outputs)

  #outputs = map(make_rates, pufiles[:4])
  #make_rates_endjob(outputs)

  j = 0
  out = make_rates(pufiles[j])
  make_rates_endjob([out])


# ______________________________________________________________________________
# Analysis: effie

elif analysis == "effie":

  # Workers
  bank = PatternBank(bankfile)
  recog = PatternRecognition(bank)
  clean = RoadCleaning()
  ptassign = PtAssignment(kerasfile)
  trkprod = TrackProducer()

  def make_effie(evt_range):
    evt = next(iter(tree))

    # Outputs
    out = {}
    for m in ["emtf", "emtf2023"]:
      out["%s_eff_vs_genpt_l1pt20" % m] = []
      out["%s_eff_vs_geneta_l1pt20" % m] = []
      out["%s_l1pt_vs_genpt" % m] = []
      out["%s_l1ptres_vs_genpt" % m] = []

    # Loop over events
    for ievt in evt_range:
      tree.GetEntry(ievt)

      part = evt.particles[0]  # particle gun
      part.invpt = np.true_divide(part.q, part.pt)

      roads = recog.run(evt.hits)
      clean_roads = clean.run(roads)
      variables = np.array([road.to_variables() for road in clean_roads], dtype=np.float32)
      variables_1, predictions, chi2_vars = ptassign.run(variables)
      emtf2023_tracks = trkprod.run(clean_roads, variables_1, predictions, chi2_vars)

      if ievt < 20 and False:
        print("evt {0} has {1} roads, {2} clean roads, {3} tracks".format(ievt, len(roads), len(clean_roads), len(emtf2023_tracks)))
        for itrk, mytrk in enumerate(emtf2023_tracks):
          y = np.true_divide(part.q, part.pt)
          y_pred = np.true_divide(mytrk.q, mytrk.xml_pt)
          print(".. {0} {1}".format(y, y_pred))

      # Fill histograms
      def fill_efficiency():
        trigger = any([select(trk) for trk in tracks])  # using scaled pT
        out[hname1].append((part.pt, trigger))
        if part.pt > 20.:
          out[hname2].append((abs(part.eta), trigger))

      def fill_resolution():
        if len(tracks) > 0:
          trk = tracks[0]
          trk.invpt = np.true_divide(trk.q, trk.xml_pt)  # using unscaled pT
          out[hname1].append((part.invpt, trk.invpt))
          out[hname2].append((abs(part.invpt), (abs(trk.invpt) - abs(part.invpt))/abs(part.invpt)))

      select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in (11,13,14,15)) and (trk.pt > 20.)
      tracks = evt.tracks
      hname1 = "emtf_eff_vs_genpt_l1pt20"
      hname2 = "emtf_eff_vs_geneta_l1pt20"
      fill_efficiency()
      hname1 = "emtf_l1pt_vs_genpt"
      hname2 = "emtf_l1ptres_vs_genpt"
      fill_resolution()

      select = lambda trk: trk and (trk.chi2/trk.ndof) <= 8.5 and (trk.pt > 20.)
      tracks = emtf2023_tracks
      hname1 = "emtf2023_eff_vs_genpt_l1pt20"
      hname2 = "emtf2023_eff_vs_geneta_l1pt20"
      fill_efficiency()
      hname1 = "emtf2023_l1pt_vs_genpt"
      hname2 = "emtf2023_l1ptres_vs_genpt"
      fill_resolution()

      continue  # end loop over events

    return out

  def make_effie_endjob(outputs):
    outputs_copy = []
    for j, out in enumerate(outputs):
      print("Retrieving job: {0}".format(j))
      outputs_copy.append(out)

    print('[INFO] Creating file: histos_tbc.root')
    with root_open('histos_tbc.root', 'recreate') as f:
      histograms_1 = {}
      for m in ["emtf", "emtf2023"]:
        for k in ["denom", "numer"]:
          hname = "%s_eff_vs_genpt_l1pt20_%s" % (m,k)
          histograms_1[hname] = Hist(eff_pt_bins, name=hname, title="; gen p_{T} [GeV]", type='F')
          hname = "%s_eff_vs_geneta_l1pt20_%s" % (m,k)
          histograms_1[hname] = Hist(26, 1.2, 2.5, name=hname, title="; gen |#eta| {gen p_{T} > 20 GeV}", type='F')
        hname = "%s_l1pt_vs_genpt" % m
        histograms_1[hname] = Hist2D(100, -0.3, 0.3, 300, -0.3, 0.3, name=hname, title="; gen 1/p_{T} [1/GeV]; 1/p_{T} [1/GeV]", type='F')
        hname = "%s_l1ptres_vs_genpt" % m
        histograms_1[hname] = Hist2D(100, -0.3, 0.3, 300, -2, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; #Delta(p_{T})/p_{T}", type='F')

      for hname in ["emtf_eff_vs_genpt_l1pt20", "emtf_eff_vs_geneta_l1pt20", "emtf2023_eff_vs_genpt_l1pt20", "emtf2023_eff_vs_geneta_l1pt20"]:
        denom = histograms_1[hname + "_denom"]
        numer = histograms_1[hname + "_numer"]
        for j, out in enumerate(outputs_copy):
          for x, trigger in out[hname]:
            denom.fill(x)
            if trigger:
              numer.fill(x)
        denom.Write()
        numer.Write()
        eff = Efficiency(numer, denom, name=hname)
        eff.SetStatisticOption(0)  # kFCP
        eff.SetConfidenceLevel(0.682689492137)  # one sigma
        eff.Write()

      for hname in ["emtf_l1pt_vs_genpt", "emtf_l1ptres_vs_genpt", "emtf2023_l1pt_vs_genpt", "emtf2023_l1ptres_vs_genpt"]:
        h = histograms_1[hname]
        for j, out in enumerate(outputs_copy):
          for x, y in out[hname]:
            h.fill(x, y)
        h.Write()
    return


  # ____________________________________________________________________________
  # Parallel processing

  #with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
  #  outputs = executor.map(make_effie, [range(j*10, (j+1)*10) for j in xrange(4)], chunksize=4)
  #  make_effie_endjob(outputs)

  j = 0
  out = make_effie(xrange(j*10000, (j+1)*10000))
  make_effie_endjob([out])

