import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, TreeChain, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from ROOT import gROOT
from datetime import datetime


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

mystate = 1

gROOT.SetBatch(1)

# Open file
if mystate == 0:
  tree_name_i = '/home/jlow/L1MuonTrigger/CRAB3/P2_9_2_3_patch1/crab_projects/crab_ntuple_SingleMuon_PositiveEndCap/results/ntuple_SingleMuon_PositiveEndCap_%i.root'
  tree = TreeChain('ntupler/tree', [(tree_name_i % i) for i in range(1,8+1)])
elif mystate == 1:
  #tree_name_i = '/home/jlow/L1MuonTrigger/CRAB3/P2_9_2_3_patch1/crab_projects_old/crab_ntuple_SingleNeutrino_PU140_170720_232549/results/ntuple_SingleNeutrino_PU140_%i.root'
  #tree = TreeChain('ntupler/tree', [(tree_name_i % i) for i in range(1,90+1)])
  tree_name_i = '/home/jlow/L1MuonTrigger/CRAB3/P2_9_2_3_patch1/crab_projects/crab_ntuple_SingleNeutrino_PU200/results/ntuple_SingleNeutrino_PU200_%i.root'
  tree = TreeChain('ntupler/tree', [(tree_name_i % i) for i in range(1,51)])
else:
  raise Exception("Unexpected state: %i" % mystate)
maxEvents = -1
#maxEvents = 2000

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM, kTT = 0, 1, 2, 3, 20

# Lambdas
deg_to_rad = lambda x: x * np.pi/180.

rad_to_deg = lambda x: x * 180./np.pi

# Functions
def delta_phi(lhs, rhs):  # in degrees
  deg = lhs - rhs
  while deg <  -180.:  deg += 360.
  while deg >= +180.:  deg -= 360.
  return deg

# Constants
INFO = '\033[1;37mINFO\033[0m]'

#pt_bins = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 28., 30., 34., 40., 48., 60., 80., 100., 140., 200., 300., 500., 1000.]
#eta_bins = [0., 0.83, 1.24, 1.64, 2.14, 2.5]

pt_bins = list(reversed([100./x for x in xrange(1,51)])) + [250.0]
eta_bins = [0.0, 0.8] + [1.2 + x * 0.2 for x in xrange(7)]

# Classes
class EfficiencyMatrix:
  def __init__(self, nbinsx=24, xmin=0., xmax=2.4, xbins=None, nbinsy=100, ymin=0., ymax=100., ybins=None, nbinsz=100, zmin=0., zmax=100., zbins=None):
    """
    x: gen_eta
    y: gen_pt
    z: l1t_pt
    """
    if xbins:  nbinsx = len(xbins)
    if ybins:  nbinsy = len(ybins)
    if zbins:  nbinsz = len(zbins)
    self._nbinsx = nbinsx
    self._xmin = xmin
    self._xmax = xmax
    self._xwidth = (xmax - xmin) / float(nbinsx)
    self._xbins = np.array(xbins, dtype=float) if xbins else None
    self._nbinsy = nbinsy
    self._ymin = ymin
    self._ymax = ymax
    self._ywidth = (ymax - ymin) / float(nbinsy)
    self._ybins = np.array(ybins, dtype=float) if ybins else None
    self._nbinsz = nbinsz
    self._zmin = zmin
    self._zmax = zmax
    self._zwidth = (zmax - zmin) / float(nbinsz)
    self._zbins = np.array(zbins, dtype=float) if zbins else None
    self._numer = np.zeros((nbinsx, nbinsy, nbinsz), dtype=int)
    self._denom = np.zeros((nbinsx, nbinsy, nbinsz), dtype=int)
    self._effie = np.zeros((nbinsx, nbinsy, nbinsz), dtype=float)

    # In addition, profile phi and eta
    self._phi_cnt = np.zeros((nbinsx, nbinsy), dtype=int)
    self._phi_mean = np.zeros((nbinsx, nbinsy), dtype=float)
    self._phi_var = np.zeros((nbinsx, nbinsy), dtype=float)
    self._eta_cnt = np.zeros((nbinsx, nbinsy), dtype=int)
    self._eta_mean = np.zeros((nbinsx, nbinsy), dtype=float)
    self._eta_var = np.zeros((nbinsx, nbinsy), dtype=float)

  def sanity_check(self):
    def is_sorted(a):
      if a is None:
        return True
      else:
        return np.all(np.less_equal(a[:-1], a[1:]))

    if self._xbins is not None:
      assert(is_sorted(self._xbins))
    if self._ybins is not None:
      assert(is_sorted(self._ybins))
    if self._zbins is not None:
      assert(is_sorted(self._zbins))

  def find_binx(self, x):
    if self._xbins is not None:
      binx = np.searchsorted(self._xbins, x)
      binx -= 1
    else:
      binx = int(np.floor((x - self._xmin) / self._xwidth))
    if binx < 0: binx = 0
    if binx >= self._nbinsx: binx = self._nbinsx - 1
    return binx

  def find_biny(self, y):
    if self._ybins is not None:
      biny = np.searchsorted(self._ybins, y)
      biny -= 1
    else:
      biny = int(np.floor((y - self._ymin) / self._ywidth))
    if biny < 0: biny = 0
    if biny >= self._nbinsy: biny = self._nbinsy - 1
    return biny

  def find_binz(self, z):
    if self._zbins is not None:
      binz = np.searchsorted(self._zbins, z)
      binz -= 1
    else:
      binz = int(np.floor((z - self._zmin) / self._zwidth))
    if binz < 0: binz = 0
    if binz >= self._nbinsz: binz = self._nbinsz - 1
    return binz

  def find_edgex(self, binx):
    if self._xbins is not None:
      edgex = self._xbins[binx]
    else:
      edgex = self._xmin + binx * self._xwidth
    return edgex

  def find_edgey(self, biny):
    if self._ybins is not None:
      edgey = self._ybins[biny]
    else:
      edgey = self._ymin + biny * self._ywidth
    return edgey

  def find_edgez(self, binz):
    if self._zbins is not None:
      edgez = self._zbins[binz]
    else:
      edgez = self._zmin + binz * self._zwidth
    return edgez

  def fill(self, gen_eta, gen_pt, l1t_pt):
    binx = self.find_binx(gen_eta)
    biny = self.find_biny(gen_pt)
    binz = self.find_binz(l1t_pt)
    #
    passed = [z <= binz for z in xrange(self._nbinsz)]
    self._denom[binx, biny, :] += 1
    self._numer[binx, biny, passed] += 1

  def profile(self, gen_eta, gen_pt, phi, eta):
    binx = self.find_binx(gen_eta)
    biny = self.find_biny(gen_pt)
    #
    self._phi_cnt[binx, biny] += 1
    delta = (phi - self._phi_mean[binx, biny])
    self._phi_mean[binx, biny] += delta / float(self._phi_cnt[binx, biny])
    self._phi_var[binx, biny] += delta * (phi - self._phi_mean[binx, biny])
    #
    self._eta_cnt[binx, biny] += 1
    delta = (eta - self._eta_mean[binx, biny])
    self._eta_mean[binx, biny] += delta / float(self._eta_cnt[binx, biny])
    self._eta_var[binx, biny] += delta * (eta - self._eta_mean[binx, biny])

  def freeze(self):
    tmp_numer = self._numer.astype(float)
    tmp_denom = self._denom.astype(float)
    tmp_denom = np.where(tmp_denom < 1e-9, 1e-9, tmp_denom)  # avoid division by zero
    np.true_divide(tmp_numer, tmp_denom, out=self._effie)
    #a = self._effie.reshape(-1)
    #for i, v in enumerate(a):
    #  if np.isnan(v):
    #    a_[i] = 0.
    #
    tmp_numer = self._phi_var.astype(float)
    tmp_denom = self._phi_cnt.astype(float)
    tmp_denom = np.where(tmp_denom < 1e-9, 1e-9, tmp_denom)  # avoid division by zero
    np.true_divide(tmp_numer, tmp_denom, out=self._phi_var)
    #
    tmp_numer = self._eta_var.astype(float)
    tmp_denom = self._eta_cnt.astype(float)
    tmp_denom = np.where(tmp_denom < 1e-9, 1e-9, tmp_denom)  # avoid division by zero
    np.true_divide(tmp_numer, tmp_denom, out=self._eta_var)

  def sitrep(self):
    print self._effie


# Book histograms
histograms = {}
histogram2Ds = {}

hname = "nevents"
histograms[hname] = Hist(5, 0, 5, name=hname, title="; count", type='F')

hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

for i in xrange(14,22+1):
  hname = "emtf_ptmin%i_qmin12_eta" % i
  histograms[hname] = Hist(10, 1.55, 2.55, name=hname, title="; |#eta|; entries", type='F')


hname = "tp_emtf_absEtaMin0_absEtaMax2.5_genpt"  # efficiency denom
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

for i in xrange(14,22+1):
  hname = "tp_emtf_absEtaMin0_absEtaMax2.5_qmin12_ptmin%i_genpt" % i  # efficiency numer
  histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

hname = "tp_emtf_absEtaMin1.65_absEtaMax2.15_genpt"  # efficiency denom
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

for i in xrange(14,22+1):
  hname = "tp_emtf_absEtaMin1.65_absEtaMax2.15_qmin12_ptmin%i_genpt" % i  # efficiency numer
  histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')


em = EfficiencyMatrix(xbins=eta_bins, ybins=pt_bins)
em.sanity_check()

print INFO, "eta_bins=%s" % repr(eta_bins), "pt_bins=%s" % repr(pt_bins)

# ______________________________________________________________________________
# Load efficiency matrix

if mystate == 1:
  a = np.loadtxt('test.out', delimiter=',')
  a = a.reshape(em._effie.shape)
  assert(em._effie.shape == a.shape)
  em._effie = a
  #
  a = np.loadtxt('test1.out', delimiter=',')
  assert(em._phi_mean.shape == a.shape)
  em._phi_mean = a
  #
  a = np.loadtxt('test2.out', delimiter=',')
  assert(em._eta_mean.shape == a.shape)
  em._eta_mean = a


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
  # Fill efficiency matrix
  if mystate == 0:
    assert len(evt.genparticles) == 1

    mypart = evt.genparticles[0]
    #mytrk = evt.tracks[0]
    try:
      # Select highest pT track from tracks that have a station 1 hit
      mytrk = max(filter(lambda x: (x.mode in [11,13,14,15]), evt.tracks), key=lambda x: x.pt)
    except ValueError, e:
      mytrk = None

    # Fill efficiency matrix
    gen_eta = mypart.eta
    gen_pt = mypart.pt
    if mytrk:
      l1t_pt = mytrk.pt
      l1t_phi = delta_phi(mytrk.phi, rad_to_deg(mypart.phi)) if mytrk.q < 0 else -delta_phi(mytrk.phi, rad_to_deg(mypart.phi))
      l1t_eta = mytrk.eta
    else:
      l1t_pt = 0.
      l1t_phi = 0.
      l1t_eta = 0.

    em.fill(abs(gen_eta), gen_pt, l1t_pt)
    if mytrk:
      em.profile(abs(gen_eta), gen_pt, l1t_phi, abs(l1t_eta))

    def doit():
      h = histograms[hname]
      if select(mytrk):
        h.fill(gen_pt)

    if (0. <= abs(gen_eta) <= 2.5):
      select = lambda trk: True
      hname = "tp_emtf_absEtaMin0_absEtaMax2.5_genpt"
      doit()

      for select_pt in xrange(14, 22+1):
        select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15]) and (trk.pt > float(select_pt))
        hname = "tp_emtf_absEtaMin0_absEtaMax2.5_qmin12_ptmin%i_genpt" % select_pt
        doit()

    if (1.65 <= abs(gen_eta) <= 2.15):
      select = lambda trk: True
      hname = "tp_emtf_absEtaMin1.65_absEtaMax2.15_genpt"
      doit()

      for select_pt in xrange(14, 22+1):
        select = lambda trk: trk and (1.65 <= abs(trk.eta) <= 2.15) and (trk.mode in [11,13,14,15]) and (trk.pt > float(select_pt))
        hname = "tp_emtf_absEtaMin1.65_absEtaMax2.15_qmin12_ptmin%i_genpt" % select_pt
        doit()


  # ____________________________________________________________________________
  # Fill efficiency matrix
  elif mystate == 1:

    #for ipart, part in enumerate(evt.genparticles):
    #  pass

    h = histograms["nevents"]
    h.fill(1.0)

    def doit():
      h = histograms[hname]
      highest_pt = -999999.
      for itrk, trk in enumerate(evt.tracks):
        if select(trk):
          if highest_pt < trk.pt:
            highest_pt = trk.pt
      if highest_pt > 0.:
        highest_pt = min(100.-1e-3, highest_pt)
        h.fill(highest_pt)

    select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15])
    hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"
    doit()

    def doit():
      h = histograms[hname]
      eta_bins = [False] * (10+2)
      for itrk, trk in enumerate(evt.tracks):
        if select(trk):
          b = h.FindFixBin(abs(trk.eta))
          eta_bins[b] = True
      for b in xrange(len(eta_bins)):
        if eta_bins[b]:
          h.fill(h.GetBinCenter(b))

    for select_pt in xrange(14, 22+1):
      select = lambda trk: trk and (trk.pt > float(select_pt)) and (trk.mode in [11,13,14,15])
      hname = "emtf_ptmin%i_qmin12_eta" % select_pt
      doit()

  continue  # end loop over event

# ______________________________________________________________________________
# Save efficiency matrix

if mystate == 0:
  em.freeze()
  em.sitrep()
  #np.savetxt('test.out', em._effie, delimiter=',')  # only works for 2D array
  np.savetxt('test.out', em._effie.reshape(-1,em._effie.shape[-1]), delimiter=',')  # reshape to 2D array
  #
  try:
    a = np.loadtxt('test.out', delimiter=',')
    a.reshape(em._effie.shape)
  except:
    raise
  #
  np.savetxt('test1.out', em._phi_mean, delimiter=',')
  np.savetxt('test2.out', em._eta_mean, delimiter=',')

# ______________________________________________________________________________
# Write histograms

tag = '{0:%y}{0:%m}{0:%d}_{0:%H}{0:%M}{0:%S}'.format(datetime.now())
with root_open('rateplots_rootpy_%s.root' % tag, 'RECREATE') as f:
  f.mkdir('trackcounting').cd()
  for k, v in histograms.iteritems():
    v.Write()
  for k, v in histogram2Ds.iteritems():
    v.Write()
