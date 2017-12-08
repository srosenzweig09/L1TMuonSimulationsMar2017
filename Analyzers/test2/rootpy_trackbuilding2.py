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
#infile = root_open('ntuple_SingleMuon_PositiveEndCap.6.root')
infile = root_open('rateplots_mc_r305310_run2_all.0.root')

tree = infile.ntupler.tree
maxEvents = -1
#maxEvents = 100
print "[INFO] Opening file: %s" % infile

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='particles', prefix='vp_', size='vp_size')

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


# Book histograms
histograms = {}
histogram2Ds = {}

for st in xrange(4):
  hname = "hit_occupancy_csc_st%i" % (st+1)
  histograms[hname] = Hist(21, -0.5, 20.5, name=hname, title="; # of hits per sector (ME%i)" % (st+1), type='F')

  hname = "hit_occupancy_rpc_st%i" % (st+1)
  histograms[hname] = Hist(21, -0.5, 20.5, name=hname, title="; # of hits per sector (RE%i)" % (st+1), type='F')

hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth1_pt"
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth2_pt"
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth3_pt"
histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

hname = "nevents"
histograms[hname] = Hist(5, 0, 5, name=hname, title="; count", type='F')




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
  # Quick study
  # - Single muon MC matching

  quick_study = False

  if quick_study:

    ntrks_pt10 = 0
    for itrk, trk in enumerate(evt.tracks):
      if trk.mode in [11,13,14,15] and trk.pt > 10.:
        ntrks_pt10 += 1
    print("evt {0} has {1} hits and {2} tracks ({3} with pT > 10)".format(ievt, len(evt.hits), len(evt.tracks), ntrks_pt10))

    for ipart, part in enumerate(evt.particles):
      print(".. part {0} {1} {2}".format(ipart, part.pt, part.eta))

      for istation in xrange(4):
        for ihit, hit in enumerate(evt.hits):
          if hit.sim_tp1 == ipart and hit.sim_tp2 == ipart and hit.neighbor == False:
            if hit.station == istation + 1 and hit.type == kCSC:
              print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.sim_tp1, hit.sim_tp2))
            #elif hit.station == istation + 1 and hit.type == kRPC:
            #  print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.sim_tp1, hit.sim_tp2))

    for itrk, trk in enumerate(evt.tracks):
      print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7} {8}".format(itrk, trk.pt, trk.eta, trk.sector, trk.mode, trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4))
      for ihit in [trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4]:
        if ihit != -1:
          hit = evt.hits[ihit]
          print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.sim_tp1, hit.sim_tp2))


  # ____________________________________________________________________________
  # Full study

  full_study = True

  if full_study:

    # Number of events
    histograms["nevents"].fill(1.0)

    # Get hit occupancy
    def doit():
      h = histograms[hname]
      for endcap in [-1,+1]:
        for sector in [1,2,3,4,5,6]:
          cnt = 0
          for ihit, hit in enumerate(evt.hits):
            if hit.endcap == endcap and hit.sector == sector and select(hit):
              cnt += 1
          h.fill(cnt)

    for st in xrange(4):
      select = lambda hit: hit and (hit.type == kCSC) and (hit.station == st+1)
      hname = "hit_occupancy_csc_st%i" % (st+1)
      doit()

      select = lambda hit: hit and (hit.type == kRPC) and (hit.station == st+1)
      hname = "hit_occupancy_rpc_st%i" % (st+1)
      doit()


    # Get track rates
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

    if True:
      # MC truth matching
      for itrk, trk in enumerate(evt.tracks):
        hit_refs = [trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4]
        hit_sim_tp1s = []
        hit_sim_tp2s = []
        for hit_ref in hit_refs:
          if hit_ref != -1:
            hit_sim_tp1s.append(evt.hits[hit_ref].sim_tp1)
            hit_sim_tp2s.append(evt.hits[hit_ref].sim_tp2)
          else:
            hit_sim_tp1s.append(-1)
            hit_sim_tp2s.append(-1)

        trk.mctruth = 0
        hit_sim_tp_set = set(tp for tp in (hit_sim_tp1s + hit_sim_tp2s) if tp != -1)

        if len(hit_sim_tp_set) == 1:  # match to the same tracking particle
          trk.mctruth = 1
          tp = next(iter(hit_sim_tp_set))
          part = evt.particles[tp]
          if abs(trk.xml_pt - part.pt) / part.pt < 0.4:  # within 40% pT resolution
            trk.mctruth = 2
          if abs(trk.xml_pt - part.pt) / part.pt < 0.2:  # within 20% pT resolution
            trk.mctruth = 3

      select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15]) and (trk.mctruth >= 1)
      hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth1_pt"
      doit()

      select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15]) and (trk.mctruth >= 2)
      hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth2_pt"
      doit()

      select = lambda trk: trk and (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15]) and (trk.mctruth >= 3)
      hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_mctruth3_pt"
      doit()


# ______________________________________________________________________________
# Write histograms

with root_open('histos_tb.root', 'RECREATE') as f:
  f.mkdir('trackcounting').cd()
  for k, v in histograms.iteritems():
    v.Write()
  for k, v in histogram2Ds.iteritems():
    v.Write()
