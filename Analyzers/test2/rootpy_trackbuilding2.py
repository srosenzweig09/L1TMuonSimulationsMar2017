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
#infile = root_open('rateplots_mc_r305310_run2_all.0.root')
infile = root_open('rateplots_mc_r305310_run2_all.1.root')

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

for cat in ["real", "fake"]:
  hname = "trk_eta_%s" % cat
  histograms[hname] = Hist(30, 1.0, 2.5, name=hname, title="; trk |#eta|", type='F')
  hname = "trk_mode_%s" % cat
  histograms[hname] = Hist(16, 0, 16, name=hname, title="; trk mode", type='F')

  hname = "trk_hit_type1_%s" % cat
  histograms[hname] = Hist(7, -2, 5, name=hname, title="; hit type (ME1)", type='F')
  hname = "trk_hit_type2_%s" % cat
  histograms[hname] = Hist(7, -2, 5, name=hname, title="; hit type (ME2)", type='F')
  hname = "trk_hit_type3_%s" % cat
  histograms[hname] = Hist(7, -2, 5, name=hname, title="; hit type (ME3)", type='F')
  hname = "trk_hit_type4_%s" % cat
  histograms[hname] = Hist(7, -2, 5, name=hname, title="; hit type (ME4)", type='F')

  hname = "trk_hit_pattern1_%s" % cat
  histograms[hname] = Hist(12, 0, 12, name=hname, title="; hit pattern (ME1)", type='F')
  hname = "trk_hit_pattern2_%s" % cat
  histograms[hname] = Hist(12, 0, 12, name=hname, title="; hit pattern (ME2)", type='F')
  hname = "trk_hit_pattern3_%s" % cat
  histograms[hname] = Hist(12, 0, 12, name=hname, title="; hit pattern (ME3)", type='F')
  hname = "trk_hit_pattern4_%s" % cat
  histograms[hname] = Hist(12, 0, 12, name=hname, title="; hit pattern (ME4)", type='F')

  hname = "trk_hit_occu1_%s" % cat
  histograms[hname] = Hist(30, 0, 30, name=hname, title="; hit occupancy in sector (ME1)", type='F')
  hname = "trk_hit_occu2_%s" % cat
  histograms[hname] = Hist(30, 0, 30, name=hname, title="; hit occupancy in sector (ME2)", type='F')
  hname = "trk_hit_occu3_%s" % cat
  histograms[hname] = Hist(30, 0, 30, name=hname, title="; hit occupancy in sector (ME3)", type='F')
  hname = "trk_hit_occu4_%s" % cat
  histograms[hname] = Hist(30, 0, 30, name=hname, title="; hit occupancy in sector (ME4)", type='F')

  hname = "trk_dphi1_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{12}", type='F')
  hname = "trk_dphi2_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{23}", type='F')
  hname = "trk_dphi3_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{34}", type='F')
  hname = "trk_dphi4_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{13}", type='F')
  hname = "trk_dphi5_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{14}", type='F')
  hname = "trk_dphi6_%s" % cat
  histograms[hname] = Hist(61, -61, 61, name=hname, title="; #Delta#phi_{24}", type='F')

  hname = "trk_dtheta1_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{12}", type='F')
  hname = "trk_dtheta2_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{23}", type='F')
  hname = "trk_dtheta3_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{34}", type='F')
  hname = "trk_dtheta4_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{13}", type='F')
  hname = "trk_dtheta5_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{14}", type='F')
  hname = "trk_dtheta6_%s" % cat
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; #Delta#theta_{24}", type='F')

  hname = "trk_bad_hit_%s" % cat
  histograms[hname] = Hist(7, 0, 7, name=hname, title="; station with bad hit", type='F')


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

      for st in xrange(4):
        for ihit, hit in enumerate(evt.hits):
          if hit.sim_tp1 == ipart and hit.sim_tp2 == ipart and hit.neighbor == False:
            if hit.station == st + 1:
              print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.emtf_phi, hit.emtf_theta, hit.sim_tp1, hit.sim_tp2))
            #elif hit.station == st + 1 and hit.type == kRPC:
            #  print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.emtf_phi, hit.emtf_theta, hit.sim_tp1, hit.sim_tp2))

    for itrk, trk in enumerate(evt.tracks):
      print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7} {8}".format(itrk, trk.pt, trk.eta, trk.sector, trk.mode, trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4))
      for ihit in [trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4]:
        if ihit != -1:
          hit = evt.hits[ihit]
          print(".... hit {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.strip, hit.wire, hit.emtf_phi, hit.emtf_theta, hit.sim_tp1, hit.sim_tp2))


  # ____________________________________________________________________________
  # Full study

  full_study = True

  if full_study:

    # Number of events
    histograms["nevents"].fill(1.0)

    if False:
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


    if True:
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
            hit_sim_tp1s.append(-2)
            hit_sim_tp2s.append(-2)

        trk.mctruth = 0
        #hit_sim_tp_set = set(tp for tp in (hit_sim_tp1s + hit_sim_tp2s) if tp != -1)
        hit_sim_tp_set = set(tp for tp in (hit_sim_tp1s + hit_sim_tp2s) if tp != -2)
        assert(-2 not in hit_sim_tp_set)

        if len(hit_sim_tp_set) == 1 and (-1 not in hit_sim_tp_set):  # match to the same tracking particle
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

      # Real vs Fake
      for itrk, trk in enumerate(evt.tracks):
        if (0. <= abs(trk.eta) <= 2.5) and (trk.mode in [11,13,14,15]) and (trk.pt >= 16.):
          if trk.mctruth >= 2:
            cat = "real"
          else:
            cat = "fake"

          def get_hit_occupancy(endcap, sector, station):
            cnt = 0
            for ihit, hit in enumerate(evt.hits):
              if hit.endcap == endcap and hit.sector == sector and hit.station == station and hit.type == kCSC:  # only count CSC
                cnt += 1
            return cnt

          hit1 = evt.hits[trk.hitref1] if trk.hitref1 != -1 else None
          hit2 = evt.hits[trk.hitref2] if trk.hitref2 != -1 else None
          hit3 = evt.hits[trk.hitref3] if trk.hitref3 != -1 else None
          hit4 = evt.hits[trk.hitref4] if trk.hitref4 != -1 else None

          hit_type1 = hit1.type if hit1 else -1
          hit_type2 = hit2.type if hit2 else -1
          hit_type3 = hit3.type if hit3 else -1
          hit_type4 = hit4.type if hit4 else -1

          hit_pattern1 = hit1.pattern if hit1 and hit1.type == kCSC else -1
          hit_pattern2 = hit2.pattern if hit2 and hit2.type == kCSC else -1
          hit_pattern3 = hit3.pattern if hit3 and hit3.type == kCSC else -1
          hit_pattern4 = hit4.pattern if hit4 and hit4.type == kCSC else -1

          hit_occu1 = get_hit_occupancy(trk.endcap, trk.sector, 1)
          hit_occu2 = get_hit_occupancy(trk.endcap, trk.sector, 2)
          hit_occu3 = get_hit_occupancy(trk.endcap, trk.sector, 3)
          hit_occu4 = get_hit_occupancy(trk.endcap, trk.sector, 4)

          dphi1 = hit2.emtf_phi - hit1.emtf_phi if hit2 and hit1 else -999
          dphi2 = hit3.emtf_phi - hit2.emtf_phi if hit3 and hit2 else -999
          dphi3 = hit4.emtf_phi - hit3.emtf_phi if hit4 and hit3 else -999
          dphi4 = hit3.emtf_phi - hit1.emtf_phi if hit3 and hit1 else -999
          dphi5 = hit4.emtf_phi - hit1.emtf_phi if hit4 and hit1 else -999
          dphi6 = hit4.emtf_phi - hit2.emtf_phi if hit4 and hit2 else -999

          dtheta1 = hit2.emtf_theta - hit1.emtf_theta if hit2 and hit1 else -999
          dtheta2 = hit3.emtf_theta - hit2.emtf_theta if hit3 and hit2 else -999
          dtheta3 = hit4.emtf_theta - hit3.emtf_theta if hit4 and hit3 else -999
          dtheta4 = hit3.emtf_theta - hit1.emtf_theta if hit3 and hit1 else -999
          dtheta5 = hit4.emtf_theta - hit1.emtf_theta if hit4 and hit1 else -999
          dtheta6 = hit4.emtf_theta - hit2.emtf_theta if hit4 and hit2 else -999

          # Fill
          hname = "trk_eta_%s" % cat
          histograms[hname].fill(abs(trk.eta))
          hname = "trk_mode_%s" % cat
          histograms[hname].fill(trk.mode)

          hname = "trk_hit_type1_%s" % cat
          histograms[hname].fill(hit_type1)
          hname = "trk_hit_type2_%s" % cat
          histograms[hname].fill(hit_type2)
          hname = "trk_hit_type3_%s" % cat
          histograms[hname].fill(hit_type3)
          hname = "trk_hit_type4_%s" % cat
          histograms[hname].fill(hit_type4)

          hname = "trk_hit_pattern1_%s" % cat
          histograms[hname].fill(hit_pattern1)
          hname = "trk_hit_pattern2_%s" % cat
          histograms[hname].fill(hit_pattern2)
          hname = "trk_hit_pattern3_%s" % cat
          histograms[hname].fill(hit_pattern3)
          hname = "trk_hit_pattern4_%s" % cat
          histograms[hname].fill(hit_pattern4)

          hname = "trk_hit_occu1_%s" % cat
          histograms[hname].fill(hit_occu1)
          hname = "trk_hit_occu2_%s" % cat
          histograms[hname].fill(hit_occu2)
          hname = "trk_hit_occu3_%s" % cat
          histograms[hname].fill(hit_occu3)
          hname = "trk_hit_occu4_%s" % cat
          histograms[hname].fill(hit_occu4)

          hname = "trk_dphi1_%s" % cat
          histograms[hname].fill(dphi1)
          hname = "trk_dphi2_%s" % cat
          histograms[hname].fill(dphi2)
          hname = "trk_dphi3_%s" % cat
          histograms[hname].fill(dphi3)
          hname = "trk_dphi4_%s" % cat
          histograms[hname].fill(dphi4)
          hname = "trk_dphi5_%s" % cat
          histograms[hname].fill(dphi5)
          hname = "trk_dphi6_%s" % cat
          histograms[hname].fill(dphi6)

          hname = "trk_dtheta1_%s" % cat
          histograms[hname].fill(dtheta1)
          hname = "trk_dtheta2_%s" % cat
          histograms[hname].fill(dtheta2)
          hname = "trk_dtheta3_%s" % cat
          histograms[hname].fill(dtheta3)
          hname = "trk_dtheta4_%s" % cat
          histograms[hname].fill(dtheta4)
          hname = "trk_dtheta5_%s" % cat
          histograms[hname].fill(dtheta5)
          hname = "trk_dtheta6_%s" % cat
          histograms[hname].fill(dtheta6)

          # Find station with bad hits
          if cat == "fake":
            for st in xrange(4):
              hit_refs_minus = [trk.hitref1, trk.hitref2, trk.hitref3, trk.hitref4]
              del hit_refs_minus[st]
              #
              hit_sim_tp1s = []
              hit_sim_tp2s = []
              for hit_ref in hit_refs_minus:
                if hit_ref != -1:
                  hit_sim_tp1s.append(evt.hits[hit_ref].sim_tp1)
                  hit_sim_tp2s.append(evt.hits[hit_ref].sim_tp2)
                else:
                  hit_sim_tp1s.append(-2)
                  hit_sim_tp2s.append(-2)
              #
              hit_sim_tp_set = set(tp for tp in (hit_sim_tp1s + hit_sim_tp2s) if tp != -2)
              assert(-2 not in hit_sim_tp_set)
              if len(hit_sim_tp_set) == 1:
                hname = "trk_bad_hit_%s" % cat
                histograms[hname].fill(st + 1)


# ______________________________________________________________________________
# Write histograms

with root_open('histos_tb.root', 'RECREATE') as f:
  f.mkdir('trackcounting').cd()
  for k, v in histograms.iteritems():
    v.Write()
  for k, v in histogram2Ds.iteritems():
    v.Write()
