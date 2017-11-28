import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, TreeChain, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from ROOT import gROOT, gStyle, gPad


# ______________________________________________________________________________
# Tree models
#   see: http://www.rootpy.org/auto_examples/tree/model.html

class Hit(TreeModel):
  pass

class Track(TreeModel):
  pass

class Particle(TreeModel):
  q = IntCol()
  pt = FloatCol()
  eta = FloatCol()
  phi = FloatCol()


# ______________________________________________________________________________
# Analyzer

# Open file
#infile = root_open('rateplots_mc.0.root')
infile = root_open('rateplots_mc_r305310_run2_all.root')
tree = infile.ntupler.tree
maxEvents = -1
#maxEvents = 50000
print "[INFO] Opening file: %s" % infile

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM, kTT = 0, 1, 2, 3, 20

# Lambdas
get_zone_hit = lambda x: ((x + (1<<4)) >> 5)

extrapolate_to_muon = lambda phi, qoverpt: phi - 1.204 * qoverpt  # phi in radians

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

def calc_eta_from_theta_rad(theta_rad):
  return -1. * np.log(np.tan(theta_rad/2.))

def calc_eta_from_theta_deg(theta_deg):
  theta_rad = np.deg2rad(theta_deg)
  return calc_eta_from_theta_rad(theta_rad)

def select_by_eta(eta):
  return 1.24 <= abs(eta) < 2.4

def select_by_bx(bx):
  return bx == 0

def select_by_vertex(vx, vy, vz):
  return np.sqrt(vx*vx + vy*vy) < 15. and abs(vz) < 50.


def find_zones(eta):
  zones = []
  endcap = 1 if eta >= 0. else -1
  theta_deg = calc_theta_deg_from_eta(eta)
  theta_int = calc_theta_int(theta_deg, endcap)
  boundaries = (0,41,49,87,127)
  overlap = 2
  if boundaries[0] <= theta_int <= boundaries[1]+overlap:
    zones.append(0)
  if boundaries[1]-overlap <= theta_int <= boundaries[2]+overlap:
    zones.append(1)
  if boundaries[2]-overlap <= theta_int <= boundaries[3]+overlap:
    zones.append(2)
  if boundaries[3]-overlap <= theta_int <= boundaries[4]:
    zones.append(3)
  assert(len(zones) <= 2)
  return zones


# Book histograms
histograms = {}
histogram2Ds = {}

for i in xrange(4):
  hname = "h_dphi_re%i" % (i+1)
  histograms[hname] = Hist(41, -20.5, 20.5, name=hname, title="; diff in phi_int", type="F")
  hname = "h_dtheta_re%i" % (i+1)
  histograms[hname] = Hist(21, -10.5, 10.5, name=hname, title="; diff in theta_int", type="F")

hname = "muon_ptmin2_dxy"
histograms[hname] = Hist(200, 0, 100, name=hname, title="; |d_{xy}| [cm]", type="F")
hname = "muon_ptmin2_dz"
histograms[hname] = Hist(200, 0, 100, name=hname, title="; |d_{z}| [cm]", type="F")
hname = "muon_ptmin2_eta"
histograms[hname] = Hist(50, 0, 2.5, name=hname, title="; |#eta|", type="F")
hname = "muon_ptmin2_bx"
histograms[hname] = Hist(11, -5.5, 5.5, name=hname, title="; BX", type="F")

for z in xrange(5):
  # zone 0-3 are the 4 zones. zone 4 is inclusive.
  hname = "muon_pt_denom_zone%i" % z
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; p_{T} [GeV]", type="F")
  hname = "muon_pt_fiducial_zone%i" % z
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; p_{T} [GeV]", type="F")
  hname = "muon_pt_trigger_zone%i" % z
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; p_{T} [GeV]", type="F")


# Tools
pattern_bank = []


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
    for ipart, part in enumerate(evt.genparticles):
      print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))
  else:
    if (ievt % 1000 == 0):  print("Processing event: {0}".format(ievt))


  # ____________________________________________________________________________
  # Quick study (1)
  # - Find the hit patterns, using EMTF tracks to flag an event

  quick_study = False

  if quick_study:
    print("evt {0} has {1} hits and {2} tracks".format(ievt, len(evt.hits), len(evt.tracks)))

    trigger = False
    for itrk, trk in enumerate(evt.tracks):
      if trk.pt > 7.:
        print(".. trk  {0} {1} {2} {3} {4}".format(itrk, trk.mode, trk.pt, trk.eta, np.deg2rad(trk.phi)))
        trigger = True
        break

    if trigger:
      for ihit, hit in enumerate(evt.hits):
        if hit.type == kCSC:
          part = Particle()
          if hit.sim_tp1 != -1:
            part = evt.genparticles[hit.sim_tp1]

          print(".. hit  {0} {1} {2} {3} {4} {5} -> part {6} {7} {8} {9}".format(ihit, hit.sector, hit.station, get_zone_hit(hit.emtf_phi), hit.emtf_theta, hit.sim_tp1, part.pt, part.eta, part.phi, extrapolate_to_muon(part.phi, np.true_divide(part.q, part.pt))))


  # ____________________________________________________________________________
  # Quick study (2)
  # - Find the hit patterns, using tracking particles to flag an event

  quick_study2 = True

  if quick_study2:
    print("evt {0} has {1} hits and {2} tracks".format(ievt, len(evt.hits), len(evt.tracks)))

    trigger = False
    for ipart, part in enumerate(evt.genparticles):
      if 1.24 <= abs(part.eta) < 2.5 and part.pt > 7.:
        print(".. part {0} {1} {2} {3} {4} {5} {6}".format(ipart, part.pt, part.eta, part.phi, part.vx, part.vy, part.vz))
        trigger = True
        break

    if trigger:
      for ihit, hit in enumerate(evt.hits):
        if hit.type == kCSC:
          part = Particle()
          if hit.sim_tp1 != -1:
            part = evt.genparticles[hit.sim_tp1]

          print(".. hit  {0} {1} {2} {3} {4} {5} -> part {6} {7} {8} {9}".format(ihit, hit.sector, hit.station, get_zone_hit(hit.emtf_phi), hit.emtf_theta, hit.sim_tp1, part.pt, part.eta, part.phi, extrapolate_to_muon(part.phi, np.true_divide(part.q, part.pt))))


  # ____________________________________________________________________________
  # Quick study (3)
  # - Check RPC hit conversion

  quick_study3 = False

  if quick_study3:
    print("evt {0} has {1} hits and {2} tracks".format(ievt, len(evt.hits), len(evt.tracks)))

    for ihit, hit in enumerate(evt.hits):
      if hit.type == kRPC:
        print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(ihit, hit.sector, hit.station, hit.ring, hit.emtf_phi, hit.emtf_theta, hit.sim_phi, hit.sim_theta, calc_phi_loc_int(hit.sim_phi, hit.sector), calc_theta_int(hit.sim_theta, hit.endcap)))

        st = hit.station
        hname = "h_dphi_re%i" % st
        histograms[hname].fill(hit.emtf_phi - calc_phi_loc_int(hit.sim_phi, hit.sector))
        hname = "h_dtheta_re%i" % st
        histograms[hname].fill(hit.emtf_theta - calc_theta_int(hit.sim_theta, hit.endcap))


  # ____________________________________________________________________________
  # Quick study (4)
  # - Check quadstrip occupancy

  quick_study4 = False

  if quick_study4:
    print("evt {0} has {1} hits and {2} tracks".format(ievt, len(evt.hits), len(evt.tracks)))



  # ____________________________________________________________________________
  # Full study

  full_study = True

  if full_study:

    nhits_csc = 0
    for ihit, hit in enumerate(evt.hits):
      if hit.type == kCSC:
        nhits_csc += 1
    print("evt {0} has {1} hits ({2} CSC) and {3} tracks".format(ievt, len(evt.hits), nhits_csc, len(evt.tracks)))

    # Make list of particles
    map_of_iparts = {}

    for ihit, hit in enumerate(evt.hits):
      if hit.type == kCSC:
        if hit.sim_tp1 != -1:
          if hit.sim_tp1 not in map_of_iparts:
            map_of_iparts[hit.sim_tp1] = 0
          istation = (hit.station - 1)
          map_of_iparts[hit.sim_tp1] |= (1<<istation)

    map_of_iparts_1 = {}

    for k, v in map_of_iparts.iteritems():
      cnt = 0
      if (v & (1<<0)):  # station 1
        cnt += 1
      if (v & (1<<1)):  # station 2
        cnt += 1
      if (v & (1<<2)) or (v & (1<<3)):  # station 3 or 4
        cnt += 1
      if cnt >= 2:  # at least one hit in two different layers
        map_of_iparts_1[k] = v

    # Get part
    highest_pt = 0.
    highest_pt_ipart = -1
    for ipart in map_of_iparts_1:
      part = evt.genparticles[ipart]
      if part.pt > 2.:  histograms["muon_ptmin2_eta"].Fill(part.eta)

      if select_by_eta(part.eta):
        if part.pt > 2.:  histograms["muon_ptmin2_bx"].Fill(part.bx)

        if select_by_bx(part.bx):
          if part.pt > 2.:  histograms["muon_ptmin2_dxy"].Fill(np.sqrt(part.vx*part.vx + part.vy*part.vy))
          if part.pt > 2.:  histograms["muon_ptmin2_dz"].Fill(abs(part.vz))

          if select_by_vertex(part.vx, part.vy, part.vz):
            if highest_pt < part.pt:
              highest_pt = part.pt
              highest_pt_ipart = ipart
      continue  # end for


    # Quick efficiency study
    for ipart, part in enumerate(evt.genparticles):
      if 1.2 <= abs(part.eta) < 2.2:  # restricted eta to avoid geometric acceptance effect
        if select_by_vertex(part.vx, part.vy, part.vz) and select_by_bx(part.bx):
          zones = find_zones(part.eta)
          for z in zones + [4]:
            histograms["muon_pt_denom_zone%i" % z].Fill(part.pt)
            if ipart in map_of_iparts_1:
              histograms["muon_pt_fiducial_zone%i" % z].Fill(part.pt)

    # Get hits
    if highest_pt_ipart != -1:
      part = evt.genparticles[highest_pt_ipart]
      highest_pt_part_phi = extrapolate_to_muon(part.phi, np.true_divide(part.q, part.pt))
      highest_pt_part_pt = part.pt
      highest_pt_part_q = part.q
      #print ".. selected part", highest_pt_ipart, highest_pt_part_phi, highest_pt_part_pt, highest_pt_part_q

      highest_pt_ihits = [-1, -1, -1, -1]

      for istation in xrange(4):
        station = istation + 1
        min_abs_dphi = 9999.
        min_abs_dphi_ihit = -1

        for ihit, hit in enumerate(evt.hits):
          if hit.type == kCSC:
            if hit.station == station:
              if hit.sim_tp1 == highest_pt_ipart:
                abs_dphi = abs(delta_phi(np.deg2rad(hit.sim_phi), highest_pt_part_phi))
                if min_abs_dphi > abs_dphi:
                  min_abs_dphi = abs_dphi
                  min_abs_dphi_ihit = ihit

        highest_pt_ihits[istation] = min_abs_dphi_ihit

      highest_pt_hit_phis = [-1, -1, -1, -1]

      for istation in xrange(4):
        ihit = highest_pt_ihits[istation]
        hit_sim_phi = -999999.
        hit_emtf_phi = -1
        hit_zone_phi = -1
        if ihit != -1:
          hit = evt.hits[highest_pt_ihits[istation]]
          hit_sim_phi = hit.sim_phi
          hit_emtf_phi = hit.emtf_phi
          hit_zone_phi = get_zone_hit(hit.emtf_phi)
        highest_pt_hit_phis[istation] = hit_zone_phi
        #print ".. selected hit", ihit, hit_sim_phi, hit_emtf_phi, hit_zone_phi

      # Save the results
      t = (highest_pt_part_pt, highest_pt_part_q, highest_pt_hit_phis, ievt)
      pattern_bank.append(t)
      pass  # end if


# ______________________________________________________________________________
# End job

full_study = True

if full_study:
  # Save histograms
  with root_open("histos_pm.root", "recreate") as f:
    for k, v in histograms.iteritems():
      v.Write()
    for k, v in histogram2Ds.iteritems():
      v.Write()

  ## Save histograms
  #for hname in ["muon_ptmin2_dxy", "muon_ptmin2_dz", "muon_ptmin2_eta", "muon_ptmin2_bx"]:
  #  outname = hname + ".root"
  #  with root_open(outname, "recreate") as f:
  #    histograms[hname].Write()

  # Print patterns
  pattern_bank_sorted = sorted(pattern_bank, key = lambda x: x[0], reverse=True)
  for x in pattern_bank_sorted:
    if x[0] > 1.:
      print x
