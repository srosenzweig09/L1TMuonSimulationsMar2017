import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from rootpy.memory import keepalive


# ______________________________________________________________________________
# Analyzer

# Open file
infile = root_open('ntuple.3.root')
tree = infile.ntupler.tree
maxEvents = -1
#maxEvents = 200000

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM = 0, 1, 2, 3

# Lambdas
deg_to_rad = lambda x: x * np.pi/180.

rad_to_deg = lambda x: x * 180./np.pi


# Book histograms
histograms = {}
histogram2Ds = {}

# GEM-CSC bend
for i in xrange(3):
  # i=0: rear, i=1: front, i=2: all
  hname = "h_bend_ge11_pt5_fr%i" % i
  histograms[hname] = Hist(200, 0, 0.03, name=hname, title="; #Delta#phi(ME1/1,GE1/1) [rad]", type='F')

  hname = "h_bend_ge11_pt20_fr%i" % i
  histograms[hname] = Hist(200, 0, 0.03, name=hname, title="; #Delta#phi(ME1/1,GE1/1) [rad]", type='F')

  hname = "h_bend_ge21_pt5_fr%i" % i
  histograms[hname] = Hist(200, 0, 0.01, name=hname, title="; #Delta#phi(ME2/1,GE2/1) [rad]", type='F')

  hname = "h_bend_ge21_pt20_fr%i" % i
  histograms[hname] = Hist(200, 0, 0.01, name=hname, title="; #Delta#phi(ME2/1,GE2/1) [rad]", type='F')


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
  # Event selection

  no_genparticles = (len(evt.genparticles) == 0)

  #no_genparticles_in_eta_range = (len([part for part in evt.genparticles if 1.7 < abs(part.eta) < 2.1]) == 0)
  no_genparticles_in_eta_range_1 = (len([part for part in evt.genparticles if 1.64 < abs(part.eta) < 2.14]) == 0)

  no_genparticles_pt20 = (len([part for part in evt.genparticles if part.pt > 20]) == 0)

  #no_tracks = (len(evt.tracks) == 0)

  #no_tracks_l1pt20 = (len([trk for trk in evt.tracks if trk.pt > 20.]) == 0)

  #no_tracks_l1qual = (len([trk for trk in evt.tracks if trk.mode in [11,13,14,15]]) == 0)

  #no_tracks_l1qual_l1pt20 = (len([trk for trk in evt.tracks if trk.mode in [11,13,14,15] and trk.pt > 20.]) == 0)


  # Skip events if no gen particles
  if no_genparticles:  continue
  assert len(evt.genparticles) == 1

  # Skip events if no gen particles in eta range
  if no_genparticles_in_eta_range_1:  continue

  mypart = evt.genparticles[0]
  #mytrk = evt.tracks[0]
  #try:
  #  # Select highest pT track from tracks that have a station 1 hit
  #  mytrk = max(filter(lambda x: (x.mode in [11,13,14,15]), evt.tracks), key=lambda x: x.pt)
  #except ValueError:
  #  ## If no tracks have a station 1 hit, just select highest pT track
  #  #mytrk = max(evt.tracks, key=lambda x: x.pt)
  #  raise ValueError(e)


  myhits_me11 = [hit for hit in evt.hits if hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]

  myhits_me21 = [hit for hit in evt.hits if hit.station == 2 and hit.ring == 1 and hit.type == kCSC]

  myhits_ge11 = [hit for hit in evt.hits if hit.station == 1 and hit.ring == 1 and hit.type == kGEM]

  myhits_ge21 = [hit for hit in evt.hits if hit.station == 2 and hit.ring == 1 and hit.type == kGEM]


  # ____________________________________________________________________________
  # Make plots

  try:
    myhit_me11 = np.random.choice(myhits_me11)
    myhit_ge11 = np.random.choice(myhits_ge11)
    mybend_ge11 = myhit_me11.sim_phi - myhit_ge11.sim_phi
    if (mypart.q > 0):  mybend_ge11 = -mybend_ge11
    myfr_ge11 = myhit_me11.fr
    is_good_ge11 = True
  except ValueError, e:
    #raise ValueError(e)
    is_good_ge11 = False

  try:
    myhit_me21 = np.random.choice(myhits_me21)
    myhit_ge21 = np.random.choice(myhits_ge21)
    mybend_ge21 = myhit_me21.sim_phi - myhit_ge21.sim_phi
    if (mypart.q > 0):  mybend_ge21 = -mybend_ge21
    myfr_ge21 = myhit_me21.fr
    is_good_ge21 = True
  except ValueError, e:
    #raise ValueError(e)
    is_good_ge21 = False

  myptbin = -1
  if ((1.0/2 - 0.02) < 1.0/mypart.pt <= (1.0/2)):
    myptbin = 0;  # 2 GeV
  elif ((1.0/3 - 0.02) < 1.0/mypart.pt <= (1.0/3)):
    myptbin = 1;  # 3 GeV
  elif ((1.0/5 - 0.01) < 1.0/mypart.pt <= (1.0/5)):
    myptbin = 2;  # 5 GeV
  elif ((1.0/10 - 0.01) < 1.0/mypart.pt <= (1.0/10)):
    myptbin = 3;  # 10 GeV
  elif ((1.0/20 - 0.01) < 1.0/mypart.pt <= (1.0/20)):
    myptbin = 4;  # 20 GeV
  elif ((1.0/50 - 0.005) < 1.0/mypart.pt <= (1.0/50)):
    myptbin = 5;  # 50 GeV

  # GEM-CSC scaled bend residuals
  if myptbin == 2:
    if is_good_ge11:
      hname = "h_bend_ge11_pt5_fr%i" % myfr_ge11
      histograms[hname].fill(deg_to_rad(mybend_ge11))
      hname = "h_bend_ge11_pt5_fr%i" % 2  # inclusive
      histograms[hname].fill(deg_to_rad(mybend_ge11))
    if is_good_ge21:
      hname = "h_bend_ge21_pt5_fr%i" % myfr_ge21
      histograms[hname].fill(deg_to_rad(mybend_ge21))
      hname = "h_bend_ge21_pt5_fr%i" % 2  # inclusive
      histograms[hname].fill(deg_to_rad(mybend_ge21))

  elif myptbin == 4:
    if is_good_ge11:
      hname = "h_bend_ge11_pt20_fr%i" % myfr_ge11
      histograms[hname].fill(deg_to_rad(mybend_ge11))
      hname = "h_bend_ge11_pt20_fr%i" % 2  # inclusive
      histograms[hname].fill(deg_to_rad(mybend_ge11))
    if is_good_ge21:
      hname = "h_bend_ge21_pt20_fr%i" % myfr_ge21
      histograms[hname].fill(deg_to_rad(mybend_ge21))
      hname = "h_bend_ge21_pt20_fr%i" % 2  # inclusive
      histograms[hname].fill(deg_to_rad(mybend_ge21))


# ______________________________________________________________________________
# Drawer

make_plots = True

if make_plots:
  from drawer import *
  mydrawer = MyDrawer()
  options = mydrawer.options

  # Print
  for hname, h in histograms.iteritems():
    h.Draw("hist")
    print("TH1: {0} (N={1})".format(hname, h.GetEntries()))
    gPad.Print(options.outdir + hname + ".png")
  for hname, h in histogram2Ds.iteritems():
    h.Draw("COLZ")
    print("TH2: {0} (N={1})".format(hname, h.GetEntries()))
    gPad.Print(options.outdir + hname + ".png")


  # Apply Sven's styles
  for i in xrange(2):
    hname = "h_bend_ge11_pt20_fr%i" % i
    h = histograms[hname]
    h.linewidth = 2
    h.linecolor = 'blue'
    h.Scale(1.0/h.Integral())
    h.SetStats(0)
    h.SetYTitle("Arbitrary unit")
    h.Draw("hist")
    #
    hname = "h_bend_ge11_pt5_fr%i" % i
    h = histograms[hname]
    h.linewidth = 2
    h.linecolor = 'red'
    h.Scale(1.0/h.Integral())
    h.Draw("same hist")
    #
    hname = "h_bend_ge11_pt99_fr%i" % i
    gPad.Print(options.outdir + hname + ".png")

  for i in xrange(2):
    hname = "h_bend_ge21_pt20_fr%i" % i
    h = histograms[hname]
    h.linewidth = 2
    h.linecolor = 'blue'
    h.Scale(1.0/h.Integral())
    h.SetStats(0)
    h.SetYTitle("Arbitrary unit")
    h.Draw("hist")
    #
    hname = "h_bend_ge21_pt5_fr%i" % i
    h = histograms[hname]
    h.linewidth = 2
    h.linecolor = 'red'
    h.Scale(1.0/h.Integral())
    h.Draw("same hist")
    #
    hname = "h_bend_ge21_pt99_fr%i" % i
    gPad.Print(options.outdir + hname + ".png")
