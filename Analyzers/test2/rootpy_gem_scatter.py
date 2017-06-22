import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from rootpy.memory import keepalive


# ______________________________________________________________________________
# Analyzer

# Open file
infile = root_open('ntuple.3.root')
tree = infile.ntupler.tree
#maxEvents = -1
maxEvents = 40000

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM = 0, 1, 2, 3

# Lambdas
unscale_pt = lambda x: x/1.4

get_common_bend = lambda x, y: (x*2.3873) if y else x

get_scaled_bend = lambda x: (x*0.27302)

# Functions
def delta_phi(lhs, rhs):  # in degrees
  deg = lhs - rhs
  while deg <  -180.:  deg += 360.
  while deg >= +180.:  deg -= 360.
  return deg


# Book histograms
histograms = {}
histogram2Ds = {}

# GEM-CSC bend
for i in xrange(3):
  # i=0: bkg, i=1: signal, i=2: all
  hname = "h2_bend_signal%i" % i
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; EMTFv5 p_{T} (x1.4) [GeV]; GEM-CSC bend [deg]", type='F')
  histograms[hname].data = []
  histograms[hname].ymin = -0.5
  histograms[hname].ymax = 3.0

  hname = "h2_common_bend_signal%i" % i
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; EMTFv5 p_{T} (x1.4) [GeV]; GEM-CSC bend [deg]", type='F')
  histograms[hname].data = []
  histograms[hname].ymin = -0.5
  histograms[hname].ymax = 3.0

  hname = "h2_scaled_bend_signal%i" % i
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; EMTFv5 p_{T} (x1.4) [GeV]; residual of (GEM 1/p_{T} - EMTFv5 1/p_{T}) [1/GeV]", type='F')
  histograms[hname].data = []
  histograms[hname].ymin = -0.5
  histograms[hname].ymax = 0.5

  hname = "h2_abs_scaled_bend_signal%i" % i
  histograms[hname] = Hist(200, 0, 50, name=hname, title="; EMTFv5 p_{T} (x1.4) [GeV]; |residual| of (GEM 1/p_{T} - EMTFv5 1/p_{T}) [1/GeV]", type='F')
  histograms[hname].data = []
  histograms[hname].ymin = 0.0
  histograms[hname].ymax = 0.5


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

  no_genparticles_in_eta_range = (len([part for part in evt.genparticles if 1.7 < abs(part.eta) < 2.1]) == 0)

  no_genparticles_pt20 = (len([part for part in evt.genparticles if part.pt > 20]) == 0)

  no_tracks = (len(evt.tracks) == 0)

  no_tracks_l1pt20 = (len([trk for trk in evt.tracks if trk.pt > 20.]) == 0)

  no_tracks_l1qual = (len([trk for trk in evt.tracks if trk.mode in [11,13,14,15]]) == 0)

  no_tracks_l1qual_l1pt20 = (len([trk for trk in evt.tracks if trk.mode in [11,13,14,15] and trk.pt > 20.]) == 0)

  # ____________________________________________________________________________
  # Filter

  # Skip events if no gen particles
  if no_genparticles:  continue
  assert len(evt.genparticles) == 1

  # Skip events if no tracks
  if no_tracks:  continue

  # Skip events if no tracks with SingleMu quality
  if no_tracks_l1qual:  continue

  mypart = evt.genparticles[0]
  #mytrk = evt.tracks[0]
  try:
    # Select highest pT track from tracks that have a station 1 hit
    mytrk = max(filter(lambda x: (x.mode in [11,13,14,15]), evt.tracks), key=lambda x: x.pt)
  except ValueError, e:
    ## If no tracks have a station 1 hit, just select highest pT track
    #mytrk = max(evt.tracks, key=lambda x: x.pt)
    raise ValueError(e)

  # Skip event if no ME1/1 hit
  myhits_me11 = [hit for hit in evt.hits if hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]
  #myhits_me11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]
  if len(myhits_me11) == 0:  continue

  # Skip event if no GE1/1 hit
  myhits_ge11 = [hit for hit in evt.hits if hit.station == 1 and hit.ring == 1 and hit.type == kGEM]
  #myhits_ge11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and hit.ring == 1 and hit.type == kGEM]
  if len(myhits_ge11) == 0:  continue

  # ____________________________________________________________________________
  # Make plots

  myhit_me11 = np.random.choice(myhits_me11)
  myhit_ge11 = np.random.choice(myhits_ge11)
  mybend = delta_phi(myhit_me11.sim_phi, myhit_ge11.sim_phi)
  if (mypart.q > 0):  mybend = -mybend
  myfr = myhit_me11.fr
  mybend_common = get_common_bend(mybend, myfr)
  mybend_scaled = get_scaled_bend(mybend_common)
  mybend_scaled_residual = mybend_scaled - 1.0/unscale_pt(mytrk.pt)

  if mypart.pt > 20:
    signal = 1
  else:
    signal = 0

  hname = "h2_bend_signal%i" % signal
  histograms[hname].data.append((mytrk.pt, mybend))
  hname = "h2_bend_signal%i" % 2  # inclusive
  histograms[hname].data.append((mytrk.pt, mybend))

  hname = "h2_common_bend_signal%i" % signal
  histograms[hname].data.append((mytrk.pt, mybend_common))
  hname = "h2_common_bend_signal%i" % 2  # inclusive
  histograms[hname].data.append((mytrk.pt, mybend_common))

  hname = "h2_scaled_bend_signal%i" % signal
  histograms[hname].data.append((mytrk.pt, mybend_scaled_residual))
  hname = "h2_scaled_bend_signal%i" % 2  # inclusive
  histograms[hname].data.append((mytrk.pt, mybend_scaled_residual))

  hname = "h2_abs_scaled_bend_signal%i" % signal
  histograms[hname].data.append((mytrk.pt, abs(mybend_scaled_residual)))
  hname = "h2_abs_scaled_bend_signal%i" % 2  # inclusive
  histograms[hname].data.append((mytrk.pt, abs(mybend_scaled_residual)))

  continue  # end loop over event


# ______________________________________________________________________________
# Drawer

make_plots = True

def make_scatter_plots_pls(sig_name, bkg_name):
  hname = sig_name
  h = histograms[hname]
  h.SetMinimum(h.ymin)
  h.SetMaximum(h.ymax)
  h.Draw()
  n = len(h.data)
  #
  hname = bkg_name
  h = histograms[hname]
  g = Graph(n)      # truncate
  data = h.data[:n] # truncate
  for i, (xx, yy) in enumerate(data):
    g.SetPoint(i, xx, yy)
  g.markerstyle = 4
  g.markercolor = 'black'
  g.Draw("p")
  keepalive(gPad.func(), g)
  #
  hname = sig_name
  h = histograms[hname]
  g = Graph(n)
  data = h.data
  for i, (xx, yy) in enumerate(data):
    g.SetPoint(i, xx, yy)
  g.markerstyle = 4
  g.markercolor = 'red'
  g.Draw("p")
  keepalive(gPad.func(), g)
  gPad.Print(options.outdir + hname + ".png")
  return

if make_plots:
  from drawer import *
  mydrawer = MyDrawer()
  options = mydrawer.options

  # Make scatter plots
  make_scatter_plots_pls("h2_bend_signal1", "h2_bend_signal0")
  make_scatter_plots_pls("h2_common_bend_signal1", "h2_common_bend_signal0")
  make_scatter_plots_pls("h2_scaled_bend_signal1", "h2_scaled_bend_signal0")
  make_scatter_plots_pls("h2_abs_scaled_bend_signal1", "h2_abs_scaled_bend_signal0")
