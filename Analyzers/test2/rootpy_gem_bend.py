import numpy as np

from rootpy.plotting import Hist, Hist2D
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
# Open file
infile = root_open('ntuple.0.root')
tree = infile.ntupler.tree

# Create histograms
h2_pt_vs_genpt = Hist2D(40, 0.0, 0.5, 40, 0.0, 0.5, name='h2_pt_vs_genpt', title="; gen 1/p_{T} [1/GeV]; EMTFv5 1/p_{T} [1/GeV]", type='F')
h2_bend_vs_genpt = Hist2D(40, 0.0, 0.5, 40, 0.0, 2.5, name='h2_bend_vs_genpt', title="; gen 1/p_{T} [1/GeV]; GEM-CSC bend [deg]", type='F')

# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM = 0, 1, 2, 3

# Loop over events
for ievt, evt in enumerate(tree):
  if (ievt % 1 == 0):  print("Processing event: {0}".format(ievt))

  # ____________________________________________________________________________
  # Filter

  # Skip events if not exactly one gen particle
  if len(evt.genparticles) != 1:  continue

  # Skip events if not exactly one track
  if len(evt.tracks) != 1:  continue

  mytrk = evt.tracks[0]
  mypart = evt.genparticles[0]

  # Skip event if no station 1 hit
  mode = mytrk.mode
  if not (mode & (1<<3)):  continue

  endcap = mytrk.endcap
  sector = mytrk.sector

  # Skip event if no ME1/1 hit
  myhits_me11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]
  if not len(myhits_me11): continue

  # Skip event if no ME1/1 hit
  myhits_ge11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and hit.ring == 1 and hit.type == kGEM]
  if not len(myhits_ge11): continue


  # ____________________________________________________________________________
  # Verbose

  # Hits
  for ihit, hit in enumerate(evt.hits):
    print(".. hit  {0} {1} {2} {3} {4} {5} {6}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sim_phi, hit.sim_theta))

  # Tracks
  for itrk, trk in enumerate(evt.tracks):
    print(".. trk  {0} {1} {2} {3} {4} {5} {6}".format(itrk, trk.pt, trk.phi, trk.eta, trk.theta, trk.q, trk.mode))

  # Gen particles
  for ipart, part in enumerate(evt.genparticles):
    print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))

  # ____________________________________________________________________________
  # Make plots

  myhit_me11 = np.random.choice(myhits_me11)
  myhit_ge11 = np.random.choice(myhits_ge11)
  mybend = abs(myhit_me11.sim_phi - myhit_ge11.sim_phi)

  h2_pt_vs_genpt.Fill(1.0/mypart.pt, 1.0/mytrk.pt)  # x=mypart, y=mytrk
  h2_bend_vs_genpt.Fill(1.0/mypart.pt, mybend)

  continue  # end loop over event

# ______________________________________________________________________________
# Draw

from ROOT import gROOT, gStyle, gPad

gROOT.SetBatch(True)
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle();")

#gStyle.SetOptStat(111110)
gStyle.SetPadRightMargin(0.05)
gStyle.SetTitleOffset(1.1, "Y")
gStyle.SetPalette(57)  # kBird
gStyle.SetNumberContours(50)

h2_pt_vs_genpt.Draw("COLZ")
gPad.Print("figures/h2_pt_vs_genpt.png")

h2_bend_vs_genpt.Draw("COLZ")
gPad.Print("figures/h2_bend_vs_genpt.png")
