import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open
from rootpy.memory import keepalive


# ______________________________________________________________________________
# Analyzer

# Open file
infile_copad1 = root_open('ntuple.2.root')
infile_copad2 = root_open('ntuple.3.root')
tree_copad1 = infile_copad1.ntupler.tree
tree_copad2 = infile_copad2.ntupler.tree

maxEvents = -1
#maxEvents = 200000

# Define collection
for tree in [tree_copad1, tree_copad2]:
   tree.define_collection(name='hits', prefix='vh_', size='vh_size')
   tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
   tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM = 0, 1, 2, 3


# Book histograms
histograms = {}
histogram2Ds = {}

# Efficiency vs gen |eta|
for k in ["denom", "numer1", "numer2"]:
  hname = "eff_vs_geneta_ge11_%s" % k
  histograms[hname] = Hist(100, 1.4, 2.4, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_geneta_ge21_%s" % k
  histograms[hname] = Hist(100, 1.4, 2.4, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_geneta_me11_ge11_%s" % k
  histograms[hname] = Hist(100, 1.4, 2.4, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()

  hname = "eff_vs_geneta_me21_ge21_%s" % k
  histograms[hname] = Hist(100, 1.4, 2.4, name=hname, title="; gen |#eta|", type='F')
  histograms[hname].Sumw2()


# Loop over trees
treeinfos = [
  (1, tree_copad1, "numer1"),
  (2, tree_copad2, "numer2"),
]

for (itree, tree, tree_label) in treeinfos:
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
    # Hit reconstruction efficiency

    no_genparticles = (len(evt.genparticles) == 0)

    #no_genparticles_in_eta_range = (len([part for part in evt.genparticles if 1.7 < abs(part.eta) < 2.1]) == 0)

    no_genparticles_pt20 = (len([part for part in evt.genparticles if part.pt > 20]) == 0)

    no_hits_me11 = (len([hit for hit in evt.hits if hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]) == 0)

    no_hits_me21 = (len([hit for hit in evt.hits if hit.station == 2 and hit.ring == 1 and hit.type == kCSC]) == 0)

    no_hits_ge11 = (len([hit for hit in evt.hits if hit.station == 1 and hit.ring == 1 and hit.type == kGEM]) == 0)

    no_hits_ge21 = (len([hit for hit in evt.hits if hit.station == 2 and hit.ring == 1 and hit.type == kGEM]) == 0)

    if not no_genparticles and not no_genparticles_pt20:
      # Denominators
      if tree_label == "numer1":
        k = "denom"
        hname = "eff_vs_geneta_ge11_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

        k = "denom"
        hname = "eff_vs_geneta_ge21_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

        k = "denom"
        hname = "eff_vs_geneta_me11_ge11_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

        k = "denom"
        hname = "eff_vs_geneta_me21_ge21_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

      # Numerators
      trigger_ge11 = not no_hits_ge11
      trigger_ge21 = not no_hits_ge21
      trigger_me11_ge11 = not no_hits_me11 and not no_hits_ge11
      trigger_me21_ge21 = not no_hits_me21 and not no_hits_ge21

      if trigger_ge11:
        k = tree_label
        hname = "eff_vs_geneta_ge11_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

      if trigger_ge21:
        k = tree_label
        hname = "eff_vs_geneta_ge21_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

      if trigger_me11_ge11:
        k = tree_label
        hname = "eff_vs_geneta_me11_ge11_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

      if trigger_me21_ge21:
        k = tree_label
        hname = "eff_vs_geneta_me21_ge21_%s" % k
        histograms[hname].fill(abs(evt.genparticles[0].eta))

  continue  # end loop over event


# ______________________________________________________________________________
# Drawer

make_plots = True

def make_efficiency_plots_pls(names, outname, colors):
  if len(names) < 2:
    raise ValueError("I need 2 or more histograms!")

  if len(names) != len(colors):
    raise ValueError("I need a color for each histogram!")

  # First
  n = len(names)
  i = 1
  hname = names[0]
  denom = histograms[hname]
  hname = names[i]
  numer = histograms[hname]
  hname = outname
  eff = Efficiency(numer, denom, name=hname)
  eff.SetStatisticOption(0)  # kFCP
  eff.SetConfidenceLevel(0.682689492137)  # one sigma
  eff.linecolor = colors[i]
  eff.linewidth = 2
  eff.markercolor = colors[i]
  eff.markerstyle = 1
  #
  frame = eff.GetCopyTotalHisto().Clone(hname+"_frame")
  frame.Reset()
  frame.SetMinimum(0)
  frame.SetMaximum(1.2)
  frame.GetYaxis().SetTitle("#varepsilon")
  frame.SetStats(0)
  frame.Draw()
  tline = TLine()
  tline.SetLineColor(1)
  xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
  tline.DrawLine(xmin, 1.0, xmax, 1.0)
  #
  eff.Draw("same p")
  keepalive(gPad.func(), eff)
  gPad.Print(options.outdir + hname + ".png")

  # Second and on
  for i in xrange(2, n):
    hname = names[0]
    denom = histograms[hname]
    hname = names[i]
    numer = histograms[hname]
    hname = outname
    eff = Efficiency(numer, denom, name=hname)
    eff.SetStatisticOption(0)  # kFCP
    eff.SetConfidenceLevel(0.682689492137)  # one sigma
    eff.linecolor = colors[i]
    eff.linewidth = 2
    eff.markercolor = colors[i]
    eff.markerstyle = 1
    eff.Draw("same p")
    keepalive(gPad.func(), eff)
    gPad.Print(options.outdir + hname + ".png")
  return

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


  # Make efficiency vs gen |eta|
  make_efficiency_plots_pls(["eff_vs_geneta_ge11_denom", "eff_vs_geneta_ge11_numer1", "eff_vs_geneta_ge11_numer2"], "eff_vs_geneta_ge11", ["whatever", "black", "blue"])
  make_efficiency_plots_pls(["eff_vs_geneta_ge21_denom", "eff_vs_geneta_ge21_numer1", "eff_vs_geneta_ge21_numer2"], "eff_vs_geneta_ge21", ["whatever", "black", "blue"])
  make_efficiency_plots_pls(["eff_vs_geneta_me11_ge11_denom", "eff_vs_geneta_me11_ge11_numer1", "eff_vs_geneta_me11_ge11_numer2"], "eff_vs_geneta_me11_ge11", ["whatever", "black", "blue"])
  make_efficiency_plots_pls(["eff_vs_geneta_me21_ge21_denom", "eff_vs_geneta_me21_ge21_numer1", "eff_vs_geneta_me21_ge21_numer2"], "eff_vs_geneta_me21_ge21", ["whatever", "black", "blue"])

