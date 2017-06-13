import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Efficiency, Legend, Canvas
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
#infile = root_open('ntuple.0.root')
infile = root_open('ntuple.1.root')
tree = infile.ntupler.tree
maxEvents = -1
#maxEvents = 20000

# Book histograms
histograms = {}
histogram2Ds = {}

# pT vs gen pT
hname = "h2_pt_vs_genpt"
histogram2Ds[hname] = Hist2D(50, 0.0, 0.5, 50, 0.0, 0.5, name=hname, title="; gen 1/p_{T} [1/GeV]; EMTFv5 1/p_{T} [1/GeV]", type='F')

# GEM-CSC bend vs gen pT
for i in xrange(3):
  # i=0: rear, i=1: front, i=2: all
  hname = "h2_bend_vs_genpt_fr%i" % i
  histogram2Ds[hname] = Hist2D(50, 0.0, 0.5, 71, -1.025, 2.525, name=hname, title="; gen 1/p_{T} [1/GeV]; GEM-CSC bend [deg]", type='F')

  hname = "h2_common_bend_vs_genpt_fr%i" % i
  histogram2Ds[hname] = Hist2D(50, 0.0, 0.5, 71, -1.025, 2.525, name=hname, title="; gen 1/p_{T} [1/GeV]; GEM-CSC bend [deg]", type='F')

# GEM-CSC bend when triggered
for i in xrange(3):
  # i=0: gen pT <= 20, i=1: gen pT > 20, i=2: all
  hname = "h_bend_l1pt20un_gen%i" % i
  histograms[hname] = Hist(71, -1.025, 2.525, name=hname, title="; GEM-CSC bend [deg]", type='F')

  # i=0: gen pT <= 3, i=1: gen pT > 3, i=2: all
  hname = "h_bend_l1pt3un_gen%i" % i
  histograms[hname] = Hist(71, -1.025, 2.525, name=hname, title="; GEM-CSC bend [deg]", type='F')

# GEM-CSC scaled bend residuals
for i in xrange(6):
  # i=0..6: 2, 3, 5, 10, 20, 50 GeV
  hname = "h_bend_s_pt%i" % i
  histograms[hname] = Hist(101, -0.505, 0.505, name=hname, title="; residual of (scaled bend - EMTFv5 1/p_{T}) [1/GeV]", type='F')

# Efficiency vs gen pT
for k in ["numer", "denom"]:
  hname = "eff_vs_genpt_l1pt20_%s" % k
  histograms[hname] = Hist(50, 0.0, 100.0, name=hname, title="; gen p_{T} [GeV]", type='F')
  hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
  histograms[hname] = Hist(50, 0.0, 100.0, name=hname, title="; gen p_{T} [GeV]", type='F')
  hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
  histograms[hname] = Hist(50, 0.0, 100.0, name=hname, title="; gen p_{T} [GeV]", type='F')


# Define collection
tree.define_collection(name='hits', prefix='vh_', size='vh_size')
tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
tree.define_collection(name='genparticles', prefix='vp_', size='vp_size')

# Enums
kDT, kCSC, kRPC, kGEM = 0, 1, 2, 3

# Lambdas
unscale_pt = lambda x: x/1.4

get_common_bend = lambda x, y: (x*2.29072) if y else x

get_scaled_bend = lambda x: (x*0.262506)


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
  # Quick efficiency studies

  no_genparticles = (len(evt.genparticles) == 0)

  no_hits_me11 = (len([hit for hit in evt.hits if hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]) == 0)

  no_hits_ge11 = (len([hit for hit in evt.hits if hit.station == 1 and hit.ring == 1 and hit.type == kGEM]) == 0)

  no_tracks = (len(evt.tracks) == 0)

  no_tracks_l1pt20 = (len([trk for trk in evt.tracks if trk.pt > 20.]) == 0)

  if not no_genparticles and not no_hits_me11 and not no_hits_ge11:
    #trigger = not no_tracks
    trigger = not no_tracks_l1pt20
    k = "denom"
    hname = "eff_vs_genpt_l1pt20_%s" % k
    histograms[hname].fill(evt.genparticles[0].pt)
    if trigger:
      k = "numer"
      hname = "eff_vs_genpt_l1pt20_%s" % k
      histograms[hname].fill(evt.genparticles[0].pt)

    k = "denom"
    hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
    histograms[hname].fill(evt.genparticles[0].pt)
    #if trigger:
    #  k = "numer"
    #  hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
    #  histograms[hname].fill(evt.genparticles[0].pt)

    k = "denom"
    hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
    histograms[hname].fill(evt.genparticles[0].pt)
    #if trigger:
    #  k = "numer"
    #  hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
    #  histograms[hname].fill(evt.genparticles[0].pt)

  # ____________________________________________________________________________
  # Filter

  # Skip events if no gen particles
  if no_genparticles:  continue
  assert len(evt.genparticles) == 1

  # Skip events if no tracks
  if no_tracks:  continue

  mypart = evt.genparticles[0]
  #mytrk = evt.tracks[0]
  mytrk = max(evt.tracks, key=lambda x: x.pt)

  # Skip event if no station 1 hit
  mode = mytrk.mode
  if not (mode & (1<<3)):  continue

  # Skip event if no ME1/1 hit
  myhits_me11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.type == kCSC]
  if not len(myhits_me11): continue

  # Skip event if no ME1/1 hit
  myhits_ge11 = [hit for hit in evt.hits if hit.endcap == mytrk.endcap and hit.sector == mytrk.sector and hit.station == 1 and hit.ring == 1 and hit.type == kGEM]
  if not len(myhits_ge11): continue

  # ____________________________________________________________________________
  # Make plots

  myhit_me11 = np.random.choice(myhits_me11)
  myhit_ge11 = np.random.choice(myhits_ge11)
  mybend = myhit_me11.sim_phi - myhit_ge11.sim_phi
  if (mypart.q > 0):  mybend = -mybend
  myfr = myhit_me11.fr
  mybend_common = get_common_bend(mybend, myfr)
  mybend_scaled = get_scaled_bend(mybend_common)

  myptbin = -1
  if ((1.0/2 - 0.05) < 1.0/mypart.pt <= (1.0/2)):
    myptbin = 0;
  elif ((1.0/3 - 0.03333) < 1.0/mypart.pt <= (1.0/3)):
    myptbin = 1;
  elif ((1.0/5 - 0.03333) < 1.0/mypart.pt <= (1.0/5)):
    myptbin = 2;
  elif ((1.0/10 - 0.02) < 1.0/mypart.pt <= (1.0/10)):
    myptbin = 3;
  elif ((1.0/20 - 0.02) < 1.0/mypart.pt <= (1.0/20)):
    myptbin = 4;
  elif ((1.0/50 - 0.01) < 1.0/mypart.pt <= (1.0/50)):
    myptbin = 5;

  # pT vs gen pT
  hname = "h2_pt_vs_genpt"
  histogram2Ds[hname].fill(1.0/mypart.pt, 1.0/unscale_pt(mytrk.pt))

  # bend vs gen pT
  hname = "h2_bend_vs_genpt_fr%i" % myfr
  histogram2Ds[hname].fill(1.0/mypart.pt, mybend)
  hname = "h2_bend_vs_genpt_fr%i" % 2  # inclusive
  histogram2Ds[hname].fill(1.0/mypart.pt, mybend)

  # common bend vs gen pT
  hname = "h2_common_bend_vs_genpt_fr%i" % myfr
  histogram2Ds[hname].fill(1.0/mypart.pt, mybend_common)
  hname = "h2_common_bend_vs_genpt_fr%i" % 2  # inclusive
  histogram2Ds[hname].fill(1.0/mypart.pt, mybend_common)

  # GEM-CSC bend when triggered
  pass_l1pt = unscale_pt(mytrk.pt) > 20.
  pass_genpt = mypart.pt > 20.
  if pass_l1pt:
    hname = "h_bend_l1pt20un_gen%i" % pass_genpt
    histograms[hname].fill(mybend_common)
    hname = "h_bend_l1pt20un_gen%i" % 2  # inclusive
    histograms[hname].fill(mybend_common)

  pass_l1pt = unscale_pt(mytrk.pt) > 3.
  pass_genpt = mypart.pt > 3.
  if pass_l1pt:
    hname = "h_bend_l1pt3un_gen%i" % pass_genpt
    histograms[hname].fill(mybend_common)
    hname = "h_bend_l1pt3un_gen%i" % 2  # inclusive
    histograms[hname].fill(mybend_common)

  # GEM-CSC scaled bend residuals
  mybend_scaled_residual = mybend_scaled - 1.0/unscale_pt(mytrk.pt)
  if myptbin != -1:
    hname = "h_bend_s_pt%i" % myptbin
    histograms[hname].fill(mybend_scaled_residual)

  # ____________________________________________________________________________
  # Quick efficiency studies (con't)

  pass_l1bend1 = abs(mybend_scaled_residual) <= 0.05
  pass_l1bend2 = abs(mybend_scaled_residual) <= 0.03

  if not no_genparticles and not no_hits_me11 and not no_hits_ge11:
    #trigger = not no_tracks
    trigger = not no_tracks_l1pt20
    trigger = trigger and pass_l1bend1
    if trigger:
      k = "numer"
      hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
      histograms[hname].fill(evt.genparticles[0].pt)

    trigger = not no_tracks_l1pt20
    trigger = trigger and pass_l1bend2
    if trigger:
      k = "numer"
      hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
      histograms[hname].fill(evt.genparticles[0].pt)

  continue  # end loop over event

# ______________________________________________________________________________
# Drawer

make_plots = True

if make_plots:
  from drawer import *
  mydrawer = MyDrawer()
  options = mydrawer.options

  # Print
  for hname, h in histograms.iteritems():
    h.Draw("COLZ")
    gPad.Print(options.outdir + hname + ".png")
  for hname, h in histogram2Ds.iteritems():
    h.Draw("COLZ")
    gPad.Print(options.outdir + hname + ".png")

  # Make ratio of bend vs gen pT
  denom_ifr = 1
  hname = "h2_bend_vs_genpt_fr%i" % denom_ifr
  h2 = histogram2Ds[hname].Clone(hname + "_tmp")
  h2.RebinX(4)
  prof = h2.ProfileX(hname + "_tmp_pfx", 1, -1, "s")
  proj = prof.ProjectionX(hname + "_tmp_px", "e")
  denom = proj
  #
  numer_ifr = 1 - denom_ifr
  hname = "h2_bend_vs_genpt_fr%i" % numer_ifr
  h2 = histogram2Ds[hname].Clone(hname + "_tmp")
  h2.RebinX(4)
  prof = h2.ProfileX(hname + "_tmp_pfx", 1, -1, "s")
  proj = prof.ProjectionX(hname + "_tmp_px", "e")
  numer = proj
  #
  numer.Divide(denom)
  numer.Draw()
  numer.Fit("pol0", "", "", 0.04, 0.32)
  hname = numer.GetName()
  gPad.Print(options.outdir + hname + ".png")

  # Make ratio of bend vs gen pT [2]
  hname = "h2_common_bend_vs_genpt_fr%i" % 2  # inclusive
  h2 = histogram2Ds[hname].Clone(hname + "_tmp")
  h2.RebinX(4)
  prof = h2.ProfileX(hname + "_tmp_pfx", 1, -1, "s")
  proj = prof.ProjectionX(hname + "_tmp_px", "e")
  denom = proj
  #
  hname = "h2_pt_vs_genpt"
  h2 = histogram2Ds[hname].Clone(hname + "_tmp")
  h2.RebinX(4)
  prof = h2.ProfileX(hname + "_tmp_pfx", 1, -1, "s")
  proj = prof.ProjectionX(hname + "_tmp_px", "e")
  numer = proj
  #
  numer.Divide(denom)
  numer.Draw()
  numer.Fit("pol0", "", "", 0.04, 0.32)
  hname = numer.GetName()
  gPad.Print(options.outdir + hname + ".png")

  # Make overlay of GEM-CSC bend when triggered
  hname = "h_bend_l1pt20un_gen%i" % 2  # inclusive
  h1a = histograms[hname]
  h1a.linecolor = 'black'
  h1a.linewidth = 2
  hname = "h_bend_l1pt20un_gen%i" % 1
  h1b = histograms[hname]
  h1b.linecolor = 'red'
  h1b.linewidth = 2
  h1a.Draw("hist")
  h1b.Draw("same hist")
  hname = "h_bend_l1pt20un_gen%i" % 99
  gPad.Print(options.outdir + hname + ".png")

  hname = "h_bend_l1pt3un_gen%i" % 2  # inclusive
  h1a = histograms[hname]
  h1a.linecolor = 'black'
  h1a.linewidth = 2
  hname = "h_bend_l1pt3un_gen%i" % 1
  h1b = histograms[hname]
  h1b.linecolor = 'red'
  h1b.linewidth = 2
  h1a.Draw("hist")
  h1b.Draw("same hist")
  hname = "h_bend_l1pt3un_gen%i" % 99
  gPad.Print(options.outdir + hname + ".png")

  # Make overlay of GEM-CSC scaled bend residuals
  labels = ["2 GeV", "3 GeV", "5 GeV", "10 GeV", "20 GeV", "50 GeV"]
  leg = Legend(len(labels), leftmargin=0.55, textfont=42, textsize=0.03, entryheight=0.04, entrysep=0.01)
  leg.SetShadowColor(0)
  leg.SetBorderSize(0)
  #
  for i in xrange(6):
    hname = "h_bend_s_pt%i" % i
    h1a = histograms[hname]
    h1a.linecolor = options.palette[i]
    h1a.linewidth = 2
    h1a.Scale(1.0/h1a.Integral())  # normalize
    if i == 0:
      ymax = 0.3
      h1a.SetMaximum(ymax)
      h1a.GetYaxis().SetTitle("(normalized)")
      h1a.Draw("hist")
    else:
      h1a.Draw("same hist")
    leg.AddEntry(h1a, labels[i], "l")
  #
  leg.Draw()
  hname = "h_bend_s_pt%i" % 99
  gPad.Print(options.outdir + hname + ".png")

  # Study GEM-CSC scaled bend residual for pT = 20
  i = 4
  hname = "h_bend_s_pt%i" % i
  h1a = histograms[hname]
  in_quantiles = np.array([0.01/2, 0.02/2, 0.05/2, 0.10/2, 0.15/2, 0.20/2, 0.25/2, 1.0-0.25/2, 1.0-0.20/2, 1.0-0.15/2, 1.0-0.10/2, 1.0-0.05/2, 1.0-0.02/2, 1.0-0.01/2], dtype=np.float64)
  quantiles = np.array([0.] * len(in_quantiles), dtype=np.float64)
  h1a.GetQuantiles(len(in_quantiles), quantiles, in_quantiles)
  print "Quantiles for %s" % hname
  for in_q, q in zip(in_quantiles, quantiles):
    print "..", in_q, q

  # Make efficiency vs gen pT
  k = "denom"
  hname = "eff_vs_genpt_l1pt20_%s" % k
  denom = histograms[hname]
  k = "numer"
  hname = "eff_vs_genpt_l1pt20_%s" % k
  numer = histograms[hname]
  hname = "eff_vs_genpt_l1pt20"
  eff = Efficiency(numer, denom, name=hname)
  eff.SetStatisticOption(0)  # kFCP
  eff.SetConfidenceLevel(0.682689492137)  # one sigma
  eff.linecolor = 'black'
  eff.linewidth = 2
  eff.markercolor = 'black'
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
  gPad.Print(options.outdir + hname + ".png")
  #
  k = "denom"
  hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
  denom = histograms[hname]
  k = "numer"
  hname = "eff_vs_genpt_l1pt20_l1bend1_%s" % k
  numer = histograms[hname]
  hname = "eff_vs_genpt_l1pt20_l1bend1"
  eff2 = Efficiency(numer, denom, name=hname)
  eff2.SetStatisticOption(0)  # kFCP
  eff2.SetConfidenceLevel(0.682689492137)  # one sigma
  eff2.linecolor = 'red'
  eff2.linewidth = 2
  eff2.markercolor = 'red'
  eff2.markerstyle = 1
  eff2.Draw("same p")
  gPad.Print(options.outdir + hname + ".png")
  #
  k = "denom"
  hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
  denom = histograms[hname]
  k = "numer"
  hname = "eff_vs_genpt_l1pt20_l1bend2_%s" % k
  numer = histograms[hname]
  hname = "eff_vs_genpt_l1pt20_l1bend2"
  eff2 = Efficiency(numer, denom, name=hname)
  eff2.SetStatisticOption(0)  # kFCP
  eff2.SetConfidenceLevel(0.682689492137)  # one sigma
  eff2.linecolor = 'blue'
  eff2.linewidth = 2
  eff2.markercolor = 'blue'
  eff2.markerstyle = 1
  eff2.Draw("same p")
  gPad.Print(options.outdir + hname + ".png")
