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
infile = root_open('ntuple_SingleMuon_Toy_add.1.root')
#infile = root_open('ntuple_SingleMuon_Toy_0T_add.1.root')
#infile = root_open('ntuple_SingleMuon_Toy_0M_add.1.root')
tree = infile.ntupler.tree
maxEvents = -1
#maxEvents = 20
print "[INFO] Opening file: %s" % infile

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
def delta_phi(lhs, rhs):  # in radians
  deg = lhs - rhs
  while deg <  -np.pi:  deg += np.pi*2
  while deg >= +np.pi:  deg -= np.pi*2
  return deg

def get_pt_bin(pt):
  ipt = -1
  if ((1.0/20 - 0.005) < 1.0/pt <= (1.0/20)):
    ipt = 4
  return ipt

def get_eta_bin(eta):
  ieta = -1
  if (1.64 < abs(eta) <= 1.74):
    ieta = 0
  elif (1.74 < abs(eta) <= 1.84):
    ieta = 1
  elif (1.84 < abs(eta) <= 1.94):
    ieta = 2
  elif (1.94 < abs(eta) <= 2.04):
    ieta = 3
  elif (2.04 < abs(eta) <= 2.14):
    ieta = 4
  return ieta


# Book histograms
histograms = {}
histogram2Ds = {}

hname = "h2_fr_vs_pt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 50, 2, 0, 2, name=hname, title="; gen p_{T} [GeV]; ME1 FR bit", type='F')
hname = "h2_fr_vs_pt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 50, 2, 0, 2, name=hname, title="; gen p_{T} [GeV]; ME2 FR bit", type='F')
hname = "h2_fr_vs_pt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 50, 2, 0, 2, name=hname, title="; gen p_{T} [GeV]; ME3 FR bit", type='F')
hname = "h2_fr_vs_pt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 50, 2, 0, 2, name=hname, title="; gen p_{T} [GeV]; ME4 FR bit", type='F')

hname = "h2_fr_vs_invpt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 2, 0, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME1 FR bit", type='F')
hname = "h2_fr_vs_invpt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 2, 0, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME2 FR bit", type='F')
hname = "h2_fr_vs_invpt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 2, 0, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME3 FR bit", type='F')
hname = "h2_fr_vs_invpt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 2, 0, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME4 FR bit", type='F')

hname = "h2_dphi_vs_pt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 50, 400, -0.5, 1.5, name=hname, title="; gen p_{T} [GeV]; ME1 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_pt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 50, 400, -0.5, 1.5, name=hname, title="; gen p_{T} [GeV]; ME2 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_pt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 50, 400, -0.5, 1.5, name=hname, title="; gen p_{T} [GeV]; ME3 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_pt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 50, 400, -0.5, 1.5, name=hname, title="; gen p_{T} [GeV]; ME4 #Delta#phi [rad]", type='F')

hname = "h2_dphi_vs_invpt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 400, -0.5, 1.5, name=hname, title="; gen 1/p_{T} [1/GeV]; ME1 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_invpt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 400, -0.5, 1.5, name=hname, title="; gen 1/p_{T} [1/GeV]; ME2 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_invpt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 400, -0.5, 1.5, name=hname, title="; gen 1/p_{T} [1/GeV]; ME3 #Delta#phi [rad]", type='F')
hname = "h2_dphi_vs_invpt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 400, -0.5, 1.5, name=hname, title="; gen 1/p_{T} [1/GeV]; ME4 #Delta#phi [rad]", type='F')

hname = "h2_k_vs_pt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 50, 150, 0, 3, name=hname, title="; gen p_{T} [GeV]; ME1 k-factor", type='F')
hname = "h2_k_vs_pt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 50, 150, 0, 3, name=hname, title="; gen p_{T} [GeV]; ME2 k-factor", type='F')
hname = "h2_k_vs_pt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 50, 150, 0, 3, name=hname, title="; gen p_{T} [GeV]; ME3 k-factor", type='F')
hname = "h2_k_vs_pt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 50, 150, 0, 3, name=hname, title="; gen p_{T} [GeV]; ME4 k-factor", type='F')

hname = "h2_k_vs_invpt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 150, 0, 3, name=hname, title="; gen 1/p_{T} [1/GeV]; ME1 k-factor", type='F')
hname = "h2_k_vs_invpt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 150, 0, 3, name=hname, title="; gen 1/p_{T} [1/GeV]; ME2 k-factor", type='F')
hname = "h2_k_vs_invpt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 150, 0, 3, name=hname, title="; gen 1/p_{T} [1/GeV]; ME3 k-factor", type='F')
hname = "h2_k_vs_invpt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 150, 0, 3, name=hname, title="; gen 1/p_{T} [1/GeV]; ME4 k-factor", type='F')

hname = "h2_ms_vs_pt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 50, 320, -0.2, 0.2, name=hname, title="; gen p_{T} [GeV]; ME1 ms-factor", type='F')
hname = "h2_ms_vs_pt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 50, 160, -0.2, 0.2, name=hname, title="; gen p_{T} [GeV]; ME2 ms-factor", type='F')
hname = "h2_ms_vs_pt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 50, 160, -0.2, 0.2, name=hname, title="; gen p_{T} [GeV]; ME3 ms-factor", type='F')
hname = "h2_ms_vs_pt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 50, 160, -0.2, 0.2, name=hname, title="; gen p_{T} [GeV]; ME4 ms-factor", type='F')

hname = "h2_ms_vs_invpt_st1"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 320, -0.2, 0.2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME1 ms-factor", type='F')
hname = "h2_ms_vs_invpt_st2"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 160, -0.2, 0.2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME2 ms-factor", type='F')
hname = "h2_ms_vs_invpt_st3"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 160, -0.2, 0.2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME3 ms-factor", type='F')
hname = "h2_ms_vs_invpt_st4"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 160, -0.2, 0.2, name=hname, title="; gen 1/p_{T} [1/GeV]; ME4 ms-factor", type='F')

hname = "h2_dinvpt_vs_invpt_mode15"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 300, -0.3, 0.3, name=hname, title="; gen 1/p_{T} [1/GeV]; #Delta(1/p_{T})", type='F')
hname = "h2_dpt_vs_invpt_mode15"
histogram2Ds[hname] = Hist2D(50, 0, 0.5, 400, -2, 2, name=hname, title="; gen 1/p_{T} [1/GeV]; #Delta(p_{T})/p_{T}", type='F')



hnames = [
  "h2_dphi_vs_invpt_st1", "h2_dphi_vs_invpt_st2", "h2_dphi_vs_invpt_st3", "h2_dphi_vs_invpt_st4",
  "h2_k_vs_invpt_st1", "h2_k_vs_invpt_st2", "h2_k_vs_invpt_st3", "h2_k_vs_invpt_st4",
]

for fr in ['f', 'r']:
  for hname in hnames:
    h = histogram2Ds[hname]
    hname1 = h.GetName() + fr
    histogram2Ds[hname1] = h.Clone(hname1)

## CSC
#for i in xrange(4):
#  for j in xrange(5):
#    hname = "dphi_csc_st%i_eta%i" % (i,j)
#    histograms[hname] = Hist(100, -360., 360., name=hname, title="; #Delta#phi [deg]", type='F')

## TT
#for i in xrange(5):
#  for j in xrange(5):
#    hname = "dphi_tt_st%i_eta%i" % (i,j)
#    histograms[hname] = Hist(100, -360., 360., name=hname, title="; #Delta#phi [deg]", type='F')


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
  # Make plots

  no_genparticles = (len(evt.genparticles) == 0)

  if no_genparticles: continue
  assert len(evt.genparticles) == 1

  mypart = evt.genparticles[0]
  ipt = get_pt_bin(mypart.pt)
  ieta = get_eta_bin(mypart.eta)

  #if ipt == -1 or ieta == -1: continue

  # Fill histograms
  for ihit, hit in enumerate(evt.hits):
    if hit.type == kCSC and not bool(hit.neighbor):

      if hit.fr == 1:
        fr_lambda = lambda x: x + 'f'
      else:
        fr_lambda = lambda x: x + 'r'

      hit_sim_phi = deg_to_rad(hit.sim_phi)

      if hit.station == 1 and (hit.ring == 1 or hit.ring == 4):
        histogram2Ds["h2_fr_vs_pt_st1"].fill(mypart.pt, hit.fr)
        histogram2Ds["h2_fr_vs_invpt_st1"].fill(1.0/mypart.pt, hit.fr)
        histogram2Ds["h2_dphi_vs_pt_st1"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_dphi_vs_invpt_st1"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds[fr_lambda("h2_dphi_vs_invpt_st1")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_k_vs_pt_st1"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds["h2_k_vs_invpt_st1"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds[fr_lambda("h2_k_vs_invpt_st1")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
      elif hit.station == 2 and hit.ring == 1:
        histogram2Ds["h2_fr_vs_pt_st2"].fill(mypart.pt, hit.fr)
        histogram2Ds["h2_fr_vs_invpt_st2"].fill(1.0/mypart.pt, hit.fr)
        histogram2Ds["h2_dphi_vs_pt_st2"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_dphi_vs_invpt_st2"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds[fr_lambda("h2_dphi_vs_invpt_st2")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_k_vs_pt_st2"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds["h2_k_vs_invpt_st2"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds[fr_lambda("h2_k_vs_invpt_st2")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
      elif hit.station == 3 and hit.ring == 1:
        histogram2Ds["h2_fr_vs_pt_st3"].fill(mypart.pt, hit.fr)
        histogram2Ds["h2_fr_vs_invpt_st3"].fill(1.0/mypart.pt, hit.fr)
        histogram2Ds["h2_dphi_vs_pt_st3"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_dphi_vs_invpt_st3"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds[fr_lambda("h2_dphi_vs_invpt_st3")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_k_vs_pt_st3"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds["h2_k_vs_invpt_st3"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds[fr_lambda("h2_k_vs_invpt_st3")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
      elif hit.station == 4 and hit.ring == 1:
        histogram2Ds["h2_fr_vs_pt_st4"].fill(mypart.pt, hit.fr)
        histogram2Ds["h2_fr_vs_invpt_st4"].fill(1.0/mypart.pt, hit.fr)
        histogram2Ds["h2_dphi_vs_pt_st4"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_dphi_vs_invpt_st4"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds[fr_lambda("h2_dphi_vs_invpt_st4")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi))
        histogram2Ds["h2_k_vs_pt_st4"].fill(mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds["h2_k_vs_invpt_st4"].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)
        histogram2Ds[fr_lambda("h2_k_vs_invpt_st4")].fill(1.0/mypart.pt, delta_phi(hit_sim_phi, mypart.phi) * mypart.pt)


  # Multiple scattering
  # (only makes sense for the 0T sample)
  me0_sim_phi, me1_sim_phi, me2_sim_phi, me3_sim_phi, me4_sim_phi = None, None, None, None, None

  for ihit, hit in enumerate(evt.hits):
    if hit.type == kCSC and not bool(hit.neighbor):
      if me0_sim_phi is None:
        # (1/2rho) / cot(theta) * (z-z0)
        magfield = 3.8
        #magfield = 0
        me0_sim_phi = mypart.phi - 0.5 * 0.003 * magfield * (mypart.q/mypart.pt) / np.sinh(mypart.eta) * 500

      if hit.station == 1 and (hit.ring == 1 or hit.ring == 4) and hit.fr == 1 and me1_sim_phi is None:
        me1_sim_phi = deg_to_rad(hit.sim_phi)
      elif hit.station == 2 and hit.ring == 1 and hit.fr == 1 and me2_sim_phi is None:
        me2_sim_phi = deg_to_rad(hit.sim_phi)
      elif hit.station == 3 and hit.ring == 1 and hit.fr == 0 and me3_sim_phi is None:
        me3_sim_phi = deg_to_rad(hit.sim_phi)
      elif hit.station == 4 and hit.ring == 1 and hit.fr == 0 and me4_sim_phi is None:
        me4_sim_phi = deg_to_rad(hit.sim_phi)

  #print me0_sim_phi, me1_sim_phi, me2_sim_phi, me3_sim_phi, me4_sim_phi
  #print delta_phi(me1_sim_phi, me0_sim_phi), delta_phi(me2_sim_phi, me1_sim_phi)

  if me0_sim_phi is not None and me1_sim_phi is not None:
    histogram2Ds["h2_ms_vs_pt_st1"].fill(mypart.pt, delta_phi(me1_sim_phi, me0_sim_phi))
    histogram2Ds["h2_ms_vs_invpt_st1"].fill(1.0/mypart.pt, delta_phi(me1_sim_phi, me0_sim_phi))
  if me1_sim_phi is not None and me2_sim_phi is not None:
    histogram2Ds["h2_ms_vs_pt_st2"].fill(mypart.pt, delta_phi(me2_sim_phi, me1_sim_phi))
    histogram2Ds["h2_ms_vs_invpt_st2"].fill(1.0/mypart.pt, delta_phi(me2_sim_phi, me1_sim_phi))
  if me2_sim_phi is not None and me3_sim_phi is not None:
    histogram2Ds["h2_ms_vs_pt_st3"].fill(mypart.pt, delta_phi(me3_sim_phi, me2_sim_phi))
    histogram2Ds["h2_ms_vs_invpt_st3"].fill(1.0/mypart.pt, delta_phi(me3_sim_phi, me2_sim_phi))
  if me3_sim_phi is not None and me4_sim_phi is not None:
    histogram2Ds["h2_ms_vs_pt_st4"].fill(mypart.pt, delta_phi(me4_sim_phi, me3_sim_phi))
    histogram2Ds["h2_ms_vs_invpt_st4"].fill(1.0/mypart.pt, delta_phi(me4_sim_phi, me3_sim_phi))


  # pT resolution
  # (only make sense for the actual full sim sample)

  for itrk, trk in enumerate(evt.tracks):
    if trk.mode == 15:
      #residual = 1.0/trk.pt - 1.0/mypart.pt
      residual = 1.0/trk.xml_pt - 1.0/mypart.pt
      histogram2Ds["h2_dinvpt_vs_invpt_mode15"].fill(1.0/mypart.pt, residual)
      histogram2Ds["h2_dpt_vs_invpt_mode15"].fill(1.0/mypart.pt, residual * mypart.pt)
      break


  continue  # end loop over event


# ______________________________________________________________________________
# Output

with root_open('histos_trackbuilding.root', 'RECREATE') as f:
  for k, v in histograms.iteritems():
    v.Write()
  for k, v in histogram2Ds.iteritems():
    v.Write()

# ______________________________________________________________________________
# Drawer

make_plots = True

if make_plots:
  from drawer import *
  mydrawer = MyDrawer()
  options = mydrawer.options

  if True:
    hnames = [
      "h2_dphi_vs_invpt_st1f", "h2_dphi_vs_invpt_st2f", "h2_dphi_vs_invpt_st3r", "h2_dphi_vs_invpt_st4r",
    ]

    for hname in hnames:
      h = histogram2Ds[hname]
      h.Draw("COLZ")
      gPad.Print(hname+".png")
      h.RebinX(2)

      h_pfx = h.ProfileX(hname+"_pfx", 1, -1, "s")
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      #h_pfx.Fit("pol1", "", "", 0.025, 0.2499)
      h_pfx.Fit("pol1")
      gPad.Print(h_pfx.GetName()+".png")

      # Apply gaussian fits
      gr1 = TGraphAsymmErrors(h.GetNbinsX())
      gr2 = TGraphAsymmErrors(h.GetNbinsX())
      for i in xrange(h.GetNbinsX()):
        h_py = h.ProjectionY("_py", i+1, i+1)
        if h_py.Integral() < 15:  continue
        #r = h_py.Fit("gaus", "SNQ")
        r = h_py.Fit("gaus", "SNQ", "", h_py.GetMean() - 0.04*4, h_py.GetMean() + 0.04*4)
        mean, sigma, meanErr, sigmaErr = r.Parameter(1), r.Parameter(2), r.ParError(1), r.ParError(2)
        gr1.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), mean)
        gr1.SetPointError(i, 0, 0, sigma, sigma)
        gr2.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), sigma)
        gr2.SetPointError(i, 0, 0, sigmaErr, sigmaErr)

      hname1 = hname.replace("h2_dphi_vs_invpt", "h2_dphi_vs_invpt_snq1")
      h_pfx = h.ProfileX(hname1+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      gr1.Draw("p")
      gr1.Fit("pol1", "", "", 0.025, 0.2499)
      #fa1 = TF1("fa1", "[0]*x")
      #gr1.Fit("fa1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")

      hname2 = hname.replace("h2_dphi_vs_invpt", "h2_dphi_vs_invpt_snq2")
      h_pfx = h.ProfileX(hname2+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(0.1)
      h_pfx.SetMinimum(0)
      h_pfx.Draw()
      gr2.Draw("p")
      gr2.Fit("pol1", "", "", 0.025, 0.2499)
      #fa1 = TF1("fa1", "[0]")
      #gr2.Fit("fa1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")

  if True:
    # Multiple scattering
    hnames = [
      "h2_ms_vs_invpt_st1", "h2_ms_vs_invpt_st2", "h2_ms_vs_invpt_st3", "h2_ms_vs_invpt_st4",
    ]

    for hname in hnames:
      h = histogram2Ds[hname]
      h.Draw("COLZ")
      gPad.Print(hname+".png")
      h.RebinX(2)

      h_pfx = h.ProfileX(hname+"_pfx", 1, -1, "s")
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      h_pfx.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")

      # Apply gaussian fits
      gr1 = TGraphAsymmErrors(h.GetNbinsX())
      gr2 = TGraphAsymmErrors(h.GetNbinsX())
      for i in xrange(h.GetNbinsX()):
        h_py = h.ProjectionY("_py", i+1, i+1)
        if h_py.Integral() < 15:  continue
        #r = h_py.Fit("gaus", "SNQ")
        r = h_py.Fit("gaus", "SNQ", "", h_py.GetMean() - 0.04*4, h_py.GetMean() + 0.04*4)
        mean, sigma, meanErr, sigmaErr = r.Parameter(1), r.Parameter(2), r.ParError(1), r.ParError(2)
        gr1.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), mean)
        gr1.SetPointError(i, 0, 0, sigma, sigma)
        gr2.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), sigma)
        gr2.SetPointError(i, 0, 0, sigmaErr, sigmaErr)

      hname1 = hname.replace("h2_ms_vs_invpt", "h2_ms_vs_invpt_snq1")
      h_pfx = h.ProfileX(hname1+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      gr1.Draw("p")
      gr1.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")

      hname2 = hname.replace("h2_ms_vs_invpt", "h2_ms_vs_invpt_snq2")
      h_pfx = h.ProfileX(hname2+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(0.1)
      h_pfx.SetMinimum(0)
      h_pfx.Draw()
      gr2.Draw("p")
      #gr2.Fit("pol1", "", "", 0.025, 0.2499)
      fa1 = TF1("fa1", "[0]*x")
      gr2.Fit("fa1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")

  if True:
    # pT resolution
    hnames = [
      "h2_dinvpt_vs_invpt_mode15", "h2_dpt_vs_invpt_mode15",
    ]

    for hname in hnames:
      h = histogram2Ds[hname]
      h.Draw("COLZ")
      gPad.Print(hname+".png")
      h.RebinX(2)

      h_pfx = h.ProfileX(hname+"_pfx", 1, -1, "s")
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      #h_pfx.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")
      #
      # Apply gaussian fits
      gr1 = TGraphAsymmErrors(h.GetNbinsX())
      gr2 = TGraphAsymmErrors(h.GetNbinsX())
      gr1_topt = TGraphAsymmErrors(h.GetNbinsX())
      gr2_topt = TGraphAsymmErrors(h.GetNbinsX())
      for i in xrange(h.GetNbinsX()):
        h_py = h.ProjectionY("_py", i+1, i+1)
        if h_py.Integral() < 15:  continue
        #r = h_py.Fit("gaus", "SNQ")
        r = h_py.Fit("gaus", "SNQ", "", h_py.GetMean() - 0.04*4, h_py.GetMean() + 0.04*4)
        mean, sigma, meanErr, sigmaErr = r.Parameter(1), r.Parameter(2), r.ParError(1), r.ParError(2)
        gr1.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), mean)
        gr1.SetPointError(i, 0, 0, sigma, sigma)
        gr2.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), sigma)
        gr2.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
        gr1_topt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), mean)
        gr1_topt.SetPointError(i, 0, 0, sigma, sigma)
        gr2_topt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), sigma)
        gr2_topt.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
      #
      hname1 = hname.replace("h2_dinvpt_vs_invpt_mode15", "h2_dinvpt_vs_invpt_mode15_snq1").replace("h2_dpt_vs_invpt_mode15", "h2_dpt_vs_invpt_mode15_snq1")
      h_pfx = h.ProfileX(hname1+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      gr1.Draw("p")
      #gr1.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")
      #
      hname2 = hname.replace("h2_dinvpt_vs_invpt_mode15", "h2_dinvpt_vs_invpt_mode15_snq2").replace("h2_dpt_vs_invpt_mode15", "h2_dpt_vs_invpt_mode15_snq2")
      h_pfx = h.ProfileX(hname2+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetMaximum(0.1)
      h_pfx.SetMinimum(0)
      h_pfx.Draw()
      gr2.Draw("p")
      #gr2.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.Print(h_pfx.GetName()+".png")
      #
      hname1 = hname.replace("h2_dinvpt_vs_invpt_mode15", "h2_dinvpt_vs_pt_mode15_snq1").replace("h2_dpt_vs_invpt_mode15", "h2_dpt_vs_pt_mode15_snq1")
      h_pfx = h.ProfileX(hname1+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetBins(50, 0, 50)
      h_pfx.GetXaxis().SetTitle("gen p_{T} [GeV]")
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      gr1_topt.Draw("p")
      #gr1_topt.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.SetLogx(1)
      gPad.Print(h_pfx.GetName()+".png")
      gPad.SetLogx(0)
      #
      hname2 = hname.replace("h2_dinvpt_vs_invpt_mode15", "h2_dinvpt_vs_pt_mode15_snq2").replace("h2_dpt_vs_invpt_mode15", "h2_dpt_vs_pt_mode15_snq2")
      h_pfx = h.ProfileX(hname2+"_pfx", 1, -1, "s")
      h_pfx.Reset()
      h_pfx.SetBins(50, 0, 50)
      h_pfx.GetXaxis().SetTitle("gen p_{T} [GeV]")
      h_pfx.SetMaximum(1.2)
      h_pfx.SetMinimum(-0.2)
      h_pfx.Draw()
      gr2_topt.Draw("p")
      #gr2_topt.Fit("pol1", "", "", 0.025, 0.2499)
      gPad.SetLogx(1)
      gPad.Print(h_pfx.GetName()+".png")
      gPad.SetLogx(0)

