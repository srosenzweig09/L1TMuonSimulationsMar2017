import numpy as np
np.random.seed(2026)

import os, sys, datetime
from six.moves import range, zip, map, filter

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency
from rootpy.tree import Tree, TreeChain
from rootpy.io import root_open
from rootpy.ROOT import gROOT
gROOT.SetBatch(True)

# Adjust matplotlib logging
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


# ______________________________________________________________________________
# Utilities

def calculate_d0(invPt, phi, xv, yv, B=3.811):
  _invPt = np.asarray(invPt, dtype=np.float64)   # needs double precision
  _invPt = np.where(np.abs(_invPt) < 1./10000, np.sign(_invPt+1e-15) * 1./10000, _invPt)
  _R = -1.0 / (0.003 * B * _invPt)               # R = -pT/(0.003 q B)  [cm]
  _xc = xv - (_R * np.sin(phi))                  # xc = xv - R sin(phi)
  _yc = yv + (_R * np.cos(phi))                  # yc = yv + R cos(phi)
  _d0 = _R - (np.sign(_R) * np.hypot(_xc, _yc))  # d0 = R - sign(R) * sqrt(xc^2 + yc^2)
  return _d0.astype(np.float32, casting='same_kind')


# ______________________________________________________________________________
# Analysis: effie

class EffieAnalysis(object):
  def run(self, omtf_input=False, run2_input=False, pileup=200):
    # Book histograms
    histograms = {}
    eff_pt_bins = (0., 0.5, 1., 1.5, 2., 3., 4., 5., 6., 7., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 34., 40., 48., 60., 80., 100., 120.)
    eff_highpt_bins = (2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 34., 40., 48., 60., 80., 100., 120., 250., 500., 1000.)

    for m in ("emtf", "emtf2026"):
      for l in (0, 5, 10, 20, 30, 40, 50, 60):
        for k in ("denom", "numer"):
          hname = "%s_eff_vs_genpt_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(eff_pt_bins, name=hname, title="; gen p_{T} [GeV]", type='F')
          hname = "%s_eff_vs_genpt_highpt_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(eff_highpt_bins, name=hname, title="; gen p_{T} [GeV]", type='F')
          hname = "%s_eff_vs_genphi_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(76, -190, 190, name=hname, title="; gen #phi {gen p_{T} > 20 GeV}", type='F')
          hname = "%s_eff_vs_genphi_lowpt_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(76, -190, 190, name=hname, title="; gen #phi {5 < gen p_{T} #leq 20 GeV}", type='F')
          hname = "%s_eff_vs_geneta_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(26, 1.2, 2.5, name=hname, title="; gen |#eta| {gen p_{T} > 20 GeV}", type='F')
          hname = "%s_eff_vs_geneta_lowpt_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(26, 1.2, 2.5, name=hname, title="; gen |#eta| {5 < gen p_{T} #leq 20 GeV}", type='F')
          hname = "%s_eff_vs_gend0_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(80, 0, 120, name=hname, title="; gen |d_{0}| [cm] {gen p_{T} > 20 GeV}", type='F')
          hname = "%s_eff_vs_gendz_l1pt%i_%s" % (m,l,k)
          histograms[hname] = Hist(80, 0, 40, name=hname, title="; gen |d_{z}| [cm] {gen p_{T} > 20 GeV}", type='F')

    # Event range
    maxEvents = -1

    # __________________________________________________________________________
    # Loop over events (EMTF)
    tree = load_minbias_batch_for_effie(jobid, pileup=pileup)

    for ievt, evt in enumerate(tree):
      if maxEvents != -1 and ievt == maxEvents:
        break

      if len(evt.particles) == 0:
        continue

      part = evt.particles[0]  # particle gun
      if sum([(part.status == 1) for part in evt.particles]) == 2:
        part = evt.particles[(ievt % 2)]  # centrally produced samples contain 2 muons, pick only one
      part.invpt = np.true_divide(part.q, part.pt)
      part.d0 = calculate_d0(part.invpt, part.phi, part.vx, part.vy)
      if part.eta >= 0.:
        sectors = [0, 1, 2, 3, 4, 5]
      else:
        sectors = [6, 7, 8, 9, 10, 11]

      # ________________________________________________________________________
      # Fill histograms

      def base_fill_efficiency(harvest_part, select_part, select_track):
        if select_part(part):
          trigger = any([select_track(trk) for trk in tracks])  # using scaled pT
          denom = histograms[hname + "_denom"]
          numer = histograms[hname + "_numer"]
          denom.fill(harvest_part(part))
          if trigger:
            numer.fill(harvest_part(part))

      def fill_efficiency_pt():
        harvest_part = lambda part: part.pt
        base_fill_efficiency(harvest_part, select_part, select_track)

      def fill_efficiency_phi():
        harvest_part = lambda part: np.rad2deg(part.phi)
        base_fill_efficiency(harvest_part, select_part, select_track)

      def fill_efficiency_eta():
        harvest_part = lambda part: abs(part.eta)
        select_part_eta = lambda part: (part.bx == 0) and (abs(part.d0) >= 0)  # no cut on eta
        base_fill_efficiency(harvest_part, select_part_eta, select_track)

      def fill_efficiency_d0():
        harvest_part = lambda part: abs(part.d0)
        select_part_d0 = lambda part: (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4)  # no cut on d0
        base_fill_efficiency(harvest_part, select_part_d0, select_track)

      def fill_efficiency_dz():
        harvest_part = lambda part: abs(part.vz)
        base_fill_efficiency(harvest_part, select_part, select_track)

      # Check various L1 pT thresholds
      for l in (0, 5, 10, 20, 30, 40, 50, 60):
        # EMTF tracks
        tracks = evt.tracks
        select_part = lambda part: (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4) and (abs(part.d0) >= 0)
        select_track = lambda trk: trk and (1.24 <= abs(trk.eta) <= 2.4) and (trk.bx == 0) and (trk.mode in (11,13,14,15)) and (trk.pt > float(l)) and (trk.eta * part.eta > 0)

        hname = "emtf_eff_vs_genpt_l1pt%i" % (l)
        fill_efficiency_pt()
        hname = "emtf_eff_vs_genpt_highpt_l1pt%i" % (l)
        fill_efficiency_pt()
        if part.pt > 20.:
          hname = "emtf_eff_vs_genphi_l1pt%i" % (l)
          fill_efficiency_phi()
          hname = "emtf_eff_vs_geneta_l1pt%i" % (l)
          fill_efficiency_eta()
          hname = "emtf_eff_vs_gend0_l1pt%i" % (l)
          fill_efficiency_d0()
          hname = "emtf_eff_vs_gendz_l1pt%i" % (l)
          fill_efficiency_dz()
        elif 5 < part.pt <= 20.:
          hname = "emtf_eff_vs_genphi_lowpt_l1pt%i" % (l)
          fill_efficiency_phi()
          hname = "emtf_eff_vs_geneta_lowpt_l1pt%i" % (l)
          fill_efficiency_eta()

    # End loop over events
    unload_tree()

    # __________________________________________________________________________
    # Loop over events (EMTF++)
    tree2026 = load_minbias_batch_for_effie_emtf0(jobid, pileup=pileup)

    for ievt, evt in enumerate(tree2026):
      if maxEvents != -1 and ievt == maxEvents:
        break

      if len(evt.particles) == 0:
        continue

      part = evt.particles[0]  # particle gun
      if sum([(part.status == 1) for part in evt.particles]) == 2:
        part = evt.particles[(ievt % 2)]  # centrally produced samples contain 2 muons, pick only one
      part.invpt = np.true_divide(part.q, part.pt)
      part.d0 = calculate_d0(part.invpt, part.phi, part.vx, part.vy)
      if part.eta >= 0.:
        sectors = [0, 1, 2, 3, 4, 5]
      else:
        sectors = [6, 7, 8, 9, 10, 11]

      # ________________________________________________________________________
      # Fill histograms

      def base_fill_efficiency(harvest_part, select_part, select_track):
        if select_part(part):
          trigger = any([select_track(trk) for trk in tracks])  # using scaled pT
          denom = histograms[hname + "_denom"]
          numer = histograms[hname + "_numer"]
          denom.fill(harvest_part(part))
          if trigger:
            numer.fill(harvest_part(part))

      def fill_efficiency_pt():
        harvest_part = lambda part: part.pt
        base_fill_efficiency(harvest_part, select_part, select_track)

      def fill_efficiency_phi():
        harvest_part = lambda part: np.rad2deg(part.phi)
        base_fill_efficiency(harvest_part, select_part, select_track)

      def fill_efficiency_eta():
        harvest_part = lambda part: abs(part.eta)
        select_part_eta = lambda part: (part.bx == 0) and (abs(part.d0) >= 0)  # no cut on eta
        base_fill_efficiency(harvest_part, select_part_eta, select_track)

      def fill_efficiency_d0():
        harvest_part = lambda part: abs(part.d0)
        select_part_d0 = lambda part: (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4)  # no cut on d0
        base_fill_efficiency(harvest_part, select_part_d0, select_track)

      def fill_efficiency_dz():
        harvest_part = lambda part: abs(part.vz)
        base_fill_efficiency(harvest_part, select_part, select_track)

      # Check various L1 pT thresholds
      for l in (0, 5, 10, 20, 30, 40, 50, 60):
        # EMTF++ tracks
        tracks = evt.tracks
        select_part = lambda part: (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4) and (abs(part.d0) >= 0)
        select_track = lambda trk: trk and (1.24 <= abs(trk.eta) <= 2.4) and (trk.pt > float(l)) and (trk.eta * part.eta > 0)
        #select_track = lambda trk: trk and (1.1 <= abs(trk.eta) <= 2.6) and (abs(1.0/trk.y_displ) > float(l)) and (abs(trk.d0_displ) >= 0)

        hname = "emtf2026_eff_vs_genpt_l1pt%i" % (l)
        fill_efficiency_pt()
        hname = "emtf2026_eff_vs_genpt_highpt_l1pt%i" % (l)
        fill_efficiency_pt()
        if part.pt > 20.:
          hname = "emtf2026_eff_vs_genphi_l1pt%i" % (l)
          fill_efficiency_phi()
          hname = "emtf2026_eff_vs_geneta_l1pt%i" % (l)
          fill_efficiency_eta()
          hname = "emtf2026_eff_vs_gend0_l1pt%i" % (l)
          fill_efficiency_d0()
          hname = "emtf2026_eff_vs_gendz_l1pt%i" % (l)
          fill_efficiency_dz()
        elif 5 < part.pt <= 20.:
          hname = "emtf2026_eff_vs_genphi_lowpt_l1pt%i" % (l)
          fill_efficiency_phi()
          hname = "emtf2026_eff_vs_geneta_lowpt_l1pt%i" % (l)
          fill_efficiency_eta()

    # End loop over events
    unload_tree()

    # __________________________________________________________________________
    # Save histograms
    outfile = 'histos_tdr_effie.root'
    if use_condor:
      outfile = outfile[:-5] + ('_%i.root' % jobid)
    print('[INFO] Creating file: %s' % outfile)
    with root_open(outfile, 'recreate') as f:
      for (k, v) in histograms.iteritems():
        v.Write()


# ______________________________________________________________________________
# Analysis: rates

class RatesAnalysis(object):
  def run(self, omtf_input=False, run2_input=False, pileup=200):
    # Book histograms
    histograms = {}
    hname = "nevents"
    histograms[hname] = Hist(5, 0, 5, name=hname, title="; count", type='F')
    hname = "nevents_emtf2026"
    histograms[hname] = Hist(5, 0, 5, name=hname, title="; count", type='F')

    for m in ("emtf", "emtf2026"):
      hname = "highest_%s_absEtaMin1.24_absEtaMax2.4_qmin12_pt" % m
      histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
      hname = "highest_%s_absEtaMin1.24_absEtaMax1.65_qmin12_pt" % m
      histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
      hname = "highest_%s_absEtaMin1.65_absEtaMax2.15_qmin12_pt" % m
      histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')
      hname = "highest_%s_absEtaMin2.15_absEtaMax2.4_qmin12_pt" % m
      histograms[hname] = Hist(100, 0., 100., name=hname, title="; p_{T} [GeV]; entries", type='F')

    # Event range
    maxEvents = -1

    # __________________________________________________________________________
    # Loop over events (EMTF)
    tree = load_minbias_batch(jobid, pileup=pileup)

    for ievt, evt in enumerate(tree):
      if maxEvents != -1 and ievt == maxEvents:
        break

      # ________________________________________________________________________
      # Fill histograms
      histograms["nevents"].fill(1.0)

      def fill_highest_pt():
        highest_pt = -999999.
        for itrk, trk in enumerate(tracks):
          if select(trk):
            if highest_pt < trk.pt:  # using scaled pT
              highest_pt = trk.pt
        if highest_pt > 0.:
          highest_pt = min(100.-1e-4, highest_pt)
          histograms[hname].fill(highest_pt)

      # EMTF tracks
      tracks = evt.tracks
      select = lambda trk: trk and (1.24 <= abs(trk.eta) <= 2.4) and (trk.bx == 0) and (trk.mode in (11,13,14,15))
      hname = "highest_emtf_absEtaMin1.24_absEtaMax2.4_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (1.24 <= abs(trk.eta) < 1.65) and (trk.bx == 0) and (trk.mode in (11,13,14,15))
      hname = "highest_emtf_absEtaMin1.24_absEtaMax1.65_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (1.65 <= abs(trk.eta) < 2.15) and (trk.bx == 0) and (trk.mode in (11,13,14,15))
      hname = "highest_emtf_absEtaMin1.65_absEtaMax2.15_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (2.15 <= abs(trk.eta) <= 2.4) and (trk.bx == 0) and (trk.mode in (11,13,14,15))
      hname = "highest_emtf_absEtaMin2.15_absEtaMax2.4_qmin12_pt"
      fill_highest_pt()

    # End loop over events
    unload_tree()

    # __________________________________________________________________________
    # Loop over events (EMTF++)
    tree2026 = load_minbias_batch_emtf0(jobid, pileup=pileup)

    for ievt, evt in enumerate(tree2026):
      if maxEvents != -1 and ievt == maxEvents:
        break

      # ________________________________________________________________________
      # Fill histograms
      histograms["nevents_emtf2026"].fill(1.0)

      def fill_highest_pt():
        highest_pt = -999999.
        for itrk, trk in enumerate(tracks):
          if select(trk):
            if highest_pt < trk.pt:  # using scaled pT
              highest_pt = trk.pt
        if highest_pt > 0.:
          highest_pt = min(100.-1e-4, highest_pt)
          histograms[hname].fill(highest_pt)

      # EMTF++ tracks
      tracks = evt.tracks
      select = lambda trk: trk and (1.24 <= abs(trk.eta) <= 2.4)
      hname = "highest_emtf2026_absEtaMin1.24_absEtaMax2.4_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (1.24 <= abs(trk.eta) <= 1.65)
      hname = "highest_emtf2026_absEtaMin1.24_absEtaMax1.65_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (1.65 <= abs(trk.eta) <= 2.15)
      hname = "highest_emtf2026_absEtaMin1.65_absEtaMax2.15_qmin12_pt"
      fill_highest_pt()
      select = lambda trk: trk and (2.15 <= abs(trk.eta) <= 2.4)
      hname = "highest_emtf2026_absEtaMin2.15_absEtaMax2.4_qmin12_pt"
      fill_highest_pt()

    # End loop over events
    unload_tree()

    # __________________________________________________________________________
    # Save histograms
    outfile = 'histos_tdr_rates.root'
    if use_condor:
      outfile = outfile[:-5] + ('_%i.root' % jobid)
    print('[INFO] Creating file: %s' % outfile)
    with root_open(outfile, 'recreate') as f:
      for (k, v) in histograms.iteritems():
        v.Write()


# ______________________________________________________________________________
# Settings

# Condor or not
# if 'CONDOR_EXEC' is defined, take 3 arguments (algo, analysis, jobid)
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  os.environ['ROOTPY_GRIDMODE'] = 'true'

# Algorithm (pick one)
algo = 'default'  # phase 2
#algo = 'run3'
#algo = 'omtf'
if use_condor:
  algo = sys.argv[1]

# Analysis mode (pick one)
#analysis = 'dummy'
#analysis = 'roads'
#analysis = 'rates'
analysis = 'effie'
#analysis = 'mixing'
#analysis = 'collusion'
#analysis = 'augmentation'
#analysis = 'images'
if use_condor:
  analysis = sys.argv[2]

# Job id
jobid = 0
if use_condor:
  jobid = int(sys.argv[3])

# Pattern bank
bankfile = 'pattern_bank_18patt.29.npz'

# NN models
kerasfile = ['model.29.json', 'model_weights.29.h5',
             'model_run3.29.json', 'model_run3_weights.29.h5',
             'model_omtf.29.json', 'model_omtf_weights.29.h5',]


# ______________________________________________________________________________
# Input files

infile_r = None  # input file handle

eos_prefix = 'root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_6_3/'

def purge_bad_files(infiles):
  good_files = []
  for infile in infiles:
    try:
      _ = TreeChain('ntupler/tree', infile)
      good_files.append(infile)
    except:
      pass
  return good_files

def define_collections(tree):
  tree.define_collection(name='hits', prefix='vh_', size='vh_size')
  tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
  tree.define_collection(name='particles', prefix='vp_', size='vp_size')
  tree.define_collection(name='evt_info', prefix='ve_', size='ve_size')
  return

def load_tree_single(infile):
  print('[INFO] Opening file: %s' % infile)
  global infile_r
  infile_r = root_open(infile)
  tree = infile_r.ntupler.tree
  define_collections(tree)
  return tree

def load_tree_multiple(infiles):
  print('[INFO] Opening file: %s' % ' '.join(infiles))
  tree = TreeChain('ntupler/tree', infiles)
  define_collections(tree)
  return tree

def unload_tree():
  global infile_r
  try:
    infile_r.close()
  except:
    pass

def load_minbias_batch(k, pileup=200):
  if pileup == 140:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU140_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145646/0000/ntuple_%i.root' % (i+1) for i in range(63)]
  elif pileup == 200:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU200_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145529/0000/ntuple_%i.root' % (i+1) for i in range(85)]
  elif pileup == 250:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU250_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145757/0000/ntuple_%i.root' % (i+1) for i in range(125)]
  elif pileup == 300:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU300_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/191002_214457/0000/ntuple_%i.root' % (i+1) for i in range(111)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  #
  infile = pufiles[k]
  return load_tree_single(infile)

def load_minbias_batch_emtf0(k, pileup=200):
  if pileup == 140:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU140_PhaseIITDRSpring19_emtf0/Nu_E10-pythia8-gun/CRAB3/191003_233321/0000/ntuple_%i.root' % (i+1) for i in range(63)]
  elif pileup == 200:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU200_PhaseIITDRSpring19_emtf0/Nu_E10-pythia8-gun/CRAB3/191003_020637/0000/ntuple_%i.root' % (i+1) for i in range(85)]
  elif pileup == 250:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU250_PhaseIITDRSpring19_emtf0/Nu_E10-pythia8-gun/CRAB3/191003_233458/0000/ntuple_%i.root' % (i+1) for i in range(125)]
  elif pileup == 300:
    pufiles = [eos_prefix + 'ntuple_SingleNeutrino_PU300_PhaseIITDRSpring19_emtf0/Nu_E10-pythia8-gun/CRAB3/191007_195557/0000/ntuple_%i.root' % (i+1) for i in range(111)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  #
  infile = pufiles[k]
  return load_tree_single(infile)

def load_minbias_batch_for_effie(k, pileup=200):
  if pileup == 0:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU0_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190925_042003/0000/ntuple_MuMu_FlatPt_PU0_%i.root' % (i+1) for i in range(5)]
  elif pileup == 200:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU200_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190925_051735/0000/ntuple_%i.root' % (i+1) for i in range(33)]
  elif pileup == 300:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU300_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190924_201214/0000/ntuple_MuMu_FlatPt_PU300_%i.root' % (i+1) for i in range(280)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  #
  infile = pufiles[k]
  return load_tree_single(infile)

def load_minbias_batch_for_effie_emtf0(k, pileup=200):
  if pileup == 0:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU0_PhaseIITDRSpring19_emtf0/Mu_FlatPt2to100-pythia8-gun/CRAB3/000000_000000/0000/ntuple_%i.root' % (i+1) for i in range(5)]
  elif pileup == 200:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU200_PhaseIITDRSpring19_emtf0/Mu_FlatPt2to100-pythia8-gun/CRAB3/191003_233638/0000/ntuple_%i.root' % (i+1) for i in range(33)]
  elif pileup == 300:
    pufiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU300_PhaseIITDRSpring19_emtf0/Mu_FlatPt2to100-pythia8-gun/CRAB3/000000_000000/0000/ntuple_%i.root' % (i+1) for i in range(280)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  #
  infile = pufiles[k]
  return load_tree_single(infile)


# ______________________________________________________________________________
# Main

if __name__ == "__main__":
  start_time = datetime.datetime.now()
  print('[INFO] Current time    : {0}'.format(start_time))
  print('[INFO] Using cmssw     : {0}'.format(os.environ['CMSSW_VERSION']))
  print('[INFO] Using condor    : {0}'.format(use_condor))
  print('[INFO] Using algo      : {0}'.format(algo))
  print('[INFO] Using analysis  : {0}'.format(analysis))
  print('[INFO] Using job id    : {0}'.format(jobid))

  if algo == 'run3':
    run2_input = True
  else:
    run2_input = False

  if algo == 'omtf':
    omtf_input = True
  else:
    omtf_input = False

  if analysis == 'effie':  # default to pileup=200
    analysis = EffieAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input)
  elif analysis == 'effie0':
    analysis = EffieAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=0)
  elif analysis == 'effie200':
    analysis = EffieAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=200)
  elif analysis == 'effie300':
    analysis = EffieAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=300)

  elif analysis == 'rates':  # default to pileup=200
    analysis = RatesAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input)
  elif analysis == 'rates140':
    analysis = RatesAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=140)
  elif analysis == 'rates200':
    analysis = RatesAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=200)
  elif analysis == 'rates250':
    analysis = RatesAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=250)
  elif analysis == 'rates300':
    analysis = RatesAnalysis()
    analysis.run(omtf_input=omtf_input, run2_input=run2_input, pileup=300)

  else:
    raise RuntimeError('Cannot recognize analysis: {0}'.format(analysis))

  stop_time = datetime.datetime.now()
  print('[INFO] Elapsed time    : {0}'.format(stop_time - start_time))
  # DONE!
