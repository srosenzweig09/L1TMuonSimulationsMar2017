"""Utilities for EMTF++."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import random

from six import iteritems
from six.moves import range, zip, map, filter

# ______________________________________________________________________________
# Utilities

# Enums
kDT, kCSC, kRPC, kGEM, kME0 = 0, 1, 2, 3, 4

kDEBUG, kINFO, kWARNING, kERROR, kFATAL = 0, 1, 2, 3, 4

# Functions
def wrap_phi_rad(rad):
  while rad < -np.pi:
    rad += np.pi*2
  while rad >= +np.pi:
    rad -= np.pi*2
  return rad

def wrap_phi_deg(deg):
  while deg < -180.:
    deg += 360.
  while deg >= +180.:
    deg -= 360.
  return deg

def wrap_theta_rad(rad):
  rad = np.abs(rad)
  while rad >= np.pi:
    rad -= np.pi
  if rad >= np.pi/2:
    rad = np.pi - rad
  return rad

def wrap_theta_deg(deg):
  deg = np.abs(deg)
  while deg >= 180.:
    deg -= 180.
  if deg >= 180./2:
    deg = 180. - deg
  return deg

def delta_phi(lhs, rhs):  # in radians
  rad = lhs - rhs
  rad = wrap_phi_rad(rad)
  return rad

def delta_theta(lhs, rhs):  # in radians
  rad = lhs - rhs
  return rad

def calc_phi_loc_deg_from_glob(glob, sector):
  # glob in deg, sector [1-6]
  glob = wrap_phi_deg(glob)
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

def calc_phi_loc_deg(bits):
  loc = float(bits) / 60. - 22.
  return loc

def calc_phi_glob_deg(loc, sector):
  # loc in deg, sector [1-6]
  glob = loc + 15. + (60. * (sector-1))
  if glob >= 180.:
    glob -= 360.
  return glob

def calc_theta_int(theta, endcap):
  # theta in deg, endcap [-1,+1]
  if endcap == -1:
    theta = 180. - theta
  theta = (theta - 8.5) * 128. / (45.0-8.5)
  theta_int = int(round(theta))
  return theta_int

def calc_theta_rad_from_eta(eta):
  # returns theta in [0-pi] rad
  theta = np.arctan2(1.0, np.sinh(eta))
  return theta

def calc_theta_deg_from_eta(eta):
  # returns theta in [0-180] deg
  return np.rad2deg(calc_theta_rad_from_eta(eta))

def calc_theta_deg_from_int(theta_int):
  theta_deg = float(theta_int) * (45.0-8.5) / 128. + 8.5;
  return theta_deg

def calc_eta_from_theta_rad(theta_rad):
  eta = -1. * np.log(np.tan(theta_rad/2.))
  return eta

def calc_eta_from_theta_deg(theta_deg, endcap):
  # theta in deg, endcap [-1,+1]
  theta_deg = wrap_theta_deg(theta_deg)
  theta_rad = np.deg2rad(theta_deg)
  eta = calc_eta_from_theta_rad(theta_rad)
  if endcap == -1:
    eta = -eta
  return eta

def find_endsec(endcap, sector):
  endsec = (sector - 1) if endcap == 1 else (sector - 1 + 6)
  return endsec

def pick_the_median(lst):  # assume sorted list
  middle = 0 if len(lst) == 0 else (len(lst)-1)//2
  return lst[middle]

def save_np_arrays(outfile, outdict):
  from numpy.compat import contextlib_nullcontext
  with contextlib_nullcontext(outfile) as f:
    np.savez_compressed(f, **outdict)

def save_root_histograms(outfile, histograms):
  from rootpy.io import root_open
  with root_open(outfile, 'recreate') as f:
    for (k, v) in histograms.iteritems():
      v.Write()

# Based on
#   https://www.tensorflow.org/guide/ragged_tensor
#   https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/ops/ragged/ragged_tensor_value.py
class RaggedTensorValue(object):
  """Represents the value of a `RaggedTensor`."""

  def __init__(self, values, row_splits):
    """Creates a `RaggedTensorValue`.
    Args:
      values: A numpy array of any type and shape; or a RaggedTensorValue.
      row_splits: A 1-D int32 or int64 numpy array.
    """
    if not (isinstance(row_splits, (np.ndarray, np.generic)) and
            row_splits.dtype in (np.int64, np.int32) and row_splits.ndim == 1):
      raise TypeError("row_splits must be a 1D int32 or int64 numpy array")
    if not isinstance(values, (np.ndarray, np.generic, RaggedTensorValue)):
      raise TypeError("values must be a numpy array or a RaggedTensorValue")
    if (isinstance(values, RaggedTensorValue) and
        row_splits.dtype != values.row_splits.dtype):
      raise ValueError("row_splits and values.row_splits must have "
                       "the same dtype")
    self._values = values
    self._row_splits = row_splits

  row_splits = property(
      lambda self: self._row_splits,
      doc="""The split indices for the ragged tensor value.""")
  values = property(
      lambda self: self._values,
      doc="""The concatenated values for all rows in this tensor.""")
  dtype = property(
      lambda self: self._values.dtype,
      doc="""The numpy dtype of values in this tensor.""")

  @property
  def shape(self):
    """A tuple indicating the shape of this RaggedTensorValue."""
    return (self._row_splits.shape[0] - 1,) + (None,) + self._values.shape[1:]

  def to_list(self):
    """Returns this ragged tensor value as a nested Python list."""
    if isinstance(self._values, RaggedTensorValue):
      values_as_list = self._values.to_list()
    else:
      values_as_list = self._values.tolist()
    return [
        values_as_list[self._row_splits[i]:self._row_splits[i + 1]]
        for i in range(self._row_splits.shape[0] - 1)
    ]

  def to_array(self):
    """Returns this ragged tensor value as a nested Numpy array."""
    arr = np.empty((self._row_splits.shape[0] - 1,), dtype=np.object)
    for i in range(self._row_splits.shape[0] - 1):
      arr[i] = self._values[self._row_splits[i]:self._row_splits[i + 1]]
    return arr

# Based on
#   https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/ops/ragged/ragged_factory_ops.py
def create_ragged_array(pylist):
  """Construct a constant RaggedTensorValue from a nested list."""

  # Ragged rank for returned value
  ragged_rank = 1

  # Build the splits for each ragged rank, and concatenate the inner values
  # into a single list.
  nested_splits = []
  values = pylist
  for dim in range(ragged_rank):
    nested_splits.append([0])
    concatenated_values = []
    for row in values:
      nested_splits[dim].append(nested_splits[dim][-1] + len(row))
      concatenated_values.extend(row)
    values = concatenated_values

  values = np.asarray(values)
  for row_splits in reversed(nested_splits):
    row_splits = np.asarray(row_splits, dtype=np.int32)
    values = RaggedTensorValue(values, row_splits)
  return values


# ______________________________________________________________________________
# Models



# ______________________________________________________________________________
# Datasets

def load_tree(infiles):
  from rootpy.tree import TreeChain
  from rootpy.ROOT import gROOT
  gROOT.SetBatch(True)

  print('[INFO] Opening files: {0}'.format(infiles))
  tree = TreeChain('ntupler/tree', infiles)
  tree.define_collection(name='hits', prefix='vh_', size='vh_size')
  tree.define_collection(name='simhits', prefix='vc_', size='vc_size')
  tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
  tree.define_collection(name='particles', prefix='vp_', size='vp_size')
  return tree

eos_prefix = 'root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_10_6_3/'

def load_pgun_test():
  #infile = '../test7/ntuple_SingleMuon_Endcap_2GeV_add.6.root'
  infile = '../test7/ntuple_SingleMuon_Displaced_2GeV_PhaseIITDRSpring19.4.root'
  return load_tree(infile)

def load_pgun_batch(k):
  my_range = np.split(np.arange(1000), 100)[k]
  infiles = []
  for i in my_range:
    infiles.append(eos_prefix + 'SingleMuon_Endcap_2GeV_PhaseIITDRSpring19/ParticleGuns/CRAB3/190923_203236/%04i/ntuple_SingleMuon_Endcap_%i.root' % ((i+1)//1000, (i+1)))
    infiles.append(eos_prefix + 'SingleMuon_Endcap2_2GeV_PhaseIITDRSpring19/ParticleGuns/CRAB3/190923_203346/%04i/ntuple_SingleMuon_Endcap2_%i.root' % ((i+1)//1000, (i+1)))
  tree = load_tree(infiles)
  return tree

def load_pgun_displ_batch(k):
  my_range = np.split(np.arange(1000), 100)[k]
  infiles = []
  for i in my_range:
    infiles.append(eos_prefix + 'SingleMuon_Displaced_2GeV_PhaseIITDRSpring19/ParticleGuns/CRAB3/190923_212343/%04i/ntuple_SingleMuon_Displaced_%i.root' % ((i+1)//1000, (i+1)))
    infiles.append(eos_prefix + 'SingleMuon_Displaced2_2GeV_PhaseIITDRSpring19/ParticleGuns/CRAB3/190923_212452/%04i/ntuple_SingleMuon_Displaced2_%i.root' % ((i+1)//1000, (i+1)))
  tree = load_tree(infiles)
  return tree

def load_mixing_batch(k):
  infiles = []
  infiles += [eos_prefix + 'ntuple_SingleNeutrino_PU140_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145646/0000/ntuple_%i.root' % (i+1) for i in range(30)]  # up to 30/63
  infiles += [eos_prefix + 'ntuple_SingleNeutrino_PU200_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145529/0000/ntuple_%i.root' % (i+1) for i in range(40)]  # up to 40/85
  #infiles += [eos_prefix + 'ntuple_SingleNeutrino_PU250_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145757/0000/ntuple_%i.root' % (i+1) for i in range(50)]  # up to 50/125
  #infiles += [eos_prefix + 'ntuple_SingleNeutrino_PU300_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/191002_214457/0000/ntuple_%i.root' % (i+1) for i in range(50)]  # up to 50/111
  infiles += [eos_prefix + 'ntuple_DoubleElectron_PU140_PhaseIITDRSpring19/DoubleElectron_FlatPt-1To100/CRAB3/190925_044352/0000/ntuple_%i.root' % (i+1) for i in range(25)]
  infiles += [eos_prefix + 'ntuple_DoubleElectron_PU200_PhaseIITDRSpring19/DoubleElectron_FlatPt-1To100/CRAB3/190925_044710/0000/ntuple_%i.root' % (i+1) for i in range(25)]
  infiles += [eos_prefix + 'ntuple_DoublePhoton_PU140_PhaseIITDRSpring19/DoublePhoton_FlatPt-1To100/CRAB3/190925_044839/0000/ntuple_%i.root' % (i+1) for i in range(25)]
  infiles += [eos_prefix + 'ntuple_DoublePhoton_PU200_PhaseIITDRSpring19/DoublePhoton_FlatPt-1To100/CRAB3/190925_044947/0000/ntuple_%i.root' % (i+1) for i in range(25)]
  infiles += [eos_prefix + 'ntuple_SingleElectron_PU200_PhaseIITDRSpring19/SingleElectron_PT2to100/CRAB3/190925_181236/0000/ntuple_%i.root' % (i+1) for i in range(15)]
  infiles += [eos_prefix + 'ntuple_SinglePhoton_PU200_PhaseIITDRSpring19/PhotonFlatPt8To150/CRAB3/190925_181357/0000/ntuple_%i.root' % (i+1) for i in range(77)]
  infile = infiles[k]
  tree = load_tree(infile)
  return tree

def load_singleneutrino_batch(k, pileup=200):
  if pileup == 140:
    infiles = [eos_prefix + 'ntuple_SingleNeutrino_PU140_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145646/0000/ntuple_%i.root' % (i+1) for i in range(63)]
  elif pileup == 200:
    infiles = [eos_prefix + 'ntuple_SingleNeutrino_PU200_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145529/0000/ntuple_%i.root' % (i+1) for i in range(85)]
  elif pileup == 250:
    infiles = [eos_prefix + 'ntuple_SingleNeutrino_PU250_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/190926_145757/0000/ntuple_%i.root' % (i+1) for i in range(125)]
  elif pileup == 300:
    infiles = [eos_prefix + 'ntuple_SingleNeutrino_PU300_PhaseIITDRSpring19/Nu_E10-pythia8-gun/CRAB3/191002_214457/0000/ntuple_%i.root' % (i+1) for i in range(111)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  infile = infiles[k]
  tree = load_tree(infile)
  return tree

def load_mumu_flatpt_batch(k, pileup=200):
  if pileup == 0:
    infiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU0_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190925_042003/0000/ntuple_MuMu_FlatPt_PU0_%i.root' % (i+1) for i in range(5)]
  elif pileup == 200:
    infiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU200_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190925_051735/0000/ntuple_%i.root' % (i+1) for i in range(33)]
  elif pileup == 300:
    infiles = [eos_prefix + 'ntuple_MuMu_FlatPt_PU300_PhaseIITDRSpring19/Mu_FlatPt2to100-pythia8-gun/CRAB3/190924_201214/0000/ntuple_MuMu_FlatPt_PU300_%i.root' % (i+1) for i in range(280)]
  else:
    raise RuntimeError('Cannot recognize pileup: {0}'.format(pileup))
  infile = infiles[k]
  tree = load_tree(infile)
  return tree
