"""Data preparation."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtfpp_algos import *


# ______________________________________________________________________________
# Classes

def find_sector_rank(alist):
  rank = np.int32(0)
  for hit in alist:
    hit_lay = find_emtf_layer(hit)
    a, b, c = hit_lay//100, (hit_lay//10)%10, hit_lay%10
    if a == 5:
      rank |= (1 << 10)
    elif a == 0 and b == 0:
      rank |= (1 << 9)
    elif a == 1 and b == 0:
      rank |= (1 << 8)
    elif a == 2 and b == 0:
      rank |= (1 << 7)
    elif a == 3 and b == 0:
      rank |= (1 << 6)
    elif a == 4 and b == 0:
      rank |= (1 << 5)
    elif a == 0 and b == 1:
      rank |= (1 << 4)
    elif a == 1 and b == 1:
      rank |= (1 << 3)
    elif a == 2 and b == 1:
      rank |= (1 << 2)
    elif a == 3 and b == 1:
      rank |= (1 << 1)
    elif a == 4 and b == 1:
      rank |= (1 << 0)
  return rank

class EMTFSectorRanking(object):
  def __init__(self, find_sector_rank=find_sector_rank):
    num_sectors = 12
    self.sectors = np.empty(num_sectors, dtype=np.object)
    for ind in np.ndindex(self.sectors.shape):
      self.sectors[ind] = []
    self.find_sector_rank = find_sector_rank

  def reset(self):
    for ind in np.ndindex(self.sectors.shape):
      del self.sectors[ind][:]

  def append(self, hit):
    endsec = find_endsec(hit.endcap, hit.sector)
    self.sectors[endsec].append(hit)

  def get_best_sector(self):
    ranks = np.zeros(self.sectors.shape, dtype=np.int32)
    for ind in np.ndindex(self.sectors.shape):
      ranks[ind] = self.find_sector_rank(self.sectors[ind])
    args = np.argsort(ranks)
    return args[-1]

class EMTFZoneArtist(object):
  def __init__(self):
    num_cols = (80*64)//emtf_strip_unit
    num_rows = 10
    num_zones = 4
    self.zones = np.empty((num_zones, num_rows, num_cols), dtype=np.object)
    for ind in np.ndindex(self.zones.shape):
      self.zones[ind] = []

  def reset(self):
    for ind in np.ndindex(self.zones.shape):
      del self.zones[ind][:]

  def append(self, hit):
    zones = find_emtf_zones(hit)
    for zone in zones:
      row = find_emtf_zone_row(hit, zone)
      col = find_emtf_zone_col(hit, zone)
      self.zones[zone, row, col].append(hit)

  def _uniquify(self):
    def get_sort_criteria(hit):
      ri = 0 if hit.ring == 1 else 1 # prefer ring 1
      fr = 0 if hit.fr == 1 else 1   # prefer front chamber
      bd = abs(hit.bend)             # prefer smaller bend
      return (ri, fr, bd)

    # Sort and select
    for ind in np.ndindex(self.zones.shape):
      if len(self.zones[ind]) == 0:
        continue

      alist = self.zones[ind][:]  # copy
      if len(alist) > 1:
        alist.sort(key=get_sort_criteria)  # in-place sorting
        self.zones[ind] = alist[0:1]  # only keep the best one

      # Pick up theta ambiguities
      assert(len(self.zones[ind]) == 1)
      ahit = self.zones[ind][0]
      ahit.emtf_theta_alt = 0
      if ahit.type == kCSC:
        for hit in alist:
          if ahit.chamber == hit.chamber and ahit.strip == hit.strip and ahit.wire != hit.wire:
            ahit.emtf_theta_alt = hit.emtf_theta

  def _transform(self, ind):
    if len(self.zones[ind]) == 0:
      return None

    # 7 members: zone, zone_row, emtf_phi, emtf_bend,
    #            emtf_theta, emtf_theta_alt, emtf_qual
    assert(len(self.zones[ind]) == 1)
    hit = self.zones[ind][0]
    arr = np.zeros(7, dtype=np.int32)
    arr[0] = ind[0]
    arr[1] = ind[1]
    arr[2] = hit.emtf_phi
    arr[3] = find_emtf_bend(hit)
    arr[4] = hit.emtf_theta
    arr[5] = hit.emtf_theta_alt
    arr[6] = find_emtf_qual(hit)
    return arr

  def squeeze(self):
    self._uniquify()

    out = []
    for ind in np.ndindex(self.zones.shape):
      arr = self._transform(ind)
      if arr is not None:
        out.append(arr)
    return out

class EMTFZoneScientist(object):
  def __init__(self):
    num_cols = 1
    num_rows = 10
    num_zones = 4
    self.zones = np.empty((num_zones, num_rows, num_cols), dtype=np.object)
    for ind in np.ndindex(self.zones.shape):
      self.zones[ind] = []

  def reset(self):
    for ind in np.ndindex(self.zones.shape):
      del self.zones[ind][:]

  def append(self, hit):
    zones = find_emtf_zones(hit)
    for zone in zones:
      row = find_emtf_zone_row(hit, zone)
      col = 0
      self.zones[zone, row, col].append(hit)

  def _transform(self, ind):
    if len(self.zones[ind]) == 0:
      return None

    # 7 members: zone, zone_row, emtf_phi, emtf_bend,
    #            emtf_theta, emtf_theta_alt, emtf_qual
    assert(len(self.zones[ind]) >= 1)
    sorted_hits = sorted(self.zones[ind], key=lambda hit: hit.layer)
    if sorted_hits[0].type == kRPC or sorted_hits[0].type == kGEM:
      sorted_hits = sorted_hits[0:1]  # only keep one layer
    hits_phi = [hit.emtf_phi for hit in sorted_hits]
    hits_theta = [hit.emtf_theta for hit in sorted_hits]
    arr = np.zeros(7, dtype=np.int32)
    arr[0] = ind[0]
    arr[1] = ind[1]
    arr[2] = pick_the_median(hits_phi)
    arr[3] = (hits_phi[-1] - hits_phi[0])
    arr[4] = pick_the_median(hits_theta)
    arr[5] = 0
    arr[6] = 0
    return arr

  def squeeze(self):
    out = []
    for ind in np.ndindex(self.zones.shape):
      arr = self._transform(ind)
      if arr is not None:
        out.append(arr)
    return out


# ______________________________________________________________________________
# Analyses

class _BaseAnalysis(object):
  """Abstract base class"""
  pass

# ______________________________________________________________________________
class SignalAnalysis(_BaseAnalysis):
  """Prepare signal data used for training.

  Description.
  """

  def run(self, algo):
    out_part = []
    out_hits = []
    out_simhits = []

    sectors = EMTFSectorRanking()
    zones = EMTFZoneArtist()
    zones_simhits = EMTFZoneScientist()

    # __________________________________________________________________________
    # Load tree
    #tree = load_pgun_batch(jobid)
    tree = load_pgun_test()
    verbosity = 0

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Particles
      part = evt.particles[0]  # particle gun

      # First, determine the best sector

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if is_emtf_legit_hit(hit):
          sectors.append(hit)

      best_sector = sectors.get_best_sector()

      # Second, fill the zones with trigger primitives

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if is_emtf_legit_hit(hit) and find_endsec(hit.endcap, hit.sector) == best_sector:
          zones.append(hit)

      # Third, fill the zones with sim hits

      # Sim hits
      for isimhit, simhit in enumerate(evt.simhits):
        simhit.emtf_phi = calc_phi_loc_int(np.rad2deg(simhit.phi), (best_sector%6) + 1)
        simhit.emtf_theta = calc_theta_int(np.rad2deg(simhit.theta), 1 if (part.eta >= 0) else -1)
        if 0 <= simhit.emtf_phi < (80*64):
          zones_simhits.append(simhit)

      # Finally, extract the particle and hits
      ievt_part = np.array([part.pt, part.eta, part.phi, part.invpt, part.d0, part.vz], dtype=np.float32)
      ievt_hits = zones.squeeze()
      ievt_simhits = zones_simhits.squeeze()

      out_part.append(ievt_part)
      out_hits.append(ievt_hits)
      out_simhits.append(ievt_simhits)

      # Debug
      if verbosity >= kINFO:
        print('Processing event: {0}'.format(ievt))
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.type, hit.station, hit.ring, find_endsec(hit.endcap, hit.sector), hit.fr, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.type == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(ihit, hit_id, hit.emtf_phi, hit.emtf_theta, hit.bend, hit.quality, hit_sim_tp, hit.strip, hit.wire))
        for isimhit, simhit in enumerate(evt.simhits):
          simhit_id = (simhit.type, simhit.station, simhit.ring, simhit.layer, simhit.chamber)
          print('.. simhit {0} {1} {2} {3}'.format(isimhit, simhit_id, simhit.phi, simhit.theta))
        print('best sector: {0}'.format(best_sector))
        print('hits:')
        for x in ievt_hits:
          print(x)
        print('simhits:')
        for x in ievt_simhits:
          print(x)

      # Reset
      sectors.reset()
      zones.reset()
      zones_simhits.reset()

    # End loop over events

    # __________________________________________________________________________
    # Output
    outfile = 'signal.npz'
    if use_condor:
      outfile = outfile[:-4] + ('_%i.npz' % jobid)
    out_part = np.asarray(out_part)
    out_hits = create_ragged_array(out_hits)
    out_simhits = create_ragged_array(out_simhits)
    print('[INFO] out_part: {0} out_hits: {1} out_simhits: {2}'.format(out_part.shape, out_hits.shape, out_simhits.shape))
    out_dict = {
      'out_part': out_part,
      'out_hits_values': out_hits.values,
      'out_hits_row_splits': out_hits.row_splits,
      'out_simhits_values': out_simhits.values,
      'out_simhits_row_splits': out_simhits.row_splits,
    }
    save_np_arrays(outfile, out_dict)
    return

# ______________________________________________________________________________
class BkgndAnalysis(_BaseAnalysis):
  """Prepare background data used for training.

  Description.
  """

  def run(self, algo):
    out_hits = []

    num_sectors = 12
    twelve_zones = [EMTFZoneArtist() for i in range(num_sectors)]

    # __________________________________________________________________________
    # Load tree
    tree = load_mixing_batch(jobid)
    verbosity = 0

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # First, apply event veto

      # Particles
      veto = False
      for ipart, part in enumerate(evt.particles):
        if (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4) and (part.pt > 10.):
          veto = True
          break
      if veto:
        continue

      # Second, fill the zones with trigger primitives

      # Trigger primitives
      for sector in range(num_sectors):
        for ihit, hit in enumerate(evt.hits):
          if is_emtf_legit_hit(hit) and find_endsec(hit.endcap, hit.sector) == sector:
            twelve_zones[sector].append(hit)

      # Finally, extract the particle and hits
      for sector in range(num_sectors):
        ievt_hits = twelve_zones[sector].squeeze()
        out_hits.append(ievt_hits)

      # Debug
      if verbosity >= kINFO:
        print('Processing event: {0}'.format(ievt))
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.type, hit.station, hit.ring, find_endsec(hit.endcap, hit.sector), hit.fr, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.type == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(ihit, hit_id, hit.emtf_phi, hit.emtf_theta, hit.bend, hit.quality, hit_sim_tp, hit.strip, hit.wire))
        print('hits:')
        for sector in range(num_sectors):
          ievt_hits = twelve_zones[sector].squeeze()
          for x in ievt_hits:
            print(sector, x)

      # Reset
      for sector in range(num_sectors):
        twelve_zones[sector].reset()

    # End loop over events

    # __________________________________________________________________________
    # Output
    outfile = 'bkgnd.npz'
    if use_condor:
      outfile = outfile[:-4] + ('_%i.npz' % jobid)
    random.shuffle(out_hits)  # shuffle
    out_hits = create_ragged_array(out_hits)
    print('[INFO] out_hits: {0}'.format(out_hits.shape))
    out_dict = {
      'out_hits_values': out_hits.values,
      'out_hits_row_splits': out_hits.row_splits,
    }
    save_np_arrays(outfile, out_dict)
    return

# ______________________________________________________________________________
class EffieAnalysis(_BaseAnalysis):
  """Prepare muon+200PU data used for evaluating efficiency.

  Description.
  """

  def run(self, algo):
    return

# ______________________________________________________________________________
class RatesAnalysis(_BaseAnalysis):
  """Prepare neutrino+200PU data used for evaluating trigger rates.

  Description.
  """

  def run(self, algo):
    return


# ______________________________________________________________________________
# Main

import os, sys, datetime

# Algorithm (pick one)
algo = 'default'  # phase 2
#algo = 'run3'

# Analysis mode (pick one)
analysis = 'signal'
#analysis = 'bkgnd'
#analysis = 'effie'
#analysis = 'rates'

# Job id (pick an integer)
jobid = 0

# Max num of events (-1 means all events)
maxevents = 100

# Condor or not
# if 'CONDOR_EXEC' is defined, overwrite the 3 arguments (algo, analysis, jobid)
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  os.environ['ROOTPY_GRIDMODE'] = 'true'
  algo = sys.argv[1]
  analysis = sys.argv[2]
  jobid = int(sys.argv[3])
  maxevents = -1

# Main function
def main():
  start_time = datetime.datetime.now()
  print('[INFO] Current time    : {0}'.format(start_time))
  print('[INFO] Using cmssw     : {0}'.format(os.environ['CMSSW_VERSION']))
  print('[INFO] Using condor    : {0}'.format(use_condor))
  print('[INFO] Using algo      : {0}'.format(algo))
  print('[INFO] Using analysis  : {0}'.format(analysis))
  print('[INFO] Using jobid     : {0}'.format(jobid))
  print('[INFO] Using maxevents : {0}'.format(maxevents))

  if analysis == 'signal':
    anna = SignalAnalysis()
    anna.run(algo=algo)

  elif analysis == 'bkgnd':
    anna = BkgndAnalysis()
    anna.run(algo=algo)

  elif analysis == 'effie':
    anna = EffieAnalysis()
    anna.run(algo=algo)

  elif analysis == 'rates':
    anna = RatesAnalysis()
    anna.run(algo=algo)

  else:
    raise RuntimeError('Cannot recognize analysis: {0}'.format(analysis))

  # DONE!
  stop_time = datetime.datetime.now()
  print('[INFO] Elapsed time    : {0}'.format(stop_time - start_time))

if __name__ == '__main__':
  main()
