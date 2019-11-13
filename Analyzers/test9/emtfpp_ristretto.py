"""Data preparation."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtfpp_algos import *


# ______________________________________________________________________________
# Classes

class EMTFSectorRanking(object):
  def __init__(self):
    self.sectors = np.zeros(num_emtf_sectors, dtype=np.int32)

  def reset(self):
    self.sectors *= 0

  def add(self, hit):
    endsec = find_endsec(hit.endcap, hit.sector)
    hit_lay = find_emtf_layer(hit)
    a, b, c = hit_lay//100, (hit_lay//10)%10, hit_lay%10  # type, station, ring

    rank = np.int32(0)
    if a == 6:              # a = (6,) b = (0,) -> bit position (11,)
      rank |= (1 << (11 + 0))
    elif a == 5:            # a = (5,) b = (0,) -> bit position (10,)
      rank |= (1 << (10 + 0))
    elif a < 5 and b == 0:  # a = (0,1,2,3,4,) b = (0,) -> bit position (9,8,7,6,5,)
      rank |= (1 << (5 + (4 - a)))
    elif a < 5 and b == 1:  # a = (0,1,2,3,4,) b = (1,) -> bit position (4,3,2,1,0,)
      rank |= (1 << (0 + (4 - a)))
    self.sectors[endsec] |= rank

  def get_best_sector(self):
    argsorted = np.argsort(self.sectors)
    best_sector = argsorted[-1]
    best_sector_rank = self.sectors[best_sector]

    cnt0 = np.sum([bool(best_sector_rank & (1 << i)) for i in range(12)])  # count any station
    cnt1 = np.sum([bool(best_sector_rank & (1 << i)) for i in range(8,12)])  # count station 1
    if not (cnt0 >= 2 and cnt1 >= 1):
      best_sector = None
    return best_sector

class EMTFZoneArtist(object):
  def __init__(self):
    num_cols = max_emtf_strip//coarse_emtf_strip
    num_rows = 10
    num_zones = 4
    self.zones = np.empty((num_zones, num_rows, num_cols), dtype=np.object)
    for ind in np.ndindex(self.zones.shape):
      self.zones[ind] = []

  def reset(self):
    for ind in np.ndindex(self.zones.shape):
      if len(self.zones[ind]) == 0:
        continue
      del self.zones[ind][:]

  def add(self, hit):
    zones = find_emtf_zones(hit)
    for zone in zones:
      row = find_emtf_zone_row(hit, zone)
      col = find_emtf_zone_col(hit)
      self.zones[zone, row, col].append(hit)

  def _uniquify(self):
    def get_sort_criteria(hit):
      ri = 0 if hit.ring == 1 else 1 # prefer ring 1
      fr = 0 if hit.fr == 1 else 1   # prefer front chamber
      bd = np.abs(hit.bend)          # prefer smaller bend
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

    # 8 members: zone, zone_row, zone_col, emtf_phi,
    #            emtf_bend, emtf_theta, emtf_theta_alt, emtf_qual
    assert(len(self.zones[ind]) == 1)
    hit = self.zones[ind][0]

    arr = np.zeros(8, dtype=np.int32)
    arr[0] = ind[0]
    arr[1] = ind[1]
    arr[2] = ind[2]
    arr[3] = hit.emtf_phi
    arr[4] = find_emtf_bend(hit)
    arr[5] = hit.emtf_theta
    arr[6] = hit.emtf_theta_alt
    arr[7] = find_emtf_qual(hit)
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
      if len(self.zones[ind]) == 0:
        continue
      del self.zones[ind][:]

  def add(self, hit):
    zones = find_emtf_zones(hit)
    for zone in zones:
      row = find_emtf_zone_row(hit, zone)
      col = 0
      self.zones[zone, row, col].append(hit)

  def _transform(self, ind):
    if len(self.zones[ind]) == 0:
      return None

    # 8 members: zone, zone_row, zone_col, emtf_phi,
    #            emtf_bend, emtf_theta, emtf_theta_alt, emtf_qual
    assert(len(self.zones[ind]) >= 1)
    sorted_hits = sorted(self.zones[ind], key=lambda hit: hit.layer)
    if sorted_hits[0].type == kRPC or sorted_hits[0].type == kGEM:
      sorted_hits = sorted_hits[0:1]  # only keep one layer
    hit = pick_the_median(sorted_hits)

    arr = np.zeros(8, dtype=np.int32)
    arr[0] = ind[0]
    arr[1] = ind[1]
    arr[2] = find_emtf_zone_col(hit)
    arr[3] = hit.emtf_phi
    arr[4] = 0
    arr[5] = hit.emtf_theta
    arr[6] = 0
    arr[7] = 0
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
    tree = load_pgun_displ_batch(jobid)
    #tree = load_pgun_test()

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      if (ievt % 1000) == 0:
        print('Processing event {0}'.format(ievt))

      # Reset
      sectors.reset()
      zones.reset()
      zones_simhits.reset()

      # Particles
      try:
        part = evt.particles[0]  # particle gun
      except:
        continue

      # First, determine the best sector

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if is_emtf_legit_hit(hit):
          sectors.add(hit)

      best_sector = sectors.get_best_sector()
      if best_sector is None:
        continue

      # Second, fill the zones with trigger primitives

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if is_emtf_legit_hit(hit) and find_endsec(hit.endcap, hit.sector) == best_sector:
          zones.add(hit)

      # Third, fill the zones with sim hits

      # Sim hits
      for isimhit, simhit in enumerate(evt.simhits):
        simhit.emtf_phi = calc_phi_loc_int(np.rad2deg(simhit.phi), (best_sector%6) + 1)
        simhit.emtf_theta = calc_theta_int(np.rad2deg(simhit.theta), 1 if ((best_sector//6) == 0) else -1)
        if 0 <= simhit.emtf_phi < max_emtf_strip:
          zones_simhits.add(simhit)

      # Finally, extract the particle and hits
      ievt_part = np.array([part.pt, part.eta, part.phi, part.invpt, part.d0, part.vx, part.vy, part.vz], dtype=np.float32)
      ievt_hits = zones.squeeze()
      ievt_simhits = zones_simhits.squeeze()

      out_part.append(ievt_part)
      out_hits.append(ievt_hits)
      out_simhits.append(ievt_simhits)

      # Debug
      if verbosity >= kINFO:
        print('Processing event {0}'.format(ievt))
        print('.. part {0} {1} {2} {3} {4} {5}'.format(0, part.pt, part.eta, part.phi, part.invpt, part.d0))
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

    sector_zones = [EMTFZoneArtist() for i in range(num_emtf_sectors)]

    # __________________________________________________________________________
    # Load tree
    tree = load_mixing_batch(jobid)

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      if (ievt % 100) == 0:
        print('Processing event {0}'.format(ievt))

      # Reset
      for zones in sector_zones:
        zones.reset()

      # First, apply event veto

      # Particles
      veto = False
      for ipart, part in enumerate(evt.particles):
        if (part.bx == 0) and (part.pt > 10.) and (1.2 <= np.abs(calc_etastar_from_eta(part.eta, part.vx, part.vy, part.vz)) <= 2.5):
          veto = True
          break
      if veto:
        continue

      # Second, fill the zones with trigger primitives

      # Trigger primitives
      for sector in range(num_emtf_sectors):
        for ihit, hit in enumerate(evt.hits):
          if is_emtf_legit_hit(hit) and find_endsec(hit.endcap, hit.sector) == sector:
            sector_zones[sector].add(hit)

      # Finally, extract the particle and hits
      for zones in sector_zones:
        ievt_hits = zones.squeeze()
        out_hits.append(ievt_hits)

      # Debug
      if verbosity >= kINFO:
        print('Processing event {0}'.format(ievt))
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.type, hit.station, hit.ring, find_endsec(hit.endcap, hit.sector), hit.fr, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.type == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(ihit, hit_id, hit.emtf_phi, hit.emtf_theta, hit.bend, hit.quality, hit_sim_tp, hit.strip, hit.wire))
        print('hits:')
        for sector in range(num_emtf_sectors):
          ievt_hits = out_hits[-(12-sector)]
          for x in ievt_hits:
            print(sector, x)

    # End loop over events

    # __________________________________________________________________________
    # Output
    outfile = 'bkgnd.npz'
    if use_condor:
      outfile = outfile[:-4] + ('_%i.npz' % jobid)
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

# Verbosity
verbosity = 1

# Condor or not
# if 'CONDOR_EXEC' is defined, overwrite the 3 arguments (algo, analysis, jobid)
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  os.environ['ROOTPY_GRIDMODE'] = 'true'
  algo = sys.argv[1]
  analysis = sys.argv[2]
  jobid = int(sys.argv[3])
  maxevents = -1
  verbosity = 0

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
