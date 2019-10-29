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

class EMTFChamberCouncil(object):
  def __init__(self):
    # CSC : 18+9+9+9 (native) 3+2+2+2 (neighbor)
    # RPC : 12+12+12+12 (native) 2+2+2+2 (neighbor)
    # iRPC: 0+0+3+3 (native) 0+0+1+1 (neighbor)
    # GEM : 6+3+0+0 (native) 1+1+0+0 (neighbor) 3 (ME0 native) 1 (ME0 neighbor)
    # DT  : 2+2+2+2 (native) 1+1+1+1 (neighbor)
    num_csc_chambers = 54
    num_rpc_chambers = 56
    num_irpc_chambers = 8
    num_gem_chambers = 15
    num_dt_chambers = 12
    num_chambers = (num_csc_chambers, num_rpc_chambers, num_irpc_chambers, num_gem_chambers, num_dt_chambers)
    num_subsystems = len(num_chambers)
    self.chambers = np.empty(num_subsystems, dtype=np.object)
    for i, vi in enumerate(self.chambers):
      self.chambers[i] = np.empty(num_chambers[i], dtype=np.object)
      for j, vj in enumerate(self.chambers[i]):
        self.chambers[i][j] = []

  def reset(self):
    for i, vi in enumerate(self.chambers):
      for j, vj in enumerate(self.chambers[i]):
        del self.chambers[i][j][:]

  def get_index(self, hit):
    pc_station = -1
    pc_chamber = -1
    my_subsystem = -1
    my_chamber = -1
    if hit.type == kCSC:
      if hit.neighbor == 0:
        # ME1 subsector: 1-2 cscid: 1-9 -> pc_station: 0-1 pc_chamber: 0-8
        if hit.station == 1:
          pc_station = hit.subsector-1
          pc_chamber = hit.cscid-1
        # ME2,3,4 station: 2-4 cscid: 1-9 -> pc_station: 2-4 pc_chamber: 0-8
        else:
          pc_station = hit.station
          pc_chamber = hit.cscid-1
      else:
        # ME1n cscid: 3,6,9 -> pc_station: 5 pc_chamber: 0-2
        if hit.station == 1:
          pc_station = 5
          pc_chamber = (hit.cscid-1)//3
        # ME2n,3n,4n cscid: 3,9 -> pc_station: 5 pc_chamber: 3-8
        else:
          pc_station = 5
          pc_chamber = (hit.station) * 2 - 1 + (0 if (hit.cscid == 3) else 1)
      my_subsystem = 0
      my_chamber = (pc_station * 9) + pc_chamber
    elif hit.type == kRPC and not ((hit.station == 3 or hit.station == 4) and hit.ring == 1):
      if hit.neighbor == 0:
        # RE1/2 subsector: 1-2 cscid: 4-6 -> pc_station: 0 pc_chamber: 0-5
        # RE1/3 subsector: 1-2 cscid: 7-9 -> pc_station: 1 pc_chamber: 0-5
        if hit.station == 1:
          pc_station = (hit.station-1) * 2 + (hit.ring-2)
          pc_chamber = (hit.subsector-1) * 3 + (hit.cscid-4) - (hit.ring-2) * 3
        # RE2,3,4 station: 2-4 ring: 2-3 cscid: 4-9 -> pc_station: 2-7 pc_chamber: 0-5
        else:
          pc_station = (hit.station-1) * 2 + (hit.ring-2)
          pc_chamber = hit.cscid-4
      else:
        # RE1n,2n,3n,4n cscid: 6,9 -> pc_station: 0-7 pc_chamber: 6
        pc_station = (hit.station-1) * 2 + (hit.ring-2)
        pc_chamber = 7-1
      my_subsystem = 1
      my_chamber = (pc_station * 7) + pc_chamber
    elif hit.type == kRPC:  # iRPC
      if hit.neighbor == 0:
        # RE3/1,4/1 station: 3-4 cscid: 1-3 -> pc_station: 0-1 pc_chamber: 0-2
        pc_station = (hit.station-3)
        pc_chamber = hit.cscid-1
      else:
        # RE3/1n,4/1n cscid: 3 -> pc_station: 0-1 pc_chamber: 3
        pc_station = (hit.station-3)
        pc_chamber = 4-1
      my_subsystem = 2
      my_chamber = (pc_station * 4) + pc_chamber
    elif hit.type == kGEM or hit.type == kME0:
      if hit.neighbor == 0:
        # GE1/1 subsector: 1-2 cscid: 1-3 -> pc_station: 0-1 pc_chamber: 0-2
        if hit.type == kGEM and hit.station == 1:
          pc_station = hit.subsector-1
          pc_chamber = hit.cscid-1
        # GE2/1 station: 2 cscid: 1-3 -> pc_station: 2 pc_chamber: 0-2
        elif hit.type == kGEM and hit.station == 2:
          pc_station = 2
          pc_chamber = hit.cscid-1
        # ME0 station: 1 cscid: 1-3 -> pc_station: 3 pc_chamber: 0-2
        elif hit.type == kME0:
          pc_station = 3
          pc_chamber = hit.cscid-1
      else:
        # GE1/1n cscid: 3 -> pc_station: 4 pc_chamber: 0
        if hit.type == kGEM and hit.station == 1:
          pc_station = 4
          pc_chamber = 0
        # GE2/1n cscid: 3 -> pc_station: 4 pc_chamber: 1
        elif hit.type == kGEM and hit.station == 2:
          pc_station = 4
          pc_chamber = 1
        # ME0n cscid: 3 -> pc_station: 4 pc_chamber: 2
        elif hit.type == kME0:
          pc_station = 4
          pc_chamber = 2
      my_subsystem = 3
      my_chamber = (pc_station * 3) + pc_chamber
    elif hit.type == kDT:
      if hit.neighbor == 0:
        # MB1,2,3,4 station: 1-4 cscid: 6,9 -> pc_station: 0-3 pc_chamber: 0-1
        pc_station = (hit.station-1)
        pc_chamber = (0 if (hit.cscid == 6) else 1)
      else:
        # MB1n,2n,3n,4n cscid: 9 -> pc_station: 0-3 pc_chamber: 2
        pc_station = (hit.station-1)
        pc_chamber = 3-1
      my_subsystem = 4
      my_chamber = (pc_station * 3) + pc_chamber

    assert(pc_station != -1)
    assert(pc_chamber != -1)
    assert(my_subsystem != -1)
    assert(my_chamber != -1)
    return (my_subsystem, my_chamber)

  def append(self, hit):
    my_subsystem, my_chamber = self.get_index(hit)
    self.chambers[my_subsystem][my_chamber].append(hit)

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
      ahit.emtf_theta_alt = ahit.emtf_theta
      if ahit.type == kCSC:
        for hit in alist:
          if ahit.chamber == hit.chamber and ahit.strip == hit.strip and ahit.wire != hit.wire:
            ahit.emtf_theta_alt = hit.emtf_theta

  def _transform(self, ind):
    if len(self.zones[ind]) == 0:
      return None

    # 8 members: zone, zone_row, emtf_phi, emtf_bend,
    #            emtf_theta, emtf_theta_alt, emtf_qual, emtf_time
    assert(len(self.zones[ind]) == 1)
    hit = self.zones[ind][0]
    arr = np.zeros(8, dtype=np.int32)
    arr[0] = ind[0]
    arr[1] = ind[1]
    arr[2] = hit.emtf_phi
    arr[3] = find_emtf_bend(hit)
    arr[4] = hit.emtf_theta
    arr[5] = hit.emtf_theta_alt
    arr[6] = find_emtf_qual(hit)
    arr[7] = find_emtf_time(hit)
    return arr

  def squeeze(self):
    self._uniquify()

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
    sectors = EMTFSectorRanking()
    chambers = EMTFChamberCouncil()
    zones = EMTFZoneArtist()

    # __________________________________________________________________________
    # Load tree
    tree = load_pgun_batch(jobid)
    verbosity = 1

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

      # End loop over trigger primitives

      best_sector = sectors.get_best_sector()
      sectors.reset()

      # Second, fill the chambers and the zones
      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if is_emtf_legit_hit(hit) and find_endsec(hit.endcap, hit.sector) == best_sector:
          chambers.append(hit)
          zones.append(hit)

      # End loop over trigger primitives

      # Finally, extract the hits
      zone_hits = zones.squeeze()


      if verbosity >= kINFO:
        print('Processing event: {0}'.format(ievt))
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.type, hit.station, hit.ring, find_endsec(hit.endcap, hit.sector), hit.fr, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.type == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(ihit, hit_id, hit.emtf_phi, hit.emtf_theta, hit.bend, hit.quality, hit_sim_tp, hit.strip, hit.wire))

        print('best sector: {0}'.format(best_sector))
        #print('chambers: {0}'.format(chambers.chambers[0]))
        #print('zones: {0}'.format(zones.zones[0]))
        print('zone hits:')
        for x in zone_hits:
          print(x)

      chambers.reset()
      zones.reset()

    # End loop over events

    # __________________________________________________________________________

    return

# ______________________________________________________________________________
class NoiseAnalysis(_BaseAnalysis):
  """Prepare noise data used for training.

  Description.
  """

  def run(self, algo):
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
#analysis = 'noise'
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

  elif analysis == 'noise':
    anna = NoiseAnalysis()
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
