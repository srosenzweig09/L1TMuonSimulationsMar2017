"""Algorithms in EMTF++."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtf_utils import *


# ______________________________________________________________________________
# Classes

num_emtf_sectors = 12

min_emtf_strip = 5*64    # 5 deg

max_emtf_strip = 80*64   # 80 deg

coarse_emtf_strip = 8*2  # 'doublestrip' unit

emtf_eta_bins = (0.8, 1.2, 1.55, 1.98, 2.5)

# Decide EMTF hit layer number
class EMTFLayer(object):
  def __init__(self):
    lut = np.zeros((5,5,5), dtype=np.int32) - 99  # (type, station, ring) -> layer
    lut[1,1,4] = (0) * 100 + (0) * 10 + (1) # ME1/1a
    lut[1,1,1] = (0) * 100 + (0) * 10 + (1) # ME1/1b
    lut[1,1,2] = (1) * 100 + (0) * 10 + (2) # ME1/2
    lut[1,1,3] = (1) * 100 + (0) * 10 + (3) # ME1/3
    lut[1,2,1] = (2) * 100 + (0) * 10 + (1) # ME2/1
    lut[1,2,2] = (2) * 100 + (0) * 10 + (2) # ME2/2
    lut[1,3,1] = (3) * 100 + (0) * 10 + (1) # ME3/1
    lut[1,3,2] = (3) * 100 + (0) * 10 + (2) # ME3/2
    lut[1,4,1] = (4) * 100 + (0) * 10 + (1) # ME4/1
    lut[1,4,2] = (4) * 100 + (0) * 10 + (2) # ME4/2
    lut[2,1,2] = (1) * 100 + (1) * 10 + (2) # RE1/2
    lut[2,1,3] = (1) * 100 + (1) * 10 + (3) # RE1/3
    lut[2,2,2] = (2) * 100 + (1) * 10 + (2) # RE2/2
    lut[2,2,3] = (2) * 100 + (1) * 10 + (3) # RE2/3
    lut[2,3,1] = (3) * 100 + (1) * 10 + (1) # RE3/1
    lut[2,3,2] = (3) * 100 + (1) * 10 + (2) # RE3/2
    lut[2,3,3] = (3) * 100 + (1) * 10 + (3) # RE3/3
    lut[2,4,1] = (4) * 100 + (1) * 10 + (1) # RE4/1
    lut[2,4,2] = (4) * 100 + (1) * 10 + (2) # RE4/2
    lut[2,4,3] = (4) * 100 + (1) * 10 + (3) # RE4/3
    lut[3,1,1] = (0) * 100 + (1) * 10 + (1) # GE1/1
    lut[3,2,1] = (2) * 100 + (1) * 10 + (1) # GE2/1
    lut[4,1,1] = (5) * 100 + (0) * 10 + (1) # ME0
    lut[0,1,1] = (6) * 100 + (0) * 10 + (1) # MB1
    lut[0,2,1] = (6) * 100 + (0) * 10 + (2) # MB2
    lut[0,3,1] = (6) * 100 + (0) * 10 + (3) # MB3
    lut[0,4,1] = (6) * 100 + (0) * 10 + (4) # MB4
    self.lut = lut

  def __call__(self, hit):
    index = (hit.type, hit.station, hit.ring)
    emtf_layer = self.lut[index]
    return emtf_layer

# Decide EMTF hit zones
class EMTFZone(object):
  def __init__(self):
    lut = np.zeros((5,5,5,4,2), dtype=np.int32) - 99  # ((type, station, ring), zone) -> (min_theta,max_theta)
    lut[1,1,4][0] = 4,26    # ME1/1a
    lut[1,1,4][1] = 24,53   # ME1/1a
    lut[1,1,1][0] = 4,26    # ME1/1b
    lut[1,1,1][1] = 24,53   # ME1/1b
    lut[1,1,2][1] = 46,54   # ME1/2
    lut[1,1,2][2] = 52,88   # ME1/2
    lut[1,1,3][3] = 98,125  # ME1/3
    #
    lut[1,2,1][0] = 4,25    # ME2/1
    lut[1,2,1][1] = 23,49   # ME2/1
    lut[1,2,2][2] = 53,88   # ME2/2
    lut[1,2,2][3] = 80,111  # ME2/2
    #
    lut[1,3,1][0] = 4,25    # ME3/1
    lut[1,3,1][1] = 23,40   # ME3/1
    lut[1,3,2][1] = 44,54   # ME3/2
    lut[1,3,2][2] = 51,88   # ME3/2
    lut[1,3,2][3] = 80,98   # ME3/2
    #
    lut[1,4,1][0] = 4,25    # ME4/1
    lut[1,4,1][1] = 23,35   # ME4/1
    lut[1,4,2][1] = 38,54   # ME4/2
    lut[1,4,2][2] = 51,88   # ME4/2
    lut[1,4,2][3] = 80,91   # ME4/2
    #
    lut[2,1,2][1] = 52,56   # RE1/2
    lut[2,1,2][2] = 52,84   # RE1/2
    lut[2,1,3][3] = 100,120 # RE1/3
    lut[2,2,2][1] = 56,56   # RE2/2
    lut[2,2,2][2] = 56,76   # RE2/2
    lut[2,2,3][2] = 88,88   # RE2/3
    lut[2,2,3][3] = 88,112  # RE2/3
    #
    lut[2,3,1][0] = 4,25    # RE3/1
    lut[2,3,1][1] = 23,36   # RE3/1
    lut[2,3,2][1] = 40,52   # RE3/2
    lut[2,3,2][2] = 48,60   # RE3/2
    lut[2,3,3][2] = 72,84   # RE3/3
    lut[2,3,3][3] = 80,92   # RE3/3
    lut[2,4,1][0] = 4,25    # RE4/1
    lut[2,4,1][1] = 23,31   # RE4/1
    lut[2,4,2][1] = 36,52   # RE4/2
    lut[2,4,2][2] = 52,52   # RE4/2
    lut[2,4,3][2] = 64,84   # RE4/3
    #
    lut[3,1,1][0] = 17,26   # GE1/1
    lut[3,1,1][1] = 24,52   # GE1/1
    lut[3,2,1][0] = 4,25    # GE2/1
    lut[3,2,1][1] = 23,46   # GE2/1
    #
    lut[4,1,1][0] = 4,23    # ME0
    #
    lut[0,1,1][3] = 92,130  # MB1
    lut[0,2,1][3] = 108,138 # MB2
    lut[0,3,1][3] = 126,144 # MB3
    self.lut = lut

  def __call__(self, hit):
    index = (hit.type, hit.station, hit.ring)
    entry = self.lut[index]
    answer = (entry[:,0] <= hit.emtf_theta) & (hit.emtf_theta <= entry[:,1])
    zones = np.nonzero(answer)
    if isinstance(zones, tuple):
      zones = zones[0]
    return zones

# Decide EMTF hit row number in a zone
class EMTFZoneRow(object):
  def __init__(self):
    lut = np.zeros((5,5,5,4), dtype=np.int32) - 99  # ((type, station, ring), zone) -> row number
    lut[1,1,4][0] = 2 # ME1/1a
    lut[1,1,4][1] = 1 # ME1/1a
    lut[1,1,1][0] = 2 # ME1/1b
    lut[1,1,1][1] = 1 # ME1/1b
    lut[1,1,2][1] = 2 # ME1/2
    lut[1,1,2][2] = 0 # ME1/2
    lut[1,1,3][3] = 3 # ME1/3
    #
    lut[1,2,1][0] = 4 # ME2/1
    lut[1,2,1][1] = 5 # ME2/1
    lut[1,2,2][2] = 3 # ME2/2
    lut[1,2,2][3] = 6 # ME2/2
    #
    lut[1,3,1][0] = 5 # ME3/1
    lut[1,3,1][1] = 6 # ME3/1
    lut[1,3,2][1] = 6 # ME3/2
    lut[1,3,2][2] = 4 # ME3/2
    lut[1,3,2][3] = 7 # ME3/2
    #
    lut[1,4,1][0] = 7 # ME4/1
    lut[1,4,1][1] = 8 # ME4/1
    lut[1,4,2][1] = 8 # ME4/2
    lut[1,4,2][2] = 6 # ME4/2
    lut[1,4,2][3] = 9 # ME4/2
    #
    lut[2,1,2][1] = 3 # RE1/2
    lut[2,1,2][2] = 1 # RE1/2
    lut[2,1,3][3] = 4 # RE1/3
    lut[2,2,2][1] = 4 # RE2/2
    lut[2,2,2][2] = 2 # RE2/2
    lut[2,2,3][2] = 2 # RE2/3
    lut[2,2,3][3] = 5 # RE2/3
    #
    lut[2,3,1][0] = 6 # RE3/1
    lut[2,3,1][1] = 7 # RE3/1
    lut[2,3,2][1] = 7 # RE3/2
    lut[2,3,2][2] = 5 # RE3/2
    lut[2,3,3][2] = 5 # RE3/3
    lut[2,3,3][3] = 8 # RE3/3
    lut[2,4,1][0] = 8 # RE4/1
    lut[2,4,1][1] = 9 # RE4/1
    lut[2,4,2][1] = 9 # RE4/2
    lut[2,4,2][2] = 7 # RE4/2
    lut[2,4,3][2] = 7 # RE4/3
    #
    lut[3,1,1][0] = 1 # GE1/1
    lut[3,1,1][1] = 0 # GE1/1
    lut[3,2,1][0] = 3 # GE2/1
    lut[3,2,1][1] = 4 # GE2/1
    #
    lut[4,1,1][0] = 0 # ME0
    #
    lut[0,1,1][3] = 0 # MB1
    lut[0,2,1][3] = 1 # MB2
    lut[0,3,1][3] = 2 # MB3
    self.lut = lut

  def __call__(self, hit, zone):
    index = (hit.type, hit.station, hit.ring, zone)
    zone_row = self.lut[index]
    return zone_row

# Decide EMTF hit column number in a zone
class EMTFZoneCol(object):
  def __call__(self, hit):
    zone_col = np.int32(hit.emtf_phi)
    zone_col -= min_emtf_strip
    zone_col //= coarse_emtf_strip
    return zone_col

# Decide EMTF hit phi (integer unit)
class EMTFPhi(object):
  def __call__(self, hit):
    emtf_phi = np.int32(hit.emtf_phi)
    return emtf_phi

# Decide EMTF hit phi (integer unit) (old version)
class EMTFOldPhi(object):
  def __init__(self):
    self.ph_pattern_corr_lut = np.array([0, 0, 5, 5, 5, 5, 2, 2, 2, 2, 0], dtype=np.int32)
    self.ph_pattern_corr_sign_lut = np.array([0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0], dtype=np.int32)
    self.ph_init_lut = np.array([
      1324,1926,2527,1332,1932,2532,1415,2015,2615,1332,1932,2532,726,732,815,732,3128,3725,4325,3132,3732,4332,3214,3815,4415,3131,3732,4331,1316,2516,
      3716,1332,1932,2532,3132,3732,4332,116,732,2580,3781,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,
      1364,1326,1923,2527,1332,1932,2532,1415,2015,2614,1331,1931,2531,725,732,815,731,3123,3724,4326,3132,3732,4332,3214,3814,4415,3131,3731,4331,1316,
      2516,3715,1332,1932,2532,3132,3731,4331,116,732,2580,3780,4980,1964,2564,3164,3763,4363,4964,1380,1364,2580,3780,4979,1964,2564,3164,3763,4363,4963,
      1380,1363,1323,1926,2525,1331,1932,2532,1415,2015,2615,1331,1931,2531,726,732,816,732,3124,3727,4325,3132,3732,4332,3215,3815,4415,3131,3731,4332,
      1316,2516,3716,1332,1932,2532,3131,3731,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1381,1364,2580,3780,4980,1964,2564,3164,3764,4364,
      4964,1380,1364,1324,1929,2531,1332,1932,2532,1416,2015,2615,1331,1932,2531,725,732,815,732,3123,3728,4327,3132,3732,4332,3215,3815,4416,3132,3731,
      4332,1316,2516,3716,1332,1932,2532,3132,3733,4332,116,732,2580,3781,4980,1964,2564,3165,3765,4365,4964,1380,1364,2580,3781,4981,1964,2564,3164,3765,
      4365,4965,1380,1364,1325,1925,2524,1332,1932,2532,1415,2015,2615,1331,1931,2531,727,732,815,731,3124,3726,4325,3132,3732,4332,3215,3815,4415,3132,
      3731,4331,1316,2516,3716,1332,1932,2532,3132,3732,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,3164,
      3764,4364,4964,1380,1364,1321,1927,2524,1332,1932,2532,1415,2015,2615,1331,1931,2532,725,732,815,731,3128,3727,4326,3133,3732,4332,3215,3815,4415,
      3131,3731,4332,1316,2516,3716,1332,1932,2532,3132,3732,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,
      3164,3764,4364,4964,1380,1364,1979,2578,3178,1972,2572,3172,1890,2489,3090,1973,2573,3173,1380,1372,1289,1375,3779,4380,4978,3772,4372,4972,3689,4289,
      4889,3772,4373,4972,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,3724,1340,
      1940,2540,3140,3740,4340,124,740,1979,2578,3179,1972,2572,3172,1889,2489,3089,1973,2573,3173,1378,1372,1289,1372,3778,4380,4982,3772,4372,4972,3689,
      4289,4890,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,3724,
      1340,1940,2540,3140,3740,4340,124,740,1977,2580,3179,1972,2572,3172,1889,2489,3089,1975,2572,3173,1382,1372,1289,1372,3779,4379,4979,3772,4372,4972,
      3688,4289,4889,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,
      3724,1340,1940,2540,3140,3740,4340,124,740,1979,2577,3180,1972,2572,3172,1889,2489,3089,1973,2573,3173,1379,1372,1289,1373,3780,4378,4979,3772,4372,
      4972,3689,4289,4889,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,
      2524,3724,1340,1940,2541,3141,3741,4341,124,740,1978,2580,3179,1972,2572,3172,1889,2489,3089,1973,2573,3173,1380,1372,1290,1373,3780,4378,4981,3772,
      4372,4972,3689,4290,4889,3775,4372,4977,2588,3787,4987,1972,2572,3172,3772,4372,4972,1388,1372,1324,2523,3723,1340,1940,2540,3140,3740,4340,124,740,
      1324,2523,3724,1341,1941,2541,3141,3741,4341,124,741,1979,2581,3178,1973,2573,3173,1890,2490,3090,1973,2575,3173,1382,1373,1290,1377,3779,4380,4981,
      3773,4373,4973,3690,4290,4890,3774,4374,4976,2589,3789,4989,1973,2573,3173,3773,4373,4973,1388,1373,1325,2525,3725,1341,1941,2541,3141,3741,4341,124,
      741,1325,2525,3725,1341,1941,2541,3141,3741,4341,124,742
    ], dtype=np.int32)
    assert(len(self.ph_init_lut) == 2*6*61)

  def __call__(self, hit):
    emtf_phi = np.int32(hit.emtf_phi)
    clct_pattern = np.int32(hit.pattern)
    if hit.type == kCSC:
      # Is this chamber mounted in reverse direction?
      # (i.e., phi vs. strip number is reversed)
      ph_reverse = (hit.endcap == 1 and hit.station >= 3) or (hit.endcap == -1 and hit.station < 3)

      # Is this 10-deg or 20-deg chamber?
      is_10degree = (hit.station == 1) or (hit.station >= 2 and hit.ring == 2)  # ME1 and ME2,3,4/2

      pc_station = -1
      pc_chamber = -1
      if hit.neighbor == 0:
        if hit.station == 1:  # ME1: 0 - 8, 9 - 17
          pc_station = hit.subsector-1
          pc_chamber = hit.cscid-1
        else:                 # ME2,3,4: 18 - 26, 27 - 35, 36 - 44
          pc_station = hit.station
          pc_chamber = hit.cscid-1
      else:
        if hit.station == 1:  # ME1n: 45 - 47
          pc_station = 5
          pc_chamber = (hit.cscid-1)//3
        else:                 # ME2n,3n,4n: 48 - 53
          pc_station = 5
          pc_chamber = (hit.station) * 2 - 1 + (0 if (hit.cscid-1 < 3) else 1)
      assert(pc_station != -1)
      assert(pc_chamber != -1)

      pc_lut_id = pc_chamber
      if pc_station == 0:    # ME1 sub 1: 0 - 11
        pc_lut_id = pc_lut_id + 9 if (hit.ring == 4) else pc_lut_id
      elif pc_station == 1:  # ME1 sub 2: 16 - 27
        pc_lut_id += 16
        pc_lut_id = pc_lut_id + 9 if (hit.ring == 4) else pc_lut_id
      elif pc_station == 2:  # ME2: 28 - 36
        pc_lut_id += 28
      elif pc_station == 3:  # ME3: 39 - 47
        pc_lut_id += 39
      elif pc_station == 4:  # ME4 : 50 - 58
        pc_lut_id += 50
      elif pc_station == 5 and pc_chamber < 3:  # neighbor ME1: 12 - 15
        pc_lut_id = pc_lut_id + 15 if (hit.ring == 4) else pc_lut_id + 12
      elif pc_station == 5 and pc_chamber < 5:  # neighbor ME2: 37 - 38
        pc_lut_id += 28 + 9 - 3
      elif pc_station == 5 and pc_chamber < 7:  # neighbor ME3: 48 - 49
        pc_lut_id += 39 + 9 - 5
      elif pc_station == 5 and pc_chamber < 9:  # neighbor ME4: 59 - 60
        pc_lut_id += 50 + 9 - 7
      assert(pc_lut_id < 61)

      fw_sector = hit.sector-1
      fw_endcap = 0 if (hit.endcap == 1) else 1
      fw_strip = hit.strip  # already starts from 0

      # Apply phi correction from CLCT pattern number
      clct_pat_corr = self.ph_pattern_corr_lut[clct_pattern]
      clct_pat_corr_sign_ = self.ph_pattern_corr_sign_lut[clct_pattern]
      clct_pat_corr_sign = 1 if (clct_pat_corr_sign_ == 0) else -1
      if fw_strip == 0 and clct_pat_corr_sign == -1:
        clct_pat_corr = 0

      if is_10degree:
        eighth_strip = fw_strip << 2  # full precision, uses only 2 bits of pattern correction
        eighth_strip += clct_pat_corr_sign * (clct_pat_corr >> 1)
      else:
        eighth_strip = fw_strip << 3  # multiply by 2, uses all 3 bits of pattern correction
        eighth_strip += clct_pat_corr_sign * (clct_pat_corr >> 0)
      assert(eighth_strip >= 0)

      # Multiplicative factor for eighth_strip
      factor = 1024
      if hit.station == 1 and hit.ring == 4:
        factor = 1707
      elif hit.station == 1 and hit.ring == 1:
        factor = 1301
      elif hit.station == 1 and hit.ring == 3:
        factor = 947

      ph_tmp = (eighth_strip * factor) >> 10
      ph_tmp_sign = 1 if (ph_reverse == 0) else -1

      endsec_pc_lut_id = (fw_endcap * 6 + fw_sector) * 61 + pc_lut_id
      fph = self.ph_init_lut[endsec_pc_lut_id] + ph_tmp_sign * ph_tmp
      assert(0 <= fph and fph < 5000)
      #
      emtf_phi = fph
    else:  # non-CSC
      pass
    return emtf_phi

# Decide EMTF hit theta (integer unit)
class EMTFTheta(object):
  def __call__(self, hit):
    assert(hit.emtf_theta > 0)
    emtf_theta = np.int32(hit.emtf_theta)
    if hit.type == kDT:
      # wire -1 means no theta SL
      # quality 0&1 are RPC digis
      if (hit.wire == -1) or (hit.quality < 2):
        if hit.station == 1:
          emtf_theta = 112
        elif hit.station == 2:
          emtf_theta = 122
        elif hit.station == 3:
          emtf_theta = 131
      else:
        pass
    else:
      pass
    return emtf_theta

# Decide EMTF hit bend
class EMTFBend(object):
  def __call__(self, hit):
    emtf_bend = np.int32(hit.bend)
    if hit.type == kCSC:
      # Special case for ME1/1a:
      # rescale the bend to the same scale as ME1/1b
      if hit.station == 1 and hit.ring == 4:
        emtf_bend = np.round(emtf_bend.astype(np.float32) * 0.026331/0.014264).astype(np.int32)
        emtf_bend = np.clip(emtf_bend, -32, 31)
      emtf_bend *= hit.endcap
      emtf_bend = np.round(emtf_bend.astype(np.float32) * 0.5001).astype(np.int32)  # from 1/32-strip unit to 1/16-strip unit
      emtf_bend = np.clip(emtf_bend, -16, 15)
    elif hit.type == kME0:
      emtf_bend = np.round(emtf_bend.astype(np.float32) * 0.5001).astype(np.int32)  # from 1/4-strip unit to 1/2-strip unit
      emtf_bend = np.clip(emtf_bend, -64, 63)
    elif hit.type == kDT:
      if hit.quality >= 4:
        emtf_bend = np.clip(emtf_bend, -512, 511)
      else:
        #emtf_bend = np.int32(0)
        emtf_bend = np.clip(emtf_bend, -512, 511)
    else:  # kRPC, kGEM
      emtf_bend = np.int32(0)
    return emtf_bend

# Decide EMTF hit bend (old version)
class EMTFOldBend(object):
  def __init__(self):
    self.lut = np.array([5, -5, 4, -4, 3, -3, 2, -2, 1, -1, 0], dtype=np.int32)

  def __call__(self, hit):
    if hit.type == kCSC:
      clct = int(hit.pattern)
      bend = self.lut[clct]
      bend *= hit.endcap
    #elif hit.type == kGEM:
    #  bend = np.int32(hit.bend)
    #  bend *= hit.endcap
    elif hit.type == kME0:
      bend = np.int32(hit.bend)
    elif hit.type == kDT:
      bend = np.int32(hit.bend)
    else:  # kRPC, kGEM
      bend = np.int32(0)
    return bend

# Decide EMTF hit quality
# Currently using F/R for CSC & ME0, time for RPC
class EMTFQuality(object):
  def __call__(self, hit):
    emtf_qual = np.int32(0)
    if hit.type == kCSC or hit.type == kME0:
      if int(hit.fr) == 1:
        emtf_qual += 1
      else:
        emtf_qual -= 1
    elif hit.type == kRPC:
      emtf_qual += hit.time
      emtf_qual = np.clip(emtf_qual, -8, 7)
    else:  # kGEM, kDT
      pass
    return emtf_qual

# Decide EMTF hit quality (old version)
class EMTFOldQuality(object):
  def __call__(self, hit):
    emtf_qual = np.int32(hit.quality)
    if hit.type == kCSC or hit.type == kME0:
      # front chamber -> +1
      # rear chamber  -> -1
      if int(hit.fr) == 1:
        emtf_qual *= +1
      else:
        emtf_qual *= -1
    elif hit.type == kRPC or hit.type == kGEM:
      emtf_qual = np.int32(0)
    else:  # kDT
      pass
    return emtf_qual

# Decide EMTF hit time (integer unit)
class EMTFTime(object):
  def __call__(self, hit):
    emtf_time = np.int32(hit.time)
    emtf_time = np.clip(emtf_time, -8, 7)
    return emtf_time

# Functionalize
find_emtf_layer = EMTFLayer()
find_emtf_zones = EMTFZone()
find_emtf_zone_row = EMTFZoneRow()
find_emtf_zone_col = EMTFZoneCol()
find_emtf_phi = EMTFPhi()
find_emtf_old_phi = EMTFOldPhi()
find_emtf_theta = EMTFTheta()
find_emtf_bend = EMTFBend()
find_emtf_old_bend = EMTFOldBend()
find_emtf_qual = EMTFQuality()
find_emtf_old_qual = EMTFOldQuality()
find_emtf_time = EMTFTime()

# Decide whether hit is very legal and very cool
def is_emtf_legit_hit(hit):
  def check_bx(hit):
    if hit.type == kCSC:
      return hit.bx in (-1,0)
    elif hit.type == kDT:
      return hit.bx in (-1,0,+1)
    else:
      return hit.bx == 0
  def check_phi(hit):
    if hit.type == kME0:
      return hit.emtf_phi > 0
    elif hit.type == kDT:
      return hit.emtf_phi > 0
    else:
      return True
  return check_bx(hit) and check_phi(hit)

# Decide whether hit is used for Run 2 EMTF
def is_valid_for_run2(hit):
  is_csc = (hit.type == kCSC)
  is_rpc = (hit.type == kRPC)
  is_irpc = (hit.type == kRPC) and ((hit.station == 3 or hit.station == 4) and hit.ring == 1)
  is_omtf = (hit.type == kRPC) and ((hit.station == 1 or hit.station == 2) and hit.ring == 3)
  return (is_csc or (is_rpc and not is_irpc and not is_omtf))

# Decide the zone which the particle is belong to
def find_particle_zone(part):
  etastar = calc_etastar_from_eta(part.eta, part.phi, part.vx, part.vy, part.vz)
  ind = np.searchsorted(emtf_eta_bins, np.abs(etastar))
  return (4 - ind)  # ind = (1,2,3,4) -> zone (3,2,1,0)
