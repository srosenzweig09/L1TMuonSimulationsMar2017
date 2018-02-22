import numpy as np
np.random.seed(2023)

from itertools import izip
from rootpy.plotting import Hist, Hist2D, Graph, Efficiency
from rootpy.tree import Tree, TreeChain, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open


# ______________________________________________________________________________
# Analyzer

# Enums
kDT, kCSC, kRPC, kGEM, kTT = 0, 1, 2, 3, 20

# ______________________________________________________________________________
# Classes

import mmap
import struct

class EMTFTrack(object):
  def __init__(self):
    pass

class EMTFPtLUT(object):
  def __init__(self):
    pass

class EMTFPtAssignment(object):
  def __init__(self):
    self._load_ptlut()

  def close(self):
    self._unload_ptlut()

  def _load_ptlut(self):
    self.ptlut_path = "LUT_v07_07June17.dat"
    self.ptlut_file = open(self.ptlut_path, "rb")
    self.ptlut_mmap = mmap.mmap(self.ptlut_file.fileno(), 0, prot=mmap.PROT_READ)

  def _unload_ptlut(self):
    self.ptlut_mmap.close()
    self.ptlut_file.close

  def convert_to_gmt_pt(self, xml_pt):
    # Check MakePtLUT::makeLUT() in
    #   L1Trigger/L1TMuonEndCap/test/tools/MakePtLUT.cc
    if xml_pt < 0.:
      pt = 1.
    else:
      pt = xml_pt
    #
    max_pt = min(20., pt)
    pt_scale = 1.2 / (1 - 0.015*max_pt)
    pt *= pt_scale
    #
    gmt_pt = (pt * 2) + 1
    gmt_pt = int(gmt_pt)
    if gmt_pt > 511:
      gmt_pt = 511
    else:
      gmt_pt = gmt_pt
    return gmt_pt

  def convert_to_xml_pt(self, gmt_pt):
    # Check MakePtLUT::makeLUT() in
    #   L1Trigger/L1TMuonEndCap/test/tools/MakePtLUT.cc
    if gmt_pt <= 0:
      pt = 0
    else:
      pt = (gmt_pt-1) * 0.5
    #
    pt_unscale = 1 / (1.2 + 0.015*pt)
    pt_unscale = max(pt_unscale, (1 - 0.015*20)/1.2)
    #
    xml_pt = pt
    xml_pt *= pt_unscale
    return xml_pt

  # ____________________________________________________________________________
  def lookup(self, ptlut_addr):
    mm = self.ptlut_mmap
    mm.seek((ptlut_addr//2)*4)  # each address contains two pt values, stored as 32-bit (4-byte)
    bs = mm.read(4)             # returns 32-bit (4-byte) binary string
    pt_word = struct.unpack("<I", bs) # little-endian (x86), unsigned int
    pt_word = pt_word[0]              # result is a tuple even if it contains exactly one item
    if (ptlut_addr%2) == 0:     # low bit of address selects value (stored in 9 bits)
      pt_value = pt_word & 0x1FF
    else:
      pt_value = (pt_word >> 9) & 0x1FF
    xml_pt = self.convert_to_xml_pt(pt_value)
    return xml_pt

  # ____________________________________________________________________________
  def make_track(self, hits):
    track = None

    for hit in hits:
      sort_code = 0
      if hit.station == 1 and (hit.ring == 1 or hit.ring == 4):
        sort_code = 50
        if hit.type == kCSC:
          sort_code += 3
        elif hit.type == kGEM:
          sort_code += 2
        elif hit.type == kRPC:
          sort_code += 1
      elif hit.station == 1 and hit.ring == 2:
        sort_code = 10
        if hit.type == kCSC:
          sort_code += 3
        elif hit.type == kGEM:
          sort_code += 2
        elif hit.type == kRPC:
          sort_code += 1
      elif hit.station == 2 or hit.station == 3 or hit.station == 4:
        sort_code = hit.station * 10
        if hit.type == kCSC:
          sort_code += 3
        elif hit.type == kGEM:
          sort_code += 2
        elif hit.type == kRPC:
          sort_code += 1
      elif hit.station == 1:
        pass
      else:
        assert(False)
      hit.sort_code = sort_code

    selected_hits = [None, None, None, None]
    for station in (1,2,3,4):
      filtered_hits = None
      if station == 1:
        filtered_hits = filter(lambda x: 50 <= x.sort_code <= 59, hits)
        if filtered_hits:
          selected = max(filtered_hits, key=lambda x: x.sort_code)
          selected_hits[station-1] = selected
        else:
          filtered_hits = filter(lambda x: 10 <= x.sort_code <= 19, hits)
          if filtered_hits:
            selected = max(filtered_hits, key=lambda x: x.sort_code)
            selected_hits[station-1] = selected
      else:
        filtered_hits = filter(lambda x: (station * 10) <= x.sort_code <= (station * 10 + 9), hits)
        if filtered_hits:
          selected = max(filtered_hits, key=lambda x: x.sort_code)
          selected_hits[station-1] = selected


    # __________________________________________________________________________
    track = EMTFTrack()
    track.hits = [hit for hit in selected_hits if hit is not None]

    if not track.hits:  # exit if no hits
      return track

    ptlut_data = EMTFPtLUT()
    ptlut_data.delta_ph = [8191,8191,8191,8191,8191,8191]
    ptlut_data.sign_ph  = [1,1,1,1,1,1]
    ptlut_data.delta_th = [127,127,127,127,127,127]
    ptlut_data.sign_th  = [1,1,1,1,1,1]
    ptlut_data.cpattern = [0,0,0,0]
    ptlut_data.fr       = [0,0,0,0]

    ipair = 0
    for ist1 in (0,1,2):
      for ist2 in (1,2,3):
        if ist2 <= ist1:
          continue
        if selected_hits[ist1] is not None:
          if selected_hits[ist2] is not None:
            hit1 = selected_hits[ist1]
            hit2 = selected_hits[ist2]
            ptlut_data.delta_ph[ipair] = abs(hit1.emtf_phi - hit2.emtf_phi)
            ptlut_data.sign_ph[ipair] = (hit1.emtf_phi <= hit2.emtf_phi)
            ptlut_data.delta_th[ipair] = abs(hit1.emtf_theta - hit2.emtf_theta)
            ptlut_data.sign_th[ipair] = (hit1.emtf_theta <= hit2.emtf_theta)
        ipair += 1

    for ist in (0,1,2,3):
      if selected_hits[ist] is not None:
        hit = selected_hits[ist]
        ptlut_data.cpattern[ist] = hit.pattern
        ptlut_data.fr[ist] = hit.fr

    ptlut_data.st1_ring2 = 0
    for hit in track.hits:
      if hit.station == 1 and (hit.ring == 2 or hit.ring == 3):
        ptlut_data.st1_ring2 = 1
    track.ptlut_data = ptlut_data

    track.mode = 0
    for hit in track.hits:
      station = hit.station
      track.mode |= (1 << (4 - station))

    track.theta = np.median([hit.emtf_theta for hit in track.hits])
    track.theta = int(track.theta)
    return track

  def getNLBdPhiBin(self, dPhi, bits, max_):
    if not hasattr(self, 'dPhiNLBMap_4bit_256Max'):
      self.dPhiNLBMap_4bit_256Max = (
        0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 25, 31, 46, 68, 136
      )
      assert(len(self.dPhiNLBMap_4bit_256Max) == 16)
    if not hasattr(self, 'dPhiNLBMap_5bit_256Max'):
      self.dPhiNLBMap_5bit_256Max = (
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        16, 17, 19, 20, 21, 23, 25, 28, 31, 34, 39, 46, 55, 68, 91, 136
      )
      assert(len(self.dPhiNLBMap_5bit_256Max) == 32)
    if not hasattr(self, 'dPhiNLBMap_7bit_512Max'):
      self.dPhiNLBMap_7bit_512Max = (
          0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
         16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
         32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
         48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
         64,  65,  66,  67,  68,  69,  71,  72,  73,  74,  75,  76,  77,  79,  80,  81,
         83,  84,  86,  87,  89,  91,  92,  94,  96,  98, 100, 102, 105, 107, 110, 112,
        115, 118, 121, 124, 127, 131, 135, 138, 143, 147, 152, 157, 162, 168, 174, 181,
        188, 196, 204, 214, 224, 235, 247, 261, 276, 294, 313, 336, 361, 391, 427, 470
      )
      assert(len(self.dPhiNLBMap_7bit_512Max) == 128)

    assert((bits == 4 and max_ == 256) or (bits == 5 and max_ == 256) or (bits == 7 and max_ == 512))
    dPhiBin_ = (1 << bits) - 1
    sign_ = 1
    if dPhi < 0:
      sign_ = -1
    dPhi = sign_ * dPhi

    if max_ == 256:
      if bits == 4:
        for edge in xrange((1 << bits) - 1):
          if (self.dPhiNLBMap_4bit_256Max[edge] <= dPhi and self.dPhiNLBMap_4bit_256Max[edge+1] > dPhi):
            dPhiBin_ = edge
            break
      elif bits == 5:
        for edge in xrange((1 << bits) - 1):
          if (self.dPhiNLBMap_5bit_256Max[edge] <= dPhi and self.dPhiNLBMap_5bit_256Max[edge+1] > dPhi):
            dPhiBin_ = edge
            break
    elif max_ == 512:
      if bits == 7:
        for edge in xrange((1 << bits) - 1):
          if (self.dPhiNLBMap_7bit_512Max[edge] <= dPhi and self.dPhiNLBMap_7bit_512Max[edge+1] > dPhi):
            dPhiBin_ = edge
            break

    assert(dPhiBin_ >= 0 and dPhiBin_ < (1 << bits))
    return dPhiBin_

  def getdTheta(self, dTheta, bits):
    assert(bits == 2 or bits == 3)

    if bits == 2:
      if   (abs(dTheta) <= 1): dTheta_ = 2
      elif (abs(dTheta) <= 2): dTheta_ = 1
      elif (dTheta <= -3)    : dTheta_ = 0
      else                   : dTheta_ = 3

    elif bits == 3:
      if   (dTheta <= -4): dTheta_ = 0
      elif (dTheta == -3): dTheta_ = 1
      elif (dTheta == -2): dTheta_ = 2
      elif (dTheta == -1): dTheta_ = 3
      elif (dTheta ==  0): dTheta_ = 4
      elif (dTheta == +1): dTheta_ = 5
      elif (dTheta == +2): dTheta_ = 6
      else               : dTheta_ = 7

    assert(dTheta_ >= 0 and dTheta_ < (1 << bits))
    return dTheta_

  def get8bMode15(self, theta, st1_ring2, endcap, sPhiAB, clctA, clctB, clctC, clctD):
    if st1_ring2:
      theta = (min(max(theta, 46), 87) - 46) / 7
    else:
      theta = (min(max(theta,  5), 52) -  5) / 6
    assert(0 <= theta < 10)

    clctA_2b = self.getCLCT(clctA, endcap, sPhiAB, 2)
    nRPC = (clctA == 0) + (clctB == 0) + (clctC == 0) + (clctD == 0)

    if st1_ring2:
      if   (nRPC >= 2 and clctA == 0 and clctB == 0): rpc_word =  0
      elif (nRPC >= 2 and clctA == 0 and clctC == 0): rpc_word =  1
      elif (nRPC >= 2 and clctA == 0 and clctD == 0): rpc_word =  2
      elif (nRPC == 1 and clctA == 0               ): rpc_word =  3
      elif (nRPC >= 2 and clctD == 0 and clctB == 0): rpc_word =  4
      elif (nRPC >= 2 and clctD == 0 and clctC == 0): rpc_word =  8
      elif (nRPC >= 2 and clctB == 0 and clctC == 0): rpc_word = 12
      elif (nRPC == 1 and clctD == 0               ): rpc_word = 16
      elif (nRPC == 1 and clctB == 0               ): rpc_word = 20
      elif (nRPC == 1 and clctC == 0               ): rpc_word = 24
      else                                          : rpc_word = 28
      rpc_clct  = rpc_word + clctA_2b
      mode15_8b = (theta*32) + rpc_clct + 64
    else:
      if   (theta >= 4 and clctD == 0): rpc_word = 0
      elif (theta >= 4 and clctC == 0): rpc_word = 1
      elif (theta >= 4               ): rpc_word = 2
      else                            : rpc_word = 3
      rpc_clct  = rpc_word*4 + clctA_2b
      mode15_8b = ((theta % 4)*16) + rpc_clct

    assert(mode15_8b >= 0 and mode15_8b < (1 << 8))
    return mode15_8b

  def get2bRPC(self, clctA, clctB, clctC):
    if   (clctA == 0): rpc_2b = 0
    elif (clctC == 0): rpc_2b = 1
    elif (clctB == 0): rpc_2b = 2
    else             : rpc_2b = 3

    assert(rpc_2b >= 0 and rpc_2b < 4)
    return rpc_2b

  def getCLCT(self, clct, endcap, dPhiSign, bits):
    assert((0 <= clct <= 10) and abs(endcap) == 1 and abs(dPhiSign) == 1 and (bits == 2 or bits == 3))

    clct_ = 0
    sign_ = -1 * endcap * dPhiSign

    if bits == 2:
      if   clct == 10: clct_ = 1
      elif clct ==  9: clct_ = np.where(sign_ > 0 , 1 , 2)
      elif clct ==  8: clct_ = np.where(sign_ > 0 , 2 , 1)
      elif clct ==  7: clct_ = np.where(sign_ > 0 , 0 , 3)
      elif clct ==  6: clct_ = np.where(sign_ > 0 , 3 , 0)
      elif clct ==  5: clct_ = np.where(sign_ > 0 , 0 , 3)
      elif clct ==  4: clct_ = np.where(sign_ > 0 , 3 , 0)
      elif clct ==  3: clct_ = np.where(sign_ > 0 , 0 , 3)
      elif clct ==  2: clct_ = np.where(sign_ > 0 , 3 , 0)
      elif clct ==  1: clct_ = np.where(sign_ > 0 , 0 , 3)
      elif clct ==  0: clct_ = 0
      else           : clct_ = 1

    elif bits == 3:
      if   clct == 10: clct_ = 4
      elif clct ==  9: clct_ = np.where(sign_ > 0 , 3 , 5)
      elif clct ==  8: clct_ = np.where(sign_ > 0 , 5 , 3)
      elif clct ==  7: clct_ = np.where(sign_ > 0 , 2 , 6)
      elif clct ==  6: clct_ = np.where(sign_ > 0 , 6 , 2)
      elif clct ==  5: clct_ = np.where(sign_ > 0 , 1 , 7)
      elif clct ==  4: clct_ = np.where(sign_ > 0 , 7 , 1)
      elif clct ==  3: clct_ = np.where(sign_ > 0 , 1 , 7)
      elif clct ==  2: clct_ = np.where(sign_ > 0 , 7 , 1)
      elif clct ==  1: clct_ = np.where(sign_ > 0 , 1 , 7)
      elif clct ==  0: clct_ = 0
      else           : clct_ = 4

    assert(clct_ >= 0 and clct_ < (1 << bits))
    return clct_

  def getTheta(self, theta, st1_ring2, bits):
    assert((5 <= theta < 128) and (st1_ring2 == 0 or st1_ring2 == 1) and (bits == 4 or bits == 5))

    # For use in mode 15
    if bits == 4:
      if st1_ring2 == 0:
        theta_ = (min(max(theta, 5), 52) - 5) / 6
      elif st1_ring2 == 1:
        theta_ = ((min(max(theta, 46), 87) - 46) / 7) + 8

    # For use in all 2- and 3-station modes (all modes except 15)
    elif bits == 5:
      if st1_ring2 == 0:
        theta_ = (max(theta, 1) - 1) / 4
      elif st1_ring2 == 1:
        theta_ = ((min(theta, 104) - 1) / 4) + 6

    assert(theta_ >= 0 and ((bits == 4 and theta_ <= 13) or (bits == 5 and theta_ < (1<<bits))))
    return theta_


  def calculate_address(self, track):
    # from L1Trigger/L1TMuonEndCap/src/PtAssignmentEngine2017.cc
    address = 0

    mode = track.mode
    theta = track.theta
    endcap = track.hits[0].endcap
    nhits = len(track.hits)
    data = track.ptlut_data
    st1_ring2 = data.st1_ring2

    # __________________________________________________________________________
    # 'A' is first station in the track, 'B' the second, etc

    # Indices for quantities by station or station pair
    if   mode == 15:
      mode_ID, iA, iB, iC, iD = 0b001, 0, 1, 2, 3
    elif mode == 14:
      mode_ID, iA, iB, iC, iD = 0b011, 0, 1, 2, 99
    elif mode == 13:
      mode_ID, iA, iB, iC, iD = 0b010, 0, 1, 3, 99
    elif mode == 11:
      mode_ID, iA, iB, iC, iD = 0b001, 0, 2, 3, 99
    elif mode ==  7:
      mode_ID, iA, iB, iC, iD = 0b001, 1, 2, 3, 99
    elif mode == 12:
      mode_ID, iA, iB, iC, iD = 0b111, 0, 1, 99, 99
    elif mode == 10:
      mode_ID, iA, iB, iC, iD = 0b110, 0, 2, 99, 99
    elif mode ==  9:
      mode_ID, iA, iB, iC, iD = 0b101, 0, 3, 99, 99
    elif mode ==  6:
      mode_ID, iA, iB, iC, iD = 0b100, 1, 2, 99, 99
    elif mode ==  5:
      mode_ID, iA, iB, iC, iD = 0b011, 1, 3, 99, 99
    elif mode ==  3:
      mode_ID, iA, iB, iC, iD = 0b010, 2, 3, 99, 99
    else:
      raise ValueError('Unexpected mode: %i' % mode)

    iAB = np.where((iA >= 0 and iB >= 0) , iA + iB - (iA == 0) , -1)
    iAC = np.where((iA >= 0 and iC >= 0) , iA + iC - (iA == 0) , -1)
    iAD = np.where((iA >= 0 and iD >= 0) , 2                   , -1)
    iBC = np.where((iB >= 0 and iC >= 0) , iB + iC             , -1)
    iCD = np.where((iC >= 0 and iD >= 0) , 5                   , -1)

    # Fill variable info from pT LUT data
    if nhits == 4:
      dPhiAB = data.delta_ph[iAB]
      dPhiBC = data.delta_ph[iBC]
      dPhiCD = data.delta_ph[iCD]
      sPhiAB = data.sign_ph[iAB]
      sPhiBC = (data.sign_ph[iBC] == sPhiAB)
      sPhiCD = (data.sign_ph[iCD] == sPhiAB)
      dTheta = data.delta_th[iAD] * np.where(data.sign_th[iAD] , 1 , -1)
      frA    = data.fr      [iA]
      clctA  = data.cpattern[iA]
      clctB  = data.cpattern[iB]
      clctC  = data.cpattern[iC]
      clctD  = data.cpattern[iD]
    elif nhits == 3:
      dPhiAB = data.delta_ph[iAB]
      dPhiBC = data.delta_ph[iBC]
      sPhiAB = data.sign_ph[iAB]
      sPhiBC = (data.sign_ph[iBC] == sPhiAB)
      dTheta = data.delta_th[iAC] * np.where(data.sign_th[iAC] , 1 , -1)
      frA    = data.fr      [iA]
      frB    = data.fr      [iB]
      clctA  = data.cpattern[iA]
      clctB  = data.cpattern[iB]
      clctC  = data.cpattern[iC]
    elif nhits == 2:
      dPhiAB = data.delta_ph[iAB]
      sPhiAB = data.sign_ph[iAB]
      dTheta = data.delta_th[iAB] * np.where(data.sign_th[iAB] , 1 , -1)
      frA    = data.fr      [iA]
      frB    = data.fr      [iB]
      clctA  = data.cpattern[iA]
      clctB  = data.cpattern[iB]
    else:
      raise ValueError('Unexpected nhits: %i' % nhits)

    # Convert variables to words for pT LUT address
    if nhits == 4:
      dPhiAB = self.getNLBdPhiBin ( dPhiAB, 7, 512 )
      dPhiBC = self.getNLBdPhiBin ( dPhiBC, 5, 256 )
      dPhiCD = self.getNLBdPhiBin ( dPhiCD, 4, 256 )
      dTheta = self.getdTheta     ( dTheta, 2 )
      mode15_8b = self.get8bMode15( theta, st1_ring2, endcap, np.where(sPhiAB == 1 , 1 , -1), clctA, clctB, clctC, clctD )
    elif nhits == 3:
      dPhiAB = self.getNLBdPhiBin( dPhiAB, 7, 512 )
      dPhiBC = self.getNLBdPhiBin( dPhiBC, 5, 256 )
      dTheta = self.getdTheta    ( dTheta, 3 )
      rpc_2b = self.get2bRPC     ( clctA, clctB, clctC ) # Have to use un-compressed CLCT words
      clctA  = self.getCLCT      ( clctA, endcap, np.where(sPhiAB == 1 , 1 , -1), 2 )
      theta  = self.getTheta     ( theta, st1_ring2, 5 )
    elif nhits == 2:
      dPhiAB = self.getNLBdPhiBin( dPhiAB, 7, 512 )
      dTheta = self.getdTheta    ( dTheta, 3 )
      clctA  = self.getCLCT      ( clctA, endcap, np.where(sPhiAB == 1 , 1 , -1), 3 )
      clctB  = self.getCLCT      ( clctB, endcap, np.where(sPhiAB == 1 , 1 , -1), 3 )
      theta  = self.getTheta     ( theta, st1_ring2, 5 )
    else:
      raise ValueError('Unexpected nhits: %i' % nhits)

    # Form the pT LUT address
    if nhits == 4:
      address |= (dPhiAB    & ((1<<7)-1)) << (0)
      address |= (dPhiBC    & ((1<<5)-1)) << (0+7)
      address |= (dPhiCD    & ((1<<4)-1)) << (0+7+5)
      address |= (sPhiBC    & ((1<<1)-1)) << (0+7+5+4)
      address |= (sPhiCD    & ((1<<1)-1)) << (0+7+5+4+1)
      address |= (dTheta    & ((1<<2)-1)) << (0+7+5+4+1+1)
      address |= (frA       & ((1<<1)-1)) << (0+7+5+4+1+1+2)
      address |= (mode15_8b & ((1<<8)-1)) << (0+7+5+4+1+1+2+1)
      address |= (mode_ID   & ((1<<1)-1)) << (0+7+5+4+1+1+2+1+8)
      assert((1 << 29) <= address < (1 << 30))
    elif nhits == 3:
      address |= (dPhiAB    & ((1<<7)-1)) << (0)
      address |= (dPhiBC    & ((1<<5)-1)) << (0+7)
      address |= (sPhiBC    & ((1<<1)-1)) << (0+7+5)
      address |= (dTheta    & ((1<<3)-1)) << (0+7+5+1)
      address |= (frA       & ((1<<1)-1)) << (0+7+5+1+3)

      bit = 0
      if mode != 7:
        address |= (frB     & ((1<<1)-1)) << (0+7+5+1+3+1)
        bit = 1

      address |= (clctA     & ((1<<2)-1)) << (0+7+5+1+3+1+bit)
      address |= (rpc_2b    & ((1<<2)-1)) << (0+7+5+1+3+1+bit+2)
      address |= (theta     & ((1<<5)-1)) << (0+7+5+1+3+1+bit+2+2)

      if mode != 7:
        address |= (mode_ID & ((1<<2)-1)) << (0+7+5+1+3+1+bit+2+2+5)
        assert((1 << 27) <= address < (1 << 29))
      else:
        address |= (mode_ID & ((1<<1)-1)) << (0+7+5+1+3+1+bit+2+2+5)
        assert((1 << 26) <= address < (1 << 27))
    elif nhits == 2:
      address |= (dPhiAB    & ((1<<7)-1)) << (0)
      address |= (dTheta    & ((1<<3)-1)) << (0+7)
      address |= (frA       & ((1<<1)-1)) << (0+7+3)
      address |= (frB       & ((1<<1)-1)) << (0+7+3+1)
      address |= (clctA     & ((1<<3)-1)) << (0+7+3+1+1)
      address |= (clctB     & ((1<<3)-1)) << (0+7+3+1+1+3)
      address |= (theta     & ((1<<5)-1)) << (0+7+3+1+1+3+3)
      address |= (mode_ID   & ((1<<3)-1)) << (0+7+3+1+1+3+3+5)
      assert((1 << 24) <= address < (1 << 26))
    else:
      raise ValueError('Unexpected nhits: %i' % nhits)

    return address

  def calculate_pt(self, address):
    xml_pt = self.lookup(address)
    return xml_pt


# ______________________________________________________________________________
# Settings

# Get number of events
#maxEvents = -1
#maxEvents = 2000000
maxEvents = 20

use_condor = False

analysis = "verbose"

infile_r = None  # input file handle

def load_pgun():
  global infile_r
  #infile = 'ntuple_SingleMuon_Toy_5GeV_add.3.root'
  infile = '/tmp/jiafu/ntuple_SingleMuon_Toy_2GeV.0.root'
  infile_r = root_open(infile)
  tree = infile_r.ntupler.tree
  #tree = TreeChain('ntupler/tree', [infile])
  print('[INFO] Opening file: %s' % infile)

  # Define collection
  tree.define_collection(name='hits', prefix='vh_', size='vh_size')
  tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
  tree.define_collection(name='particles', prefix='vp_', size='vp_size')
  return tree

def unload_tree():
  global infile_r
  infile_r.close()


# ______________________________________________________________________________
# Analysis: verbose
if analysis == "verbose":
  tree = load_pgun()

  # Workers
  emtfptassign = EMTFPtAssignment()
  #for i in xrange(10):
  #  v = emtfptassign.lookup(i)
  #  print i, v
  #emtfptassign.close()

  # Loop over events
  for ievt, evt in enumerate(tree):
    if maxEvents != -1 and ievt == maxEvents:
      break

    print("Processing event: {0}".format(ievt))

    # Hits
    for ihit, hit in enumerate(evt.hits):
      print(".. hit  {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(ihit, hit.bx, hit.type, hit.station, hit.ring, hit.sector, hit.fr, hit.sim_phi, hit.sim_theta, hit.sim_tp1, hit.sim_tp2))
    # Tracks
    for itrk, trk in enumerate(evt.tracks):
      print(".. trk  {0} {1} {2} {3} {4} {5} {6} {7}".format(itrk, trk.sector, trk.mode, trk.pt, trk.phi, trk.eta, trk.theta, trk.q))
      trk_xml_pt = emtfptassign.convert_to_xml_pt(emtfptassign.convert_to_gmt_pt(trk.xml_pt))
      trk_xml_pt_check = emtfptassign.lookup(trk.address)
      #print(".. .. {0} {1} {2}".format(trk.xml_pt, trk_xml_pt, trk_xml_pt_check))
      print(".. ..  {0} {1} {2:10d} {3:030b} {4}".format(itrk, trk.mode, trk.address, trk.address, trk_xml_pt_check))
      assert(np.isclose(trk_xml_pt, trk_xml_pt_check))
    # Gen particles
    for ipart, part in enumerate(evt.particles):
      print(".. part {0} {1} {2} {3} {4} {5}".format(ipart, part.pt, part.phi, part.eta, part.theta, part.q))


    # __________________________________________________________________________
    # Make track

    track = emtfptassign.make_track(evt.hits)
    if track.hits:
      address = emtfptassign.calculate_address(track)
      xml_pt = emtfptassign.calculate_pt(address)
      print(".. ..  {0} {1} {2:10d} {3:030b} {4}".format(0, track.mode, address, address, xml_pt))
    else:
      print(".. no new trk")


  # End loop over events
  unload_tree()

  # Close workers
  emtfptassign.close()
