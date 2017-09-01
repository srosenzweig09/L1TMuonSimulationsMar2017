from ROOT import TFile, TTree, TH1F
from itertools import izip, count


# ______________________________________________________________________________
# Configurables

filename = "L1Ntuple_2017.0.root"

maxEvents = -1

tfile = TFile.Open(filename)
tree_event = tfile.Get("l1EventTree/L1EventTree")
tree_data = tfile.Get("l1UpgradeTree/L1UpgradeTree")
tree_emul = tfile.Get("l1UpgradeEmuTree/L1UpgradeTree")
tree_data_tfmuon = tfile.Get("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree")
tree_emul_tfmuon = tfile.Get("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree")
assert(tree_event)
assert(tree_data)
assert(tree_emul)
assert(tree_data_tfmuon)
assert(tree_emul_tfmuon)


# ______________________________________________________________________________
# Functions

to_list = lambda x : [xx for xx in x]

class MyMuon:
  def __init__(self, bx, qual, mtf, pt, eta, phi):
    self.bx = bx
    self.qual = qual
    self.mtf = mtf
    self.pt = pt
    self.eta = eta
    self.phi = phi

  def __repr__(self):
    return "({0} {1} {2} {3} {4} {5})".format(self.bx, self.qual, self.mtf, self.pt, self.eta, self.phi)

  def __eq__(self, other):
    b = \
        self.bx == other.bx and \
        self.qual == other.qual and \
        self.mtf == other.mtf and \
        self.pt == other.pt and \
        self.eta == other.eta and \
        self.phi == other.phi
    return b

def make_muons(evt):
  muons = []
  for (bx, qual, mtf, pt, eta, phi) in izip(evt.muonBx, evt.muonQual, evt.muonTfMuonIdx, evt.muonIEt, evt.muonIEta, evt.muonIPhi):
    m = MyMuon(bx, qual, mtf, pt, eta, phi)
    muons.append(m)
  return muons

def filter_muons(muons):
  f = lambda m : (m.bx == 0) and (m.qual != 0) and ((0 <= m.mtf <= 17) or (90 <= m.mtf <= 107))  # EMTF
  #f = lambda m : True
  return [m for m in muons if f(m)]

def debug_gmt(evt, label):
  for (
      bx, qual, mtf, pt, eta, phi,
      q, iso, fpt, feta, fphi, vtxeta, vtxphi, deta, dphi, i
  ) in izip(
      evt.muonBx, evt.muonQual, evt.muonTfMuonIdx, evt.muonIEt, evt.muonIEta, evt.muonIPhi,
      evt.muonChg, evt.muonIso, evt.muonEt, evt.muonEta, evt.muonPhi, evt.muonIEtaAtVtx, evt.muonIPhiAtVtx, evt.muonIDEta, evt.muonIDPhi, count(),
  ):
    print "....", label, bx, qual, mtf, pt, eta, phi, "/", q, iso, fpt, feta, fphi, vtxeta, vtxphi, deta, dphi
  return

def debug_emtf(evt, label):
  emtf = evt.L1UpgradeEmtfMuon
  for (
      bx, qual, mtf, pt, eta, phi,
      q, vq, hf, glbphi, link, proc, i,
  ) in izip(
      emtf.tfMuonBx, emtf.tfMuonHwQual, emtf.tfMuonTrackFinderType, emtf.tfMuonHwPt, emtf.tfMuonHwEta, emtf.tfMuonHwPhi,
      emtf.tfMuonHwSign, emtf.tfMuonHwSignValid, emtf.tfMuonHwHF, emtf.tfMuonGlobalPhi, emtf.tfMuonLink, emtf.tfMuonProcessor, count(),
  ):
    trAdd0 = 0 if (emtf.tfMuonWh[i] >= 0) else 1
    trAdd1 = abs(emtf.tfMuonWh[i])
    trAdd2 = emtf.tfMuonTrAdd[4*i+0]
    trAdd3 = emtf.tfMuonTrAdd[4*i+1]
    trAdd4 = emtf.tfMuonTrAdd[4*i+2]
    trAdd5 = emtf.tfMuonTrAdd[4*i+3]
    print "....", label, bx, qual, mtf, pt, eta, phi, "/", q, vq, hf, glbphi, link, proc, trAdd0, trAdd1, trAdd2, trAdd3, trAdd4, trAdd5
  return


# ______________________________________________________________________________
# Event loop

summary = {}

print "[INFO] I am only checking (bx, qual, mtf, pt, eta, phi)"

for (ievt, (evt0, evt1, evt2, evt3, evt4)) in enumerate(izip(tree_event, tree_data, tree_emul, tree_data_tfmuon, tree_emul_tfmuon)):
  run_number, event_number, lumi_number = evt0.Event.run, evt0.Event.event, evt0.Event.lumi

  if (maxEvents != -1) and (ievt >= maxEvents):
    break

  muons1 = make_muons(evt1)  # GMT muons from unpacker
  muons2 = make_muons(evt2)  # GMT muons from emulator

  muons1_f = filter_muons(muons1)
  muons2_f = filter_muons(muons2)

  def compare(lhs, rhs):
    results = []
    for l, r in izip(lhs, rhs):
      diff = 0
      # not using diff = 1. reserved for more than one different track
      if l.bx != r.bx:
        diff = 2
      elif l.qual != r.qual:
        diff = 3
      elif l.mtf != r.mtf:
        diff = 4
      elif l.pt != r.pt:
        diff = 5
      elif l.eta != r.eta:
        diff = 6
      elif l.phi != r.phi:
        diff = 7
      results.append(diff)

    diff2 = 0
    for diff in results:
      if diff != 0:
        diff2 += 1

    diff3 = 0
    if diff2 > 1:
      # More than one different track
      diff3 = 1
    else:
      # Only one different track
      for diff in results:
        if diff != 0:
          diff3 = diff
    return diff3

  the_diff = compare(muons1_f, muons2_f)

  # Counter
  summary[the_diff] = summary.get(the_diff, 0) + 1
  summary[-1] = summary.get(-1, 0) + 1

  # Debug
  print ("Processing evt %i, Run %i, Event %i, LumiSection %i." % (ievt+1, run_number, event_number, lumi_number)), ("# muons: %i, %i" % (len(muons1_f), len(muons2_f)))
  print "..", muons1_f
  print "..", muons2_f

  if the_diff:
    print "[ERROR] diff code: %i" % the_diff

    # Debug
    debug_gmt(evt1, "gmt_unp")
    debug_gmt(evt2, "gmt_emu")

    # Debug
    debug_emtf(evt3, "emtf_unp")
    debug_emtf(evt4, "emtf_emu")

  print "-" * 60

# ______________________________________________________________________________
# Summary

print "Total difference: %i/%i" % (summary[-1] - summary[0], summary[-1])
for k, v in summary.iteritems():
  print "diff code %i: %i" % (k, v)

