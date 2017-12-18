#!/usr/bin/env python

from ROOT import TFile, TH1, TH1F, TH2F, TProfile, TGraphAsymmErrors, TEfficiency, TColor, TLine, TLegend, TLatex, TCanvas, TPad, TF1, gROOT, gStyle, gPad
from math import sqrt


# ______________________________________________________________________________
# Globals

palette = map(lambda x: TColor.GetColor(x), ("#8E8037", "#6299F4", "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00"))

donotdelete = []

# Legend
tlegend = TLegend(0.46,0.14,0.96,0.34)
tlegend.SetFillStyle(0)
tlegend.SetLineColor(0)
tlegend.SetShadowColor(0)
tlegend.SetBorderSize(0)

# Text
tlatex = TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(42)
tlatex.SetTextSize(0.035)

# Lines
tline = TLine()
tline.SetLineColor(1)

# Functions
def make_ptcut(h):
  use_overflow = True
  binsum = 0
  binerr2 = 0
  for ib in xrange(h.GetNbinsX()+2-1, 0-1, -1):
    if (not use_overflow) and (ib == 0 or ib == h.GetNbinsX()+1):
      continue
    binsum += h.GetBinContent(ib)
    binerr2 += h.GetBinError(ib)**2
    h.SetBinContent(ib, binsum)
    h.SetBinError(ib, sqrt(binerr2))
  return

def make_rate(h, nevents):
  orbitFreq = 11245.6
  nCollBunches = 1866
  nZeroBiasEvents = nevents
  convFactorToHz = orbitFreq * nCollBunches / nZeroBiasEvents
  h.Scale(convFactorToHz / 1000.)
  return


# ______________________________________________________________________________
# Main function
def main(histos, options):
  pass


# ______________________________________________________________________________
if __name__ == "__main__":
  import argparse
  import os

  # Setup argument parser
  parser = argparse.ArgumentParser()
  parser.add_argument("--infile", default="histos.root", help="input file (default: %(default)s)")
  parser.add_argument("--outdir", default="plots/", help="output directory (default: %(default)s)")
  parser.add_argument("-v", "--verbose", action="count", default=0, help="verbosity (default: %(default)s)")
  options = parser.parse_args()
  if not options.outdir.endswith("/"):
    options.outdir += "/"
  if not os.path.exists(options.outdir):
    os.makedirs(options.outdir)
  #if not os.path.isfile(options.infile):
  #  raise Exception("File does not exist: %s" % options.infile)
  #options.tfile = TFile.Open(options.infile)

  # Setup basic drawer
  gROOT.LoadMacro("tdrstyle.C")
  gROOT.ProcessLine("setTDRStyle();")
  gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
  #gROOT.SetBatch(True)
  gStyle.SetNdivisions(510, "Y")
  gStyle.SetEndErrorSize(0)
  gStyle.SetPadGridX(True)
  gStyle.SetPadGridY(True)
  gStyle.SetPadRightMargin(0.05)
  gStyle.SetTitleOffset(1.5, "Y")

  histos = {}
  #for k in options.tfile.GetListOfKeys():
  #  h = k.ReadObj()
  #  if h.ClassName() in ["TH1F", "TH2F", "TProfile", "TEfficiency"]:
  #    hname = h.GetName()
  #    histos[hname] = h

  # Call the main function
  #main(histos, options)


  # ____________________________________________________________________________
  if True:
    #infile = "rateplots_mc_r305310_run2_all.0.root"
    infile = "rateplots_mc_r305310_run2_bxWindow_all.root"
    #infile = "rateplots_mc_r305310_run2_newPatt_all.root"
    hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"

    tfile = TFile.Open(infile)
    h_nevents = tfile.Get("trackcounting/nevents")
    assert h_nevents != None, "Cannot get nevents"
    nevents = h_nevents.GetBinContent(2)
    print("nevents: {0}".format(nevents))
    h = tfile.Get("trackcounting/"+hname)
    assert h != None, "Cannot get %s" % hname
    h = h.Clone(hname + "_1")
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    h.SetLineColor(602)  # kBlue+2
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
    h.GetYaxis().SetTitle("Trigger rate (kHz)")
    h.GetYaxis().SetTitleOffset(1.3)
    histos[hname] = h

    # Draw
    h.Draw()
    gPad.SetLogy()
    if "bxWindow" in infile:
      gPad.Print(hname + "_bxWindow.png")
    elif "newPatt" in infile:
      gPad.Print(hname + "_newPatt.png")
    else:
      gPad.Print(hname + ".png")

    print h.GetBinContent(h.FindBin(0)), h.GetBinContent(h.FindBin(16)), h.GetBinContent(h.FindBin(22))


  # ____________________________________________________________________________
  if False:
    infile = "rateplots_mc_r305310_run2_all.1.root"
    #hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_mc_pt"
    hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_fake_pt"

    tfile = TFile.Open(infile)
    h_nevents = tfile.Get("minbiasmuonanalyzer/nevents")
    assert h_nevents != None, "Cannot get nevents"
    nevents = h_nevents.GetBinContent(2)
    print("nevents: {0}".format(nevents))
    h = tfile.Get("minbiasmuonanalyzer/"+hname)
    assert h != None, "Cannot get %s" % hname
    h = h.Clone(hname + "_1")
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    h.SetLineColor(602)  # kBlue+2
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
    h.GetYaxis().SetTitle("Trigger rate (kHz)")
    h.GetYaxis().SetTitleOffset(1.3)
    histos[hname] = h

    # Draw
    h.Draw()
    gPad.SetLogy()
    #gPad.Print(hname + ".png")

    print h.GetBinContent(h.FindBin(0)), h.GetBinContent(h.FindBin(16)), h.GetBinContent(h.FindBin(22))

