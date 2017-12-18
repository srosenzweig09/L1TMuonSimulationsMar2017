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


# ______________________________________________________________________________
# Main function
def main(histos, options):
  pass


# ______________________________________________________________________________
if __name__ == '__main__':
  import argparse
  import os

  # Setup argument parser
  parser = argparse.ArgumentParser()
  parser.add_argument("--infile", default="histos_pm.root", help="input file (default: %(default)s)")
  parser.add_argument("--outdir", default="plots_pm/", help="output directory (default: %(default)s)")
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
  if False:
    # muon_ptmin2_dxy
    tfile = TFile.Open("plots_pm/muon_ptmin2_dxy.root")
    h1 = tfile.Get("muon_ptmin2_dxy")
    h1.SetLineWidth(2)
    h1.SetLineColor(600)  # kBlue
    h1.Draw("hist")
    h1.SetAxisRange(0., 50., "X")
    gPad.SetLogy()
    gPad.Print("plots_pm/muon_ptmin2_dxy.png")
    tfile.Close()

    # muon_ptmin2_dz
    tfile = TFile.Open("plots_pm/muon_ptmin2_dz.root")
    h1 = tfile.Get("muon_ptmin2_dz")
    h1.SetLineWidth(2)
    h1.SetLineColor(600)  # kBlue
    h1.Draw("hist")
    h1.SetAxisRange(0., 50., "X")
    gPad.SetLogy()
    gPad.Print("plots_pm/muon_ptmin2_dz.png")
    tfile.Close()

    # muon_ptmin2_eta
    tfile = TFile.Open("plots_pm/muon_ptmin2_eta.root")
    h1 = tfile.Get("muon_ptmin2_eta")
    h1.SetLineWidth(2)
    h1.SetLineColor(600)  # kBlue
    h1.Draw("hist")
    h1.SetAxisRange(1.0, 2.5, "X")
    gPad.SetLogy(0)
    gPad.Print("plots_pm/muon_ptmin2_eta.png")
    tfile.Close()

    # muon_ptmin2_bx
    tfile = TFile.Open("plots_pm/muon_ptmin2_bx.root")
    h1 = tfile.Get("muon_ptmin2_bx")
    h1.SetLineWidth(2)
    h1.SetLineColor(600)  # kBlue
    h1.Draw("hist")
    #h1.SetAxisRange(0., 50., "X")
    gPad.SetLogy(0)
    gPad.Print("plots_pm/muon_ptmin2_bx.png")
    tfile.Close()

  if True:
    tfile = TFile.Open("histos_pm.root")
    for i in xrange(5):
      denom = tfile.Get("muon_pt_denom_zone%i" % i)
      numer = tfile.Get("muon_pt_fiducial_zone%i" % i)
      ratio = numer.Clone("ratio_zone%i" % i)
      rebin = 2
      if rebin != 1:
        denom.Rebin(rebin)
        numer.Rebin(rebin)
        ratio.Rebin(rebin)
      ratio.Divide(numer, denom, 1, 1, "b")
      ratio.SetStats(0)
      ratio.SetMinimum(0)
      ratio.SetMaximum(1.1)
      ratio.Draw()
      ratio.SetAxisRange(0., 20., "X")
      gPad.Print("plots_pm/muon_pt_fiducial_div_denom_zone%i.png" % i)
    tfile.Close()
