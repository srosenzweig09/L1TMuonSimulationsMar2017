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
  global donotdelete

  # ____________________________________________________________________________
  # Rates


  infiles = [
      "rateplots_mc_r305310_run2_bxWindow_all.root",
      "rateplots_mc_r305310_run2_dTheta_all.root",
      "rateplots_mc_r305310_run2_newPatt2_all.root",
      "rateplots_mc_r305310_run2_cscOnly_all.root",
  ]

  hnames = [
      "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt",
      "highest_emtf_absEtaMin0_absEtaMax2.5_qmin8_pt",
      "highest_emtf_absEtaMin0_absEtaMax2.5_qmin4_pt",
  ]

  for infile in infiles:
    for hname in hnames:

      #infile = "rateplots_mc_r305310_run2_all.0.root"
      #infile = "rateplots_mc_r305310_run2_bxWindow_all.root"
      #infile = "rateplots_mc_r305310_run2_dTheta_all.root"
      #infile = "rateplots_mc_r305310_run2_newPatt_all.root"
      #hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt"

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
      h.SetLineColor(632)  # kRed
      h.SetLineWidth(2)
      h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
      h.GetYaxis().SetTitle("Trigger rate (kHz)")
      h.GetYaxis().SetTitleOffset(1.3)
      histos[hname] = h
      tfile.Close()
      print h.GetBinContent(h.FindBin(0)), h.GetBinContent(h.FindBin(16)), h.GetBinContent(h.FindBin(22))

      # reference
      infile_ref = "rateplots_mc_r305310_run2_default_all.root"
      tfile = TFile.Open(infile_ref)
      h_nevents = tfile.Get("trackcounting/nevents")
      assert h_nevents != None, "Cannot get nevents"
      nevents = h_nevents.GetBinContent(2)
      print("nevents: {0}".format(nevents))
      h = tfile.Get("trackcounting/"+hname)
      assert h != None, "Cannot get %s" % hname
      h = h.Clone(hname + "_2")
      h.Sumw2()
      make_ptcut(h)
      make_rate(h, nevents)
      h.SetLineColor(632)  # kRed
      h.SetLineWidth(2)
      h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
      h.GetYaxis().SetTitle("Trigger rate (kHz)")
      h.GetYaxis().SetTitleOffset(1.3)
      histos[hname + "_ref"] = h
      tfile.Close()
      print h.GetBinContent(h.FindBin(0)), h.GetBinContent(h.FindBin(16)), h.GetBinContent(h.FindBin(22))

      # rate reduction
      cc1 = TCanvas("cc1", "cc1", 600, 700)
      cc1.Divide(1,2)
      cc1_1 = cc1.GetPad(1)
      cc1_1.SetPad(0.01,0.25,0.99,0.99)
      cc1_1.SetBottomMargin(0.01)
      cc1_1.SetGrid()
      cc1_1.SetLogy()
      cc1_2 = cc1.GetPad(2)
      cc1_2.SetPad(0.01,0.01,0.99,0.25)
      cc1_2.SetTopMargin(0.01)
      cc1_2.SetBottomMargin(0.43)
      cc1_2.SetGrid()


      h = histos[hname + "_ref"]
      denom = h.Clone("denom")
      h = histos[hname]
      numer = h.Clone("numer")

      denom.SetMinimum(0.02)
      denom.SetMaximum(5e3)
      denom.SetMarkerColor(1)  # kBlack
      denom.SetLineColor(1)  # kBlack

      ratio = numer.Clone("ratio")
      ratio.Divide(numer, denom, 1, 1, "")
      ratio.SetMinimum(0.5)
      ratio.SetMaximum(1.5)
      ratio.GetYaxis().SetTitle("ratio")

      denom.SetLabelSize(0.0);
      denom.GetXaxis().SetTitleSize(0.00)
      denom.GetYaxis().SetLabelSize(0.05)
      denom.GetYaxis().SetTitleSize(0.06)
      denom.GetYaxis().SetTitleOffset(1.10)
      ratio.GetXaxis().SetLabelSize(0.15)
      ratio.GetXaxis().SetTitleSize(0.18)
      ratio.GetXaxis().SetTitleOffset(1.10)
      ratio.GetYaxis().SetLabelSize(0.14)
      ratio.GetYaxis().SetTitleSize(0.18)
      ratio.GetYaxis().SetTitleOffset(0.37)
      ratio.GetYaxis().SetNdivisions(502)
      ratio.GetYaxis().SetLabelOffset(0.01)
      draw_option = "hist"

      cc1_1.cd()
      denom.Draw(draw_option)
      numer.Draw(draw_option + " same")

      cc1_2.cd()
      ratio.Draw(draw_option)
      xmin, xmax = ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax()
      tline.DrawLine(xmin, 1.0, xmax, 1.0)
      #ratio.Draw("same")
      ratio.Draw(draw_option + " same")

      cc1.cd()
      imgname = "rate_reduction"
      if "qmin12" in hname:
        imgname += "_qmin12"
      elif "qmin8" in hname:
        imgname += "_qmin8"
      elif "qmin4" in hname:
        imgname += "_qmin4"

      if "bxWindow" in infile:
        imgname += "_bxWindow"
      elif "dTheta" in infile:
        imgname += "_dTheta"
      elif "newPatt" in infile:
        imgname += "_newPatt"
      elif "cscOnly" in infile:
        imgname += "_cscOnly"
      gPad.Print("%s.png" % (options.outdir+imgname))

      donotdelete += [denom, numer, ratio]



# ______________________________________________________________________________
if __name__ == "__main__":
  import argparse
  import os

  # Setup argument parser
  parser = argparse.ArgumentParser()
  parser.add_argument("--infile", default="histos.root", help="input file (default: %(default)s)")
  parser.add_argument("--outdir", default="plots_compare3/", help="output directory (default: %(default)s)")
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
  gROOT.SetBatch(True)
  gStyle.SetNdivisions(510, "Y")
  gStyle.SetEndErrorSize(0)
  gStyle.SetPadGridX(True)
  gStyle.SetPadGridY(True)
  gStyle.SetPadRightMargin(0.05)
  gStyle.SetTitleOffset(1.5, "Y")
  TH1.AddDirectory(False)

  histos = {}
  #for k in options.tfile.GetListOfKeys():
  #  h = k.ReadObj()
  #  if h.ClassName() in ["TH1F", "TH2F", "TProfile", "TEfficiency"]:
  #    hname = h.GetName()
  #    histos[hname] = h

  # Call the main function
  main(histos, options)
