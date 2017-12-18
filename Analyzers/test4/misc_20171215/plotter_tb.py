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
  parser.add_argument("--infile", default="histos_tb.root", help="input file (default: %(default)s)")
  parser.add_argument("--outdir", default="plots_tb/", help="output directory (default: %(default)s)")
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
  gStyle.SetOptStat(1111110)

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
    hnames = [
      "hit_occupancy_csc_st1",
      "hit_occupancy_csc_st2",
      "hit_occupancy_csc_st3",
      "hit_occupancy_csc_st4",
    ]
    for hname in hnames:
      tfile = TFile.Open("histos_tb.root")
      h1 = tfile.Get("trackcounting/"+hname)
      h1.SetLineWidth(2)
      h1.SetLineColor(418)  # kGreen+2
      h1.SetMinimum(1.0)
      h1.Draw("hist")
      gPad.SetLogy()
      gPad.Print("plots_tb/"+hname+".png")
      tfile.Close()

    hnames = [
      "hit_occupancy_rpc_st1",
      "hit_occupancy_rpc_st2",
      "hit_occupancy_rpc_st3",
      "hit_occupancy_rpc_st4",
    ]
    for hname in hnames:
      tfile = TFile.Open("histos_tb.root")
      h1 = tfile.Get("trackcounting/"+hname)
      h1.SetLineWidth(2)
      h1.SetLineColor(600)  # kBlue
      h1.SetMinimum(1.0)
      h1.Draw("hist")
      gPad.SetLogy()
      gPad.Print("plots_tb/"+hname+".png")
      tfile.Close()


  if True:
    hnames = [
      "trk_eta_%s",
      "trk_mode_%s",
      "trk_hit_type1_%s",
      "trk_hit_type2_%s",
      "trk_hit_type3_%s",
      "trk_hit_type4_%s",
      "trk_hit_pattern1_%s",
      "trk_hit_pattern2_%s",
      "trk_hit_pattern3_%s",
      "trk_hit_pattern4_%s",
      "trk_hit_occu1_%s",
      "trk_hit_occu2_%s",
      "trk_hit_occu3_%s",
      "trk_hit_occu4_%s",
      "trk_dphi1_%s",
      "trk_dphi2_%s",
      "trk_dphi3_%s",
      "trk_dphi4_%s",
      "trk_dphi5_%s",
      "trk_dphi6_%s",
      "trk_dtheta1_%s",
      "trk_dtheta2_%s",
      "trk_dtheta3_%s",
      "trk_dtheta4_%s",
      "trk_dtheta5_%s",
      "trk_dtheta6_%s",
      "trk_bad_hit_%s",
    ]

    def get_signal_line():
      #return TColor.GetColor("#0000ee")
      return TColor.GetColor("#ff0000")

    def get_signal_fill():
      #return TColor.GetColor("#7d99d1")
      return TColor.GetColor("#ff0000")

    def get_background_line():
      #return TColor.GetColor("#ff0000")
      return TColor.GetColor("#000000")

    def get_background_fill():
      #return TColor.GetColor("#ff0000")
      return TColor.GetColor("#000000")

    for hname in hnames:
      tfile = TFile.Open("histos_tb.root")
      h1 = tfile.Get("trackcounting/" + (hname % "real"))
      h1.SetLineWidth(2)
      h1.SetLineColor(get_signal_line())
      h2 = tfile.Get("trackcounting/" + (hname % "fake"))
      h2.SetLineWidth(2)
      h2.SetLineColor(get_background_line())
      h1.SetStats(0)
      h1.SetMaximum(h1.GetMaximum() * 1.4)
      if hname == "trk_bad_hit_%s":  h1.SetMaximum(h2.GetMaximum() * 1.4)
      h1.Draw("hist")
      h2.Draw("hist same")
      gPad.Print("plots_tb/" + (hname % "real") + ".png")
      tfile.Close()
