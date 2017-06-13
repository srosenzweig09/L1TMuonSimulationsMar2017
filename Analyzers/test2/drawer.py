#!/usr/bin/env python

from ROOT import TFile, TH1, TH1F, TH2F, TProfile, TGraphAsymmErrors, TEfficiency, TColor, TLine, TLegend, TLatex, gROOT, gStyle, gPad

import os
import argparse


class MyDrawer(object):

  def __init__(self):
    # Setup argument parser
    self.parser = argparse.ArgumentParser()
    #self.parser.add_argument("--infile", default="histos.root", help="input file (default: %(default)s)")
    self.parser.add_argument("--outdir", default="figures/", help="output directory (default: %(default)s)")
    self.parser.add_argument("-v", "--verbose", action="count", default=0, help="verbosity (default: %(default)s)")
    self.options = self.parser.parse_args()
    if not self.options.outdir.endswith("/"):
      self.options.outdir += "/"
    if not os.path.exists(self.options.outdir):
      os.makedirs(self.options.outdir)
    #if not os.path.isfile(self.options.infile):
    #  raise Exception("File does not exist: %s" % self.options.infile)
    #self.options.tfile = TFile.Open(self.options.infile)

    self.options.palette = ("#333333", "#377eb8", "#e41a1c", "#984ea3", "#ff7f00", "#4daf4a")

    # Setup basic drawer
    gROOT.LoadMacro("tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle();")
    #gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    gROOT.SetBatch(True)
    gStyle.SetNdivisions(510, "Y")
    gStyle.SetEndErrorSize(0)
    gStyle.SetPadGridX(True)
    gStyle.SetPadGridY(True)
    #gStyle.SetOptStat(111110)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetTitleOffset(1.1, "Y")
    gStyle.SetPalette(57)  # kBird
    gStyle.SetNumberContours(50)
    
