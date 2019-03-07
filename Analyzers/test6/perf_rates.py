#!/usr/bin/env python

hname2026_f = lambda hname: "highest_emtf2026_" + hname[13:]

donotdelete = []

infile = "histos_tbb_add.25.root"

# Functions
from math import sqrt

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
  orbitFreq = 11246.
  nCollBunches = 2808  # for lumi=8e34, PU=200, xsec_pp=80mb
  nZeroBiasEvents = nevents
  convFactorToHz = orbitFreq * nCollBunches / nZeroBiasEvents
  h.Scale(convFactorToHz / 1000.)
  return


# ______________________________________________________________________________
if __name__ == '__main__':
  from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TH1F, TH2F, TPolyLine, TLatex, TColor, TEfficiency, TLine, TLatex, TGraph, TGraphAsymmErrors

  # ____________________________________________________________________________
  # Setup basic drawer
  gROOT.LoadMacro("tdrstyle.C")
  gROOT.ProcessLine("setTDRStyle();")
  #gROOT.SetBatch(True)
  gStyle.SetPadGridX(True)
  gStyle.SetPadGridY(True)
  gStyle.SetMarkerStyle(1)
  gStyle.SetEndErrorSize(0)
  gROOT.ForceStyle()

  tline = TLine()
  tline.SetLineColor(920+2)  # kGray+2
  tline.SetLineStyle(2)

  tlatexCMS1 = TLatex()
  tlatexCMS1.SetNDC()
  tlatexCMS1.SetTextFont(61)
  tlatexCMS1.SetTextSize(0.75*0.05)

  tlatexCMS2 = TLatex()
  tlatexCMS2.SetNDC()
  tlatexCMS2.SetTextFont(52)
  tlatexCMS2.SetTextSize(0.60*0.05)

  tlatexCMS3 = TLatex()
  tlatexCMS3.SetNDC()
  tlatexCMS3.SetTextFont(42)
  tlatexCMS3.SetTextSize(0.60*0.05)
  tlatexCMS3.SetTextAlign(11)

  tfile = TFile.Open(infile)
  h_nevents = tfile.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents = h_nevents.GetBinContent(2)


  # ____________________________________________________________________________
  # highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt

  hname_pairs = [
    ("highest_emtf_absEtaMin1.24_absEtaMax2.4_qmin12_pt", "emtf2026_rate_reduction"),
    ("highest_emtf_absEtaMin1.24_absEtaMax1.65_qmin12_pt", "emtf2026_rate_reduction_1"),
    ("highest_emtf_absEtaMin1.65_absEtaMax2.15_qmin12_pt", "emtf2026_rate_reduction_2"),
    ("highest_emtf_absEtaMin2.15_absEtaMax2.4_qmin12_pt", "emtf2026_rate_reduction_3"),
    ("highest_emtf_absEtaMin0.8_absEtaMax2.4_qmin12_pt", "emtf2026_rate_reduction_4"),
    ("highest_emtf_absEtaMin0.8_absEtaMax1.24_qmin12_pt", "emtf2026_rate_reduction_5"),
  ]

  for hname, imgname in hname_pairs:
    h = tfile.Get(hname)
    h = h.Clone(h.GetName() + "_clone")
    h.Sumw2()
    x0 = h.Integral()
    x1 = h.Integral(h.FindBin(20), h.FindBin(100))
    make_ptcut(h)
    make_rate(h, nevents)
    x2 = h.GetBinContent(h.FindBin(20))
    h.SetFillColor(632)  # kRed
    h.SetFillStyle(3003)
    h.SetLineColor(632)  # kRed
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
    h.GetYaxis().SetTitle("Trigger rate [kHz]")
    h.GetYaxis().SetTitleOffset(1.3)
    denom = h.Clone("denom")
    denom.SetMaximum(2e4)
    denom.SetMinimum(0.2)

    h = tfile.Get(hname2026_f(hname))
    h = h.Clone(h.GetName() + "_clone")
    h.Sumw2()
    x3 = h.Integral()
    x4 = h.Integral(h.FindBin(20), h.FindBin(100))
    make_ptcut(h)
    make_rate(h, nevents)
    x5 = h.GetBinContent(h.FindBin(20))
    h.SetFillColor(600)  # kBlue
    h.SetFillStyle(3003)
    h.SetLineColor(600)  # kBlue
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle("p_{T} threshold [GeV]")
    h.GetYaxis().SetTitle("Trigger rate [kHz]")
    h.GetYaxis().SetTitleOffset(1.3)
    numer = h.Clone("numer")

    ratio = numer.Clone("ratio")
    ratio.Divide(numer, denom, 1, 1, "")
    x6 = ratio.GetBinContent(ratio.FindBin(20))
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
    ratio.GetYaxis().SetTitle("ratio")
    print (x0, x1, x2), (x3, x4, x5), x6

    # ____________________________________________________________________________
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
    draw_err_option = "e3"

    cc1_1.cd()
    denom.SetStats(0)
    if imgname != "emtf2026_rate_reduction_5":
      denom.Draw(draw_err_option)
      numer.Draw(draw_err_option + " same")
      denom_no_fill = denom.Clone("denom_no_fill")
      denom_no_fill.SetFillStyle(0)
      denom_no_fill.Draw(draw_option + " same")
      numer_no_fill = numer.Clone("numer_no_fill")
      numer_no_fill.SetFillStyle(0)
      numer_no_fill.Draw(draw_option + " same")
    else:
      numer.Draw(draw_err_option)
      numer_no_fill = numer.Clone("numer_no_fill")
      numer_no_fill.SetFillStyle(0)
      numer_no_fill.Draw(draw_option + " same")

    cc1_2.cd()
    if imgname != "emtf2026_rate_reduction_5":
      ratio.SetStats(0)
      ratio.Draw(draw_err_option)
      ratio_no_fill = ratio.Clone("ratio_no_fill")
      ratio_no_fill.SetFillStyle(0)
      ratio_no_fill.Draw(draw_option + " same")
    else:
      for b in xrange(0, ratio.GetNbinsX()+2):
        ratio.SetBinContent(b, 1)
        ratio.SetBinError(b, 1e-6)
      ratio.SetStats(0)
      ratio.Draw(draw_err_option)
      ratio_no_fill = ratio.Clone("ratio_no_fill")
      ratio_no_fill.SetFillStyle(0)
      ratio_no_fill.Draw(draw_option + " same")

    xmin, xmax = ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax()
    tline.DrawLine(xmin, 1.0, xmax, 1.0)

    def draw_cms_lumi_1():
      tlatexCMS1.SetTextSize(0.75*0.05*1.1)
      tlatexCMS2.SetTextSize(0.60*0.05*1.1)
      tlatexCMS3.SetTextSize(0.60*0.05*1.1)
      tlatexCMS1.DrawLatex(0.164, 0.965, 'CMS')
      tlatexCMS2.DrawLatex(0.256, 0.965, 'Phase-2 Simulation')
      tlatexCMS3.DrawLatex(0.865, 0.965, '<PU>=200')

    cc1_1.cd()
    draw_cms_lumi_1()

    cc1.cd()
    gPad.Print("figures_perf/" + imgname + ".png")
    gPad.Print("figures_perf/" + imgname + ".pdf")

    if False:
      if imgname == "emtf2026_rate_reduction":
        outfile = TFile.Open("figures_perf/" + imgname + ".root", "RECREATE")
        for obj in [denom, numer, ratio, cc1]:
          obj.Write()
        outfile.Close()
