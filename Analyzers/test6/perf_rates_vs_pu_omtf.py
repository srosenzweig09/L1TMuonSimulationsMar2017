#!/usr/bin/env python

hname2026_f = lambda hname: "highest_emtf2026_" + hname[13:]

donotdelete = []

infile = "histos_tbb_add.25.root"
infile140 = "histos_tbb_140_add.25.root"
infile250 = "histos_tbb_250_add.25.root"
infile300 = "histos_tbb_300_add.25.root"

from perf_rates import make_ptcut, make_rate


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

  tlatex = TLatex()
  tlatex.SetTextFont(42)
  tlatex.SetTextSize(0.60*0.05)

  tfile200 = TFile.Open(infile)
  h_nevents = tfile200.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents200 = h_nevents.GetBinContent(2)

  tfile140 = TFile.Open(infile140)
  h_nevents = tfile140.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents140 = h_nevents.GetBinContent(2)

  tfile250 = TFile.Open(infile250)
  h_nevents = tfile250.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents250 = h_nevents.GetBinContent(2)

  tfile300 = TFile.Open(infile300)
  h_nevents = tfile300.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents300 = h_nevents.GetBinContent(2)


  # ____________________________________________________________________________
  # emtf_ptmin20_qmin12_pu

  #hname = "highest_emtf_absEtaMin1.24_absEtaMax2.4_qmin12_pt"
  hname = "highest_emtf_absEtaMin0.8_absEtaMax1.24_qmin12_pt"

  pileup_list = [140, 200, 250, 300]
  tfile_list = [tfile140, tfile200, tfile250, tfile300]
  nevents_list = [nevents140, nevents200, nevents250, nevents300]

  rates = []

  for pileup, tfile, nevents in zip(pileup_list, tfile_list, nevents_list):
    h = tfile.Get(hname)
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    emtf_rate = h.GetBinContent(h.FindBin(20))
    emtf_rate_err = h.GetBinError(h.FindBin(20))
    #
    h = tfile.Get(hname2026_f(hname))
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    emtf2026_rate = h.GetBinContent(h.FindBin(20))
    emtf2026_rate_err = h.GetBinError(h.FindBin(20))
    #
    rates.append((float(pileup), emtf_rate, emtf_rate_err, emtf2026_rate, emtf2026_rate_err))

  print rates

  gr_denom = TGraphAsymmErrors(len(rates))
  gr_numer = TGraphAsymmErrors(len(rates))

  for i, x in enumerate(rates):
    pu, emtf_rate, emtf_rate_err, emtf2026_rate, emtf2026_rate_err = x
    gr_denom.SetPoint(i, pu, emtf_rate)
    gr_denom.SetPointError(i, 0, 0, emtf_rate_err, emtf_rate_err)
    gr_numer.SetPoint(i, pu, emtf2026_rate)
    gr_numer.SetPointError(i, 0, 0, emtf2026_rate_err, emtf2026_rate_err)

  gr_denom.SetMarkerStyle(20)
  gr_denom.SetMarkerSize(1.4)
  gr_denom.SetMarkerColor(632)  # kRed
  gr_denom.SetLineWidth(2)
  gr_denom.SetLineColor(632)  # kRed

  gr_numer.SetMarkerStyle(20)
  gr_numer.SetMarkerSize(1.4)
  gr_numer.SetMarkerColor(600)  # kBlue
  gr_numer.SetLineWidth(2)
  gr_numer.SetLineColor(600)  # kBlue

  frame = TH1F("frame", "; PU; Trigger rate [kHz]", 100, 0, 350)
  frame.SetMinimum(0)
  frame.SetMaximum(25)
  frame.Draw()

  #gr_denom.Draw("P")
  gr_numer.Draw("P")

  for i, x in enumerate(rates):
    pu, emtf_rate, emtf_rate_err, emtf2026_rate, emtf2026_rate_err = x
    #tlatex.DrawLatex(pu - 3, emtf_rate + 4, "%.1f" % emtf_rate)
    tlatex.DrawLatex(pu - 3, emtf2026_rate + 2, "%.1f" % emtf2026_rate)

  def draw_cms_lumi_1():
    tlatexCMS1.SetTextSize(0.75*0.05)
    tlatexCMS2.SetTextSize(0.60*0.05)
    tlatexCMS3.SetTextSize(0.60*0.05)
    tlatexCMS1.DrawLatex(0.164, 0.965, 'CMS')
    tlatexCMS2.DrawLatex(0.252, 0.965, 'Phase-2 Simulation')
    tlatexCMS3.DrawLatex(0.675, 0.965, '<PU>=140, 200, 250, 300')

  draw_cms_lumi_1()
  imgname = "emtf_ptmin20_qmin12_pu"
  gPad.Print("figures_perf/" + imgname + "_omtf" + ".png")
  gPad.Print("figures_perf/" + imgname + "_omtf" + ".pdf")
  donotdelete.append([frame, gr_denom, gr_numer])
