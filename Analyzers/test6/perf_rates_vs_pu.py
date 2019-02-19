#!/usr/bin/env python

hname2026_f = lambda hname: "highest_emtf2026_" + hname[13:]

donotdelete = []

infile = "histos_tbb_add.24.root"
infile140 = "histos_tbb_140_add.24.root"
infile250 = "histos_tbb_250_add.24.root"
infile300 = "histos_tbb_300_add.24.root"

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

  hname = "highest_emtf_absEtaMin1.24_absEtaMax2.4_qmin12_pt"
  #hname = "highest_emtf_absEtaMin0.8_absEtaMax1.24_qmin12_pt"

  rates = []

  for pileup, tfile, nevents in zip([140, 200, 250, 300], [tfile140, tfile200, tfile250, tfile300], [nevents140, nevents200, nevents250, nevents300]):
    h = tfile.Get(hname)
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    emtf_rate = h.GetBinContent(h.FindBin(20))
    #
    h = tfile.Get(hname2026_f(hname))
    h.Sumw2()
    make_ptcut(h)
    make_rate(h, nevents)
    emtf2026_rate = h.GetBinContent(h.FindBin(20))
    #
    rates.append((float(pileup), emtf_rate, emtf2026_rate))

  print rates

  gr_denom = TGraph(len(rates)+1)
  gr_denom.SetPoint(0, 0, 0)
  for i, rate in enumerate(rates):
    gr_denom.SetPoint(i+1, rate[0], rate[1])

  gr_denom_lin = TGraph(3)
  gr_denom_lin.SetPoint(0, 0, 0)
  ref_rate = rates[1]  # PU200
  gr_denom_lin.SetPoint(1, ref_rate[0], ref_rate[1])
  gr_denom_lin.SetPoint(2, 400, 400/ref_rate[0]*ref_rate[1])

  gr_numer = TGraph(len(rates)+1)
  gr_numer.SetPoint(0, 0, 0)
  for i, rate in enumerate(rates):
    gr_numer.SetPoint(i+1, rate[0], rate[2])

  gr_numer_lin = TGraph(3)
  gr_numer_lin.SetPoint(0, 0, 0)
  ref_rate = rates[1]  # PU200
  gr_numer_lin.SetPoint(1, ref_rate[0], ref_rate[2])
  gr_numer_lin.SetPoint(2, 400, 400/ref_rate[0]*ref_rate[2])

  gr_denom.SetMarkerStyle(20)
  gr_denom.SetMarkerSize(1.4)
  gr_denom.SetMarkerColor(632)  # kRed
  gr_denom_lin.SetLineStyle(2)
  gr_denom_lin.SetLineWidth(2)
  gr_denom_lin.SetLineColor(632)  # kRed

  gr_numer.SetMarkerStyle(20)
  gr_numer.SetMarkerSize(1.4)
  gr_numer.SetMarkerColor(600)  # kBlue
  gr_numer_lin.SetLineStyle(2)
  gr_numer_lin.SetLineWidth(2)
  gr_numer_lin.SetLineColor(600)  # kBlue

  frame = TH1F("frame", "; PU; Trigger rate [kHz]", 100, 0, 350)
  frame.SetMinimum(0)
  frame.SetMaximum(88)
  frame.Draw()

  gr_denom.Draw("P")
  gr_numer.Draw("P")
  gr_denom_lin.Draw("C")
  gr_numer_lin.Draw("C")

  def draw_cms_lumi_1():
    tlatexCMS1.SetTextSize(0.75*0.05)
    tlatexCMS2.SetTextSize(0.60*0.05)
    tlatexCMS3.SetTextSize(0.60*0.05)
    tlatexCMS1.DrawLatex(0.164, 0.965, 'CMS')
    tlatexCMS2.DrawLatex(0.252, 0.965, 'Phase-2 Simulation')
    tlatexCMS3.DrawLatex(0.675, 0.965, '<PU>=140, 200, 250, 300')

  draw_cms_lumi_1()
  imgname = "emtf_ptmin20_qmin12_pu"
  gPad.Print("figures_perf/" + imgname + ".png")
  gPad.Print("figures_perf/" + imgname + ".pdf")
  donotdelete.append([frame, gr_denom, gr_numer, gr_denom_lin, gr_numer_lin])
