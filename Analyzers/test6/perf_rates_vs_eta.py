#!/usr/bin/env python

hname2026_f = lambda hname: "emtf2026_" + hname[5:]

donotdelete = []

infile = "histos_tbb_add.25.root"

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

  def draw_cms_lumi():
    tlatexCMS1.DrawLatex(0.164, 0.965, 'CMS')
    tlatexCMS2.DrawLatex(0.252, 0.965, 'Phase-2 Simulation')
    tlatexCMS3.DrawLatex(0.865, 0.965, '<PU>=200')

  tfile = TFile.Open(infile)
  h_nevents = tfile.Get("nevents")
  assert h_nevents != None, "Cannot get nevents"
  nevents = h_nevents.GetBinContent(2)


  # ____________________________________________________________________________
  # emtf_ptmin20_qmin12_eta

  hname = "emtf_ptmin20_qmin12_eta"
  h1a = tfile.Get(hname)
  h1b = tfile.Get(hname2026_f(hname))

  make_rate(h1a, nevents)
  make_rate(h1b, nevents)

  h1a.SetMarkerColor(632)  # kRed
  h1a.SetLineColor(632)  # kRed
  h1a.SetLineWidth(2)

  h1b.SetMarkerColor(600)  # kBlue
  h1b.SetLineColor(600)  # kBlue
  h1b.SetFillColor(600)  # kBlue
  h1b.SetLineWidth(2)

  h1b_clone = h1b.Clone(h1b.GetName() + "_clone")
  h1b_clone.SetLineColor(1)  # kBlack

  h1a.GetXaxis().SetTitle("L1 muon #eta")
  h1a.GetYaxis().SetTitle("Trigger rate per unit #eta [kHz]")

  h1a.Draw("hist")
  h1b.Draw("same hist")
  h1b_clone.Draw("same e")

  draw_cms_lumi()
  gPad.Print("figures_perf/" + hname + ".png")
  gPad.Print("figures_perf/" + hname + ".pdf")
  donotdelete.append([h1a, h1b, h1b_clone])
