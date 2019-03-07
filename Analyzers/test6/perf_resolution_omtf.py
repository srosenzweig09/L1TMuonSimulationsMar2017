#!/usr/bin/env python

hname2026_f = lambda hname: "emtf2026_" + hname[5:]

donotdelete = []

infile = "histos_tbc_omtf_add.25.root"


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
  #gStyle.SetPadRightMargin(0.05)
  #gStyle.SetTitleOffset(0.9, "Y")
  #gStyle.SetTitleOffset(1.2, "Y")
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
    tlatexCMS3.DrawLatex(0.885, 0.965, '<PU>=0')

  # Open file
  tfile = TFile.Open(infile)


  # ____________________________________________________________________________
  # emtf_l1pt_vs_genpt
  hname = "emtf_l1pt_vs_genpt"
  h2a = tfile.Get(hname)
  h2a.Draw("COLZ")
  draw_cms_lumi()
  gPad.Print("figures_perf/" + hname + "_omtf" + ".png")
  gPad.Print("figures_perf/" + hname + "_omtf" + ".pdf")

  hname = hname2026_f(hname)
  h2b = tfile.Get(hname)
  h2b.Draw("COLZ")
  draw_cms_lumi()
  gPad.Print("figures_perf/" + hname + "_omtf" + ".png")
  gPad.Print("figures_perf/" + hname + "_omtf" + ".pdf")

  donotdelete.append([h2a, h2b])

  # ____________________________________________________________________________
  # emtf_l1ptres_vs_genpt
  def doit():
    frame = h.ProfileX(hname+"_frame", 1, -1, "s")
    gr1 = TGraphAsymmErrors(h.GetNbinsX())
    gr2 = TGraphAsymmErrors(h.GetNbinsX())
    gr1_aspt = TGraphAsymmErrors(h.GetNbinsX())
    gr2_aspt = TGraphAsymmErrors(h.GetNbinsX())
    # Apply gaussian fits
    for i in xrange(h.GetNbinsX()):
      h_py = h.ProjectionY("_py", i+1, i+1)

      if 50 <= i <= 60:  # high pT, not enough entries (300 bins -> 150)
        h_py.Rebin(2)
      elif i >= 78:      # low pT, resolution affected by finite bin width
        h_py = h.ProjectionY("_py", i+1, i+2)  # merge i & (i+1) entries
        if i == 82:      # even lower pT, resolution affected by finite bin width
          h_py = h.ProjectionY("_py", i+1, i+4)  # merge i & (i+4) entries
        elif i >= 82:
          continue

      if h_py.Integral() < 20:  continue
      r = h_py.Fit("gaus", "SNQ", "", -1, 1.2)
      #r = h_py.Fit("gaus", "SNQ", "", h_py.GetMean() - 0.04*5, h_py.GetMean() + 0.04*5)
      mean, sigma, meanErr, sigmaErr = r.Parameter(1), r.Parameter(2), r.ParError(1), r.ParError(2)
      gr1.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), mean)
      gr1.SetPointError(i, 0, 0, sigma, sigma)
      gr2.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), sigma)
      gr2.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
      gr1_aspt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), mean)
      gr1_aspt.SetPointError(i, 0, 0, sigma, sigma)
      gr2_aspt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), sigma)
      gr2_aspt.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
    # Draw
    h.Draw("COLZ")
    gPad.SetLogx(0)
    #draw_cms_lumi()
    #gPad.Print("figures_perf/" + hname + "_omtf" + ".png")
    #gPad.Print("figures_perf/" + hname + "_omtf" + ".pdf")
    #
    frame.Reset()
    frame.SetBins(50, 0, 50)
    frame.GetXaxis().SetTitle("gen p_{T} [GeV]")
    frame.GetYaxis().SetTitle("#Delta(p_{T})/p_{T} bias")
    frame.SetMaximum(0.5)
    frame.SetMinimum(-0.5)
    frame.SetStats(0)
    frame.Draw()
    gr1_aspt.SetLineColor(col)
    gr1_aspt.SetMarkerColor(col)
    gr1_aspt.Draw("p")
    gPad.SetLogx()
    draw_cms_lumi()
    gPad.Print("figures_perf/" + hname + "_bias" + "_omtf" + ".png")
    gPad.Print("figures_perf/" + hname + "_bias" + "_omtf" + ".pdf")
    #
    frame.GetXaxis().SetTitle("gen p_{T} [GeV]")
    frame.GetYaxis().SetTitle("#Delta(p_{T})/p_{T} resolution")
    frame.SetMaximum(0.6)
    frame.SetMinimum(0.0)
    frame.SetStats(0)
    frame.Draw()
    gr2_aspt.SetLineColor(col)
    gr2_aspt.SetMarkerColor(col)
    gr2_aspt.Draw("p")
    #gr2_aspt.Fit("pol1", "", "", 10, 40)
    gPad.SetLogx()
    draw_cms_lumi()
    gPad.Print("figures_perf/" + hname + "_res" + "_omtf" + ".png")
    gPad.Print("figures_perf/" + hname + "_res" + "_omtf" + ".pdf")
    #
    h.cache = [frame, gr1, gr2, gr1_aspt, gr2_aspt]

  hname = "emtf_l1ptres_vs_genpt"
  h = tfile.Get(hname)
  col = 632  # kRed
  doit()

  hname = hname2026_f(hname)
  h = tfile.Get(hname)
  col = 600  # kBlue
  doit()

  donotdelete.append([h2a, h2b])
